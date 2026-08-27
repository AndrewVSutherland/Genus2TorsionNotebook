//////////////////////////////////////////////////////////////////////
// Exact verifier for CRT-reconstructed direct-contact candidates.
//////////////////////////////////////////////////////////////////////
SetColumns(0);SetSeed(1);
if not assigned candidate_file then
 candidate_file:="results/target_22224_direct_contact_deep13_crt_seed22224.tsv";
end if;
if not assigned log_file then
 log_file:="results/target_22224_direct_contact_deep13_verify.log";
end if;
SetLogFile(log_file:Overwrite:=true);
Q:=Rationals();Z:=Integers();P<X>:=PolynomialRing(Q);

function QR(n,d)return Q!StringToInteger(n)/StringToInteger(d);end function;
function DirectData(u,t,w)
 q:=X^2+u*X+t^2;h:=X^3+(1+3*u)/2*X^2+w*X+t^3;
 F:=ExactQuotient(h^2-q^3,X);return F,q,h;
end function;
function SquareQ(r)
 if r lt 0 then return false,Q!0;end if;
 a,sa:=IsSquare(Z!Numerator(r));b,sb:=IsSquare(Z!Denominator(r));
 return a and b,(a and b) select Q!sa/sb else Q!0;
end function;
function Branches(F)
 fac:=Factorization(F);if #fac ne 4 or not &and[Degree(z[1]) eq 1 and z[2] eq 1:z in fac] then return false,[];end if;
 ans:=[];for z in fac do rr:=-Coefficient(z[1],0)/Coefficient(z[1],1);ok,a:=SquareQ(-rr);if not ok or a eq 0 then return false,[];end if;Append(~ans,a);end for;
 return #Set(ans) eq 4,ans;
end function;
function Rads(a,b,c,d)
 return [a*b*c*d,a*(a+b)*(a+c)*(a+d),b*(b+a)*(b+c)*(b+d),c*(c+a)*(c+b)*(c+d)];
end function;
function IntegralModel(f)
 den:=LCM([Denominator(Coefficient(f,i)):i in [0..Degree(f)]]);return P!(den^2*f),den;
end function;
function SimpleWitness(f)
 C:=HyperellipticCurve(f);R<T>:=PolynomialRing(Q);
 for p in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97] do
  try fp:=ChangeRing(f,GF(p));if Discriminant(fp) eq 0 then continue;end if;chi:=R!LPolynomial(HyperellipticCurve(fp));if not IsIrreducible(chi) then continue;end if;K<pi>:=NumberField(chi);if &and[Degree(MinimalPolynomial(pi^n)) eq 4:n in [2..12]] then return true,p,chi;end if;catch e;end try;
 end for;return false,0,R!0;
end function;

lines:=Split(Read(candidate_file),"\n");rows:=0;split:=0;cover:=0;jtests:=0;hits:=0;seen:={};
print "TARGET_22224_DIRECT_CONTACT_VERIFY_START",candidate_file,"lines",#lines;
for i in [2..#lines] do if #lines[i] eq 0 then continue;end if;z:=Split(lines[i],"\t");if #z lt 6 then continue;end if;
 rows+:=1;u:=QR(z[1],z[2]);t:=QR(z[3],z[4]);w:=QR(z[5],z[6]);if t eq 0 then continue;end if;
 F,q,h:=DirectData(u,t,w);if h^2-q^3 ne X*F then error "identity";end if;
 ok,mags:=Branches(F);if not ok then continue;end if;split+:=1;f:=X*F;if Discriminant(f) eq 0 then continue;end if;
 for mask in [0..15] do
  vals:=[((mask div 2^(j-1)) mod 2 eq 0) select mags[j] else -mags[j]:j in [1..4]];
  a0,b0,c0,d0:=Explode(vals);arr:=[<a0,b0,c0,d0>,<a0,c0,b0,d0>,<a0,d0,b0,c0>];
  for pairing in [1..3] do ar:=arr[pairing];a:=ar[1];b:=ar[2];c:=ar[3];d:=ar[4];
   rr:=Rads(a,b,c,d);oks:=true;for r in rr do sq,sr:=SquareQ(r);if not sq or r eq 0 then oks:=false;break;end if;end for;if not oks then continue;end if;
   cover+:=1;key:=Sprint(<f,a,b,c,d,q,h>);if key in seen then continue;end if;Include(~seen,key);
   fI,sc:=IntegralModel(f);J:=Jacobian(HyperellipticCurve(fI));D3:=J![q,sc*(h mod q)];
   g4:=(X-a*b)*(X-c*d)-X*(a+b)*(c+d);L4:=(X-a*b)*(c+d)+(a+b)*(X-c*d);
   D4:=J![g4/LeadingCoefficient(g4),sc*((-X*L4) mod g4)];jtests+:=1;
   if Order(D3) ne 3 or Order(D4) ne 4 then continue;end if;
   divisible,H:=IsDivisibleBy(D4,2);if not divisible then continue;end if;
   T,mp:=TorsionSubgroup(J);inv:=Invariants(T);simple,pcert,chi:=SimpleWitness(fI);hits+:=1;
   print "TARGET_22224_DIRECT_CONTACT_HIT","u",u,"t",t,"w",w,"pairing",pairing,
         "branches",[a,b,c,d],"torsion",inv,"half_order",Order(H),
         "simple",simple,"pcert",pcert,"chi",chi,"curve",fI;
  end for;
 end for;
end for;
// Positive control: the order-96 record is on the direct contact chart and
// must reach the exact split stage, although its marked D4 is not divisible.
u0:=Q!289578289/21999628800;t0:=Q!289578289/43999257600;
w0:=Q!11494656094062061921/645311556450385920000;
F0,q0,h0:=DirectData(u0,t0,w0);ok0,m0:=Branches(F0);
print "RECORD_CONTROL","identity",h0^2-q0^3 eq X*F0,"split",ok0,"branches",m0;
print "TARGET_22224_DIRECT_CONTACT_VERIFY_DONE","rows",rows,"split",split,
      "cover_presentations",cover,"jacobian_tests",jtests,"hits",hits;
UnsetLogFile();quit;
