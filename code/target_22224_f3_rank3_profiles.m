//////////////////////////////////////////////////////////////////////
// Local Mordell--Weil coefficient profiles for the rank-three (1,2)
// elliptic quotient of F3 = (-612,34,289,-338*t^2).
//
// At p=13 the corrected primitive infinity-chart calculation forces
// t to have denominator divisible by 13^2.  The finite special-fibre
// mask is empty, so this profile retains only points mapping to t=infinity;
// exact candidates are subsequently checked for valuation >= 2.
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(3);
Q:=Rationals(); Z:=Integers(); R<T>:=PolynomialRing(Q);
if not assigned PrimeList then PrimeList:=[11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97];
elif Type(PrimeList) eq MonStgElt then PrimeList:=[StringToInteger(s):s in Split(PrimeList,",")]; end if;
if not assigned output_file then output_file:="results/target_22224_f3_rank3_profiles.tsv"; end if;
if not assigned log_file then log_file:="results/target_22224_f3_rank3_profiles.log"; end if;
SetLogFile(log_file:Overwrite:=true);

fixed:=[Q!-612,Q!34,Q!289]; d0:=Q!-338; pair:=[1,2]; other:=3;
Rx:=(fixed[pair[1]]+d0*T^2)/(fixed[pair[1]]+d0);
Ry:=(fixed[pair[2]]+d0*T^2)/(fixed[pair[2]]+d0);
C:=HyperellipticCurve(Rx*Ry); P0:=C![Q!1,Q!1,Q!1];
Eraw,phi:=EllipticCurve(C,P0); Einv:=Inverse(phi);
E,minmap:=MinimalModel(Eraw); minmapinv:=Inverse(minmap);
gens:=Generators(E); free:=[g:g in gens|Order(g) eq 0]; assert #free eq 3;
TG,tmap:=TorsionSubgroup(E); tors:=[tmap(g):g in TG];
assert Invariants(TG) eq [2];
rawMap:=DefiningPolynomials(minmapinv); curveMap:=DefiningPolynomials(Einv);

function ReduceRat(q,F)
  q:=Q!q; ell:=Z!Characteristic(F);
  if Denominator(q) mod ell eq 0 then return false,F!0; end if;
  return true,F!Numerator(q)/F!Denominator(q);
end function;
function ReducePoint(PQ,EF)
  if PQ eq Parent(PQ)!0 then return true,EF!0; end if;
  cs:=[Q!PQ[i]:i in [1..3]]; den:=LCM([Denominator(z):z in cs]);
  ints:=[Z!(z*den):z in cs]; g:=GCD([Abs(z):z in ints]);
  if g ne 0 then ints:=[z div g:z in ints]; end if;
  try return true,EF![BaseRing(EF)!z:z in ints];
  catch e return false,EF!0; end try;
end function;
function ReduceMap(polys,F)
  out:=[];
  try
    for h in polys do Append(~out,ChangeRing(h,F)); end for;
    return true,out;
  catch e
    allcoeffs:=&cat[Coefficients(h):h in polys];
    den:=LCM([Denominator(Q!c):c in allcoeffs]); scaled:=[den*h:h in polys];
    coeffs:=[Z!c:c in &cat[Coefficients(h):h in scaled]];
    cont:=GCD([Abs(c):c in coeffs]); if cont eq 0 then return false,[]; end if;
    scaled:=[h/cont:h in scaled];
    out:=[];
    try for h in scaled do Append(~out,ChangeRing(h,F)); end for;
    catch e2 return false,[]; end try;
  end try;
  return true,out;
end function;
function EvalMap(polys,coords)
  vals:=[Evaluate(h,coords):h in polys];
  if &and[x eq 0:x in vals] then return false,vals; end if;
  return true,vals;
end function;
function TFromPoint(PF,rawF,curveF,F)
  xyz:=[PF[i]:i in [1..3]]; ok,raw:=EvalMap(rawF,xyz);
  if not ok then return true,true,F!0; end if;
  ok,cp:=EvalMap(curveF,raw);
  if not ok or cp[3] eq 0 then return true,true,F!0; end if;
  return true,false,cp[1]/cp[3];
end function;
function ContactKey(vals,F)
  ell:=Z!Characteristic(F);
  if &and[F!z eq 0:z in vals] then return false,0; end if;
  best:=ell^4;
  for u in F do if u eq 0 then continue; end if;
    sq:=Sort([Z!((u*(F!z))^2):z in vals]);
    code:=sq[1]*ell^3+sq[2]*ell^2+sq[3]*ell+sq[4];
    if code lt best then best:=code; end if;
  end for;
  return true,best;
end function;
function DirectContact(vals,F)
  ell:=Z!Characteristic(F); ok,key:=ContactKey(vals,F);
  if not ok then return false; end if;
  if ell eq 11 then
    return key in {1,12,133,136,137,159,1464,1488,1489,1490,1505,1512,1516};
  elif ell eq 13 then
    return key in {1,14,183,185,186,191,192,194,220,2380,2408,2414,2422,2428,2534};
  end if;
  return true;
end function;
// 0=fail; 2=open and 3 divides #J; 3=boundary.
function StateFinite(tt,F)
  ell:=Z!Characteristic(F); P<X>:=PolynomialRing(F);
  vals:=[F!fixed[1],F!fixed[2],F!fixed[3],F!d0*tt^2];
  if ell in {11,13} and not DirectContact(vals,F) then return 0,-1; end if;
  den:=F!(fixed[other]+d0);
  if den eq 0 then return 3,-1; end if;
  r3:=(F!fixed[other]+F!d0*tt^2)/den;
  if not IsSquare(r3) then return 0,-1; end if;
  f:=X*&*[X+z^2:z in vals];
  if Degree(f) ne 5 or Discriminant(f) eq 0 or r3 eq 0 then return 3,-1; end if;
  nj:=Z!Evaluate(LPolynomial(HyperellipticCurve(f)),1);
  return (nj mod 3 eq 0) select 2 else 0,nj;
end function;

print "F3_RANK3_PROFILES_START","E",E,"rank_bounds",RankBounds(E),"gens",gens;
out:=Open(output_file,"w");
fprintf out,"prime\ttorsion_coset\tm\tn\tk\to1\to2\to3\tstate\tt\n";
for ell in PrimeList do
  F:=GF(ell); ainv:=[]; good:=true;
  for q in aInvariants(E) do ok,z:=ReduceRat(q,F); if not ok then good:=false;break;end if;Append(~ainv,z);end for;
  if not good then print "F3_PROFILE_SKIP",ell,"E_den";continue;end if;
  try EF:=EllipticCurve(ainv);catch e print "F3_PROFILE_SKIP",ell,"E_bad";continue;end try;
  gf:=[];
  for G in free do ok,g:=ReducePoint(G,EF);if not ok then good:=false;break;end if;Append(~gf,g);end for;
  if not good then print "F3_PROFILE_SKIP",ell,"gen";continue;end if;
  tf:=[];
  for t in tors do ok,g:=ReducePoint(t,EF);if not ok then good:=false;break;end if;Append(~tf,g);end for;
  if not good then print "F3_PROFILE_SKIP",ell,"tors";continue;end if;
  okr,rawF:=ReduceMap(rawMap,F);okc,curveF:=ReduceMap(curveMap,F);
  if not okr or not okc then print "F3_PROFILE_SKIP",ell,"map",okr,okc;continue;end if;
  oo:=[Order(g):g in gf]; total:=#tf*&*oo; allowed:=0; infs:=0; opens:=0; bounds:=0;
  cache:=[-1:i in [1..ell]];
  for z in [0..ell-1] do st,nj:=StateFinite(F!z,F);cache[z+1]:=st;end for;
  for ti in [1..#tf] do
    for m in [0..oo[1]-1] do for n in [0..oo[2]-1] do for k in [0..oo[3]-1] do
      ep:=m*gf[1]+n*gf[2]+k*gf[3]+tf[ti];
      ok,isinf,tt:=TFromPoint(ep,rawF,curveF,F);
      if not ok or isinf then st:=3; ts:="inf";infs+:=1;
      else st:=cache[Z!tt+1];ts:=Sprint(Z!tt);end if;
      if st in {2,3} then
        allowed+:=1;if st eq 2 then opens+:=1;else bounds+:=1;end if;
        fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",
          ell,ti,m,n,k,oo[1],oo[2],oo[3],st,ts;
      end if;
    end for;end for;end for;
  end for;
  print "F3_PROFILE","p",ell,"orders",oo,"total",total,"map_inf",infs,
        "open3",opens,"boundary",bounds,"allowed",allowed;
end for;
delete out; print "F3_RANK3_PROFILES_DONE",output_file;
UnsetLogFile();quit;
