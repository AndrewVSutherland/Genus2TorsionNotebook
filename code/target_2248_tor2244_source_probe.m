// Reverse-source probe: Richelot neighbors of the A(2,2,4,4) tuple bank.
SetColumns(0); SetSeed(1);
if not assigned start_row then start_row:=1; elif Type(start_row) eq MonStgElt then start_row:=StringToInteger(start_row); end if;
if not assigned max_rows then max_rows:=25; elif Type(max_rows) eq MonStgElt then max_rows:=StringToInteger(max_rows); end if;
if not assigned bank_file then bank_file:="data/tor2244_bank.txt"; end if;
if not assigned write_log then write_log:=true;
elif Type(write_log) eq MonStgElt then write_log:=write_log in {"true","True","1","yes"}; end if;
if not assigned log_file then log_file:="results/target_2248_tor2244_source_probe.log"; end if;
if not assigned verbose_candidates then verbose_candidates:=false;
elif Type(verbose_candidates) eq MonStgElt then verbose_candidates:=verbose_candidates in {"true","True","1","yes"}; end if;
if write_log then SetLogFile(log_file : Overwrite:=true); end if;
Q:=Rationals(); Z:=Integers(); P<x>:=PolynomialRing(Q);

function ReadTupleFile(filename)
    rows:=Split(Read(filename),"\n"); out:=[];
    for raw in rows do
        n:=#raw; if n gt 0 and raw[n] eq "\r" then n-:=1; end if;
        if n ge 2 and raw[1] eq "[" and raw[n] eq "]" then
            ss:=Split(raw[2..n-1],",");
            if #ss eq 4 then Append(~out,[StringToInteger(s):s in ss]); end if;
        end if;
    end for; return out;
end function;
function Scale(f) d:=LCM([Denominator(Coefficient(f,i)):i in [0..Degree(f)]]); return P!(d^2*f); end function;
function Norm(C)
 f,h:=HyperellipticPolynomials(C); F:=h eq 0 select P!f else P!(h^2+4*f);
 D:=SimplifiedModel(HyperellipticCurve(Scale(F))); fD,hD:=HyperellipticPolynomials(D); return D,P!fD,P!hD;
end function;
function CF(C) f,h:=HyperellipticPolynomials(C); return h eq 0 select P!f else P!(h^2+4*f); end function;
function FD(C) return Sort(&cat[[Degree(z[1]):j in [1..z[2]]]:z in Factorization(CF(C))]); end function;
function One(C) F:=CF(C); ds:=FD(C); return (Degree(F) eq 5 and ds eq [1,1,1,2]) or (Degree(F) eq 6 and ds eq [1,1,1,1,2]); end function;
function Gate(C,v)
 F:=CF(C); used:=0; g:=0; rows:=[];
 for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47] do
  try fp:=ChangeRing(F,GF(p)); if Degree(fp) ne Degree(F) or Discriminant(fp) eq 0 then continue; end if;
   L:=LPolynomial(ChangeRing(C,GF(p))); n:=Z!Evaluate(L,1); g:=used eq 0 select n else GCD(g,n);
   used+:=1; Append(~rows,<p,n,Valuation(n,2)>); if Valuation(g,2) lt v then return false,g,rows; end if;
   if used ge 2 then return true,g,rows; end if;
  catch e continue; end try;
 end for; return false,g,rows;
end function;

rows:=ReadTupleFile(bank_file);
last:=Min(#rows,start_row+max_rows-1); edges:=0; one:=0; pass:=0; exact:=0; hits:=0;
print "TARGET_2248_TOR2244_SOURCE_PROBE","bank_file",bank_file,
      "start",start_row,"last",last;
for ri in [start_row..last] do
 t:=rows[ri]; a,b,c,d:=Explode(t); f:=x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2); C,fC,hC:=Norm(HyperellipticCurve(f));
 Rs:=RichelotIsogenousSurfaces(Jacobian(C));
 for j in [1..#Rs] do
  edges+:=1; if Type(Rs[j]) ne JacHyp then continue; end if; D,fD,hD:=Norm(Curve(Rs[j]));
  if not One(D) then continue; end if; one+:=1; ok,g,rr:=Gate(D,6); if not ok then continue; end if; pass+:=1;
  exact+:=1; G,mp:=TorsionSubgroup(Jacobian(D)); inv:=Invariants(G);
  print "CANDIDATE","row",ri,"tuple",t,"edge",j,"fd",FD(D),
        "gcd",g,"red",rr,"torsion",inv,"order",#G;
  if verbose_candidates then print "CANDIDATE_CURVE","f",fD,"h",hD; end if;
  if inv eq [2,4,8] then hits+:=1; print "SOURCE_248_HIT","row",ri,"edge",j,"f",fD,"h",hD; end if;
 end for;
end for;
print "SUMMARY","rows",last-start_row+1,"edges",edges,"one_split",one,"reduction_pass",pass,"exact",exact,"hits",hits;
quit;
