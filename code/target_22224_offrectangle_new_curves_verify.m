//////////////////////////////////////////////////////////////////////
// Exact torsion and geometric-simplicity verification for the lowest-height
// new pair-scaling-fiber points.
//////////////////////////////////////////////////////////////////////

SetColumns(0); SetSeed(1);
Q:=Rationals(); Z:=Integers(); P<x>:=PolynomialRing(Q);
if not assigned input_file then
    input_file:="results/target_22224_offrectangle_new_curves_pairfiber_small.tsv";
end if;
if not assigned HeightCap then HeightCap:=10^18; end if;
if Type(HeightCap) eq MonStgElt then HeightCap:=StringToInteger(HeightCap); end if;
if not assigned log_file then
    log_file:="results/target_22224_offrectangle_new_curves_verify.log";
end if;
SetLogFile(log_file:Overwrite:=true);

function StrictSimpleWitness(f)
    for ell in PrimesUpTo(199) do
        try
            fp:=ChangeRing(f,GF(ell));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
            lp:=LPolynomial(HyperellipticCurve(fp)); R<t>:=PolynomialRing(Q); chi:=R!lp;
            if not IsIrreducible(chi) then continue; end if;
            K<pi>:=NumberField(chi);
            if &and[Degree(MinimalPolynomial(pi^n)) eq 4:n in [2..12]] then
                return true,ell,chi;
            end if;
        catch e continue;
        end try;
    end for;
    R<t>:=PolynomialRing(Q); return false,0,R!0;
end function;

rows:=[];
for line in Split(Read(input_file),"\n") do
    s:=Split(line,"\t");
    if #s lt 14 or s[1] eq "fiber" then continue; end if;
    known:=StringToInteger(s[12]);
    vals:=[StringToInteger(s[i]):i in [8..11]];
    raw:=[Q!StringToRational(s[i]):i in [4..7]];
    if known eq 0 and Maximum(vals) le HeightCap then Append(~rows,<s[1],vals,raw>); end if;
end for;

print "OFFRECTANGLE_NEW_CURVES_VERIFY_START","rows",#rows,"HeightCap",HeightCap;
for row in rows do
    name:=row[1]; vals:=row[2]; raw:=row[3]; a,b,c,d:=Explode(raw);
    rad:=[a*b*c*d,a*(a+b)*(a+c)*(a+d),b*(a+b)*(b+c)*(b+d),c*(a+c)*(b+c)*(c+d)];
    assert &and[IsSquare(Q!z):z in rad];
    f:=x*&*[x+(Q!z)^2:z in vals]; C:=HyperellipticCurve(f);
    G,mp:=TorsionSubgroup(Jacobian(C)); inv:=Invariants(G);
    simple,sp,chi:=StrictSimpleWitness(f);
    print "EXACT_NEW_POINT",name,"tuple",vals,"torsion",inv,
          "positive_signed_presentation",&and[z gt 0:z in raw],
          "simple",simple,"certificate_prime",sp,"chi",chi;
end for;
print "OFFRECTANGLE_NEW_CURVES_VERIFY_DONE","rows",#rows;
UnsetLogFile(); quit;
