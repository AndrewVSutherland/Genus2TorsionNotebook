//////////////////////////////////////////////////////////////////////
// Exact verifier for local-contact survivors from the transverse box.
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(1);
if not assigned input_file then
    input_file := "results/target_22224_a2228_deep_p5_box_B1000_allsigns.tsv";
end if;
if not assigned log_file then
    log_file := "results/target_22224_a2228_deep_p5_box_B1000_verify.log";
end if;
SetLogFile(log_file : Overwrite:=true);
Q:=Rationals(); P<x>:=PolynomialRing(Q);

function SimpleWitness(f)
    for p in PrimesUpTo(199) do
        try
            fp:=ChangeRing(f,GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
            R<T>:=PolynomialRing(Q); chi:=R!LPolynomial(HyperellipticCurve(fp));
            if not IsIrreducible(chi) then continue; end if;
            K<pi>:=NumberField(chi);
            if &and[Degree(MinimalPolynomial(pi^n)) eq 4:n in [2..12]] then
                return true,p,chi;
            end if;
        catch e
            continue;
        end try;
    end for;
    R<T>:=PolynomialRing(Q); return false,0,R!0;
end function;

rows:=0; hits:=0;
print "A2228_DEEP_P5_BOX_VERIFY_START",input_file;
for line in Split(Read(input_file),"\n") do
    if #line eq 0 or line[1] notin {"1","2","3","4","5","6","7","8","9"} then continue; end if;
    s:=Split(line,"\t");
    if #s lt 14 or s[12] ne "1" then continue; end if;
    vals:=[StringToInteger(s[i]):i in [6..9]];
    rows+:=1;
    f:=x*&*[x+(Q!z)^2:z in vals];
    G,mp:=TorsionSubgroup(Jacobian(HyperellipticCurve(f)));
    inv:=Invariants(G); simple,sp,chi:=SimpleWitness(f);
    print "EXACT_LOCAL_SURVIVOR","tuple",vals,"torsion",inv,
          "order",#G,"simple",simple,"certificate_prime",sp,"chi",chi;
    if #G mod 3 eq 0 then
        hits+:=1;
        print "TARGET_22224_HIT","tuple",vals,"torsion",inv,"simple",simple;
    end if;
end for;
print "A2228_DEEP_P5_BOX_VERIFY_DONE","rows",rows,"hits",hits;
UnsetLogFile(); quit;
