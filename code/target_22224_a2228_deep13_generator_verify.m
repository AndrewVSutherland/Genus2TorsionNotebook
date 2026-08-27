//////////////////////////////////////////////////////////////////////
// Exact/reduction verifier for the corrected resumable A2228 generator.
//
// It rechecks the integral tuple (after primitive cancellation), so rows
// passing through a base point of a homogeneous family profile cannot create
// false positives.  A row is rejected at the first good prime for which
// 3 does not divide #J(F_p).  Only full reduction survivors reach Magma's
// exact TorsionSubgroup computation.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

if not assigned input_file then
    input_file := "results/target_22224_a2228_deep13_generator_candidates_h2000.tsv";
end if;
if not assigned log_file then
    log_file := "results/target_22224_a2228_deep13_generator_verify_h2000.log";
end if;
if not assigned hits_file then
    hits_file := "results/target_22224_a2228_deep13_generator_hits_h2000.tsv";
end if;
if not assigned prime_bound then prime_bound := 499;
elif Type(prime_bound) eq MonStgElt then prime_bound := StringToInteger(prime_bound);
end if;

SetLogFile(log_file : Overwrite := true);
Q := Rationals(); Z := Integers(); P<x> := PolynomialRing(Q);

function ReadCandidates(filename)
    out := [];
    for row in Split(Read(filename),"\n") do
        if #row eq 0 or row[1] notin {"1","2","3"} then continue; end if;
        s := Split(row,"\t");
        if #s lt 8 then continue; end if;
        Append(~out,<StringToInteger(s[1]),StringToInteger(s[3]),
                     StringToInteger(s[4]),
                     [StringToInteger(s[i]):i in [5..8]]>);
    end for;
    return out;
end function;

function CoverRadicands(vals)
    a,b,c,d := Explode(vals);
    return [a*b*c*d,
            a*(a+b)*(a+c)*(a+d),
            b*(b+a)*(b+c)*(b+d),
            c*(c+a)*(c+b)*(c+d)];
end function;

function CurvePolynomial(vals)
    return x*&*[x+(Q!z)^2:z in vals];
end function;

function StrictSimpleWitness(f)
    for p in PrimesUpTo(199) do
        try
            fp := ChangeRing(f,GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
            lp := LPolynomial(HyperellipticCurve(fp));
            R<T> := PolynomialRing(Q); chi := R!lp;
            if not IsIrreducible(chi) then continue; end if;
            K<pi> := NumberField(chi);
            if &and[Degree(MinimalPolynomial(pi^n)) eq 4:n in [2..12]] then
                return true,p,chi;
            end if;
        catch e
            continue;
        end try;
    end for;
    R<T> := PolynomialRing(Q);
    return false,0,R!0;
end function;

rows := ReadCandidates(input_file);
primes := [p:p in PrimesUpTo(prime_bound)];
elim := AssociativeArray(); for p in primes do elim[p] := 0; end for;
survivors := 0; exact_hits := 0;
hout := Open(hits_file,"w");
fprintf hout,"family\tm\tn\ta\tb\tc\td\ttorsion\tsimple\tcertificate_prime\n";

print "A2228_DEEP13_VERIFY_START","input",input_file,"rows",#rows,
      "prime_bound",prime_bound;
for i in [1..#rows] do
    fam,m,n,vals := Explode(rows[i]);
    assert &and[z ne 0:z in vals];
    assert #Set([z^2:z in vals]) eq 4;
    assert &and[IsSquare(Q!z):z in CoverRadicands(vals)];
    f := CurvePolynomial(vals);
    killed := false; kp := 0; kn := 0;
    for p in primes do
        try
            fp := ChangeRing(f,GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
            nord := Z!Evaluate(LPolynomial(HyperellipticCurve(fp)),1);
            if nord mod 3 ne 0 then
                killed := true; kp := p; kn := nord; elim[p] +:= 1; break;
            end if;
        catch e
            continue;
        end try;
    end for;
    if killed then
        if i le 20 then
            print "KILLED","row",i,"family",fam,"t",m,"/",n,
                  "prime",kp,"Jorder",kn;
        end if;
        continue;
    end if;
    survivors +:= 1;
    G,mp := TorsionSubgroup(Jacobian(HyperellipticCurve(f)));
    inv := Invariants(G);
    print "REDUCTION_SURVIVOR","row",i,"family",fam,"t",m,"/",n,
          "tuple",vals,"torsion",inv;
    if #G mod 3 eq 0 then
        exact_hits +:= 1;
        simple,sp,chi := StrictSimpleWitness(f);
        print "TARGET_22224_HIT","row",i,"family",fam,"t",m,"/",n,
              "tuple",vals,"torsion",inv,"simple",simple,
              "certificate_prime",sp,"chi",chi;
        fprintf hout,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",
            fam,m,n,vals[1],vals[2],vals[3],vals[4],Sprint(inv),
            simple select 1 else 0,sp;
    end if;
end for;

delete hout;
print "ELIMINATIONS",[<p,elim[p]>:p in primes|elim[p] gt 0];
print "A2228_DEEP13_VERIFY_DONE","rows",#rows,
      "reduction_survivors",survivors,"exact_hits",exact_hits,
      "hits_file",hits_file;
UnsetLogFile();
quit;
