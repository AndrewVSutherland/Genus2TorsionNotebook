//////////////////////////////////////////////////////////////////////
// Rational-curve search on the full A(2,2,2,8) cover for 3-torsion.
//
// These are genuine full-cover curves, not the special K3 locus used by
// m2228_three_torsion_search.m.  We test Filip's first two explicit
// one-parameter families and the recovered Adam family.  Rational
// 3-torsion together with [2,2,2,8] gives [2,2,2,24].
//
// Run from torsion_jac:
//
//   magma -b height:=100 code/target_22224_a2228_curves_plus3.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

if not assigned height then height := 100;
elif Type(height) eq MonStgElt then height := StringToInteger(height);
end if;
if not assigned log_file then
    log_file := Sprintf("results/target_22224_a2228_curves_plus3_h%o.log",height);
end if;
if not assigned progress_interval then progress_interval := 10000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

SetLogFile(log_file : Overwrite := true);

Q := Rationals(); Z := Integers();
P<X> := PolynomialRing(Q);

function HeightRationals(H)
    vals := {}; 
    for den in [1..H] do
        for num in [-H..H] do
            if GCD(num,den) eq 1 then Include(~vals,Q!num/den); end if;
        end for;
    end for;
    return Sort(Setseq(vals));
end function;

function Filip1(t)
    return [
        -(t^2+t+1)^2,
        -4*t*(t+1)^2,
        4*t*(t+1),
        4*t^2*(t+1)
    ];
end function;

function Filip2(t)
    den := t^4+2*t^3-t^2-2*t+1;
    if t eq 0 or den eq 0 then return false,[]; end if;
    return true,[
        -(t^4-2*t^3-t^2+2*t+1)/den,
        -1/t,
        Q!1,
        t
    ];
end function;

function Adam(t)
    if t in {Q!-1,Q!0,Q!1} then return false,[]; end if;
    return true,[
        -t*(t+1)/(t-1),
        Q!1,
        -t^2,
        t*(t-1)/(t+1)
    ];
end function;

function CoverRadicands(vals)
    a,b,c,d := Explode(vals);
    return [
        a*b*c*d,
        a*(a+b)*(a+c)*(a+d),
        b*(b+a)*(b+c)*(b+d),
        c*(c+a)*(c+b)*(c+d)
    ];
end function;

function RationalSquare(r)
    if r lt 0 then return false; end if;
    return IsSquare(Z!Numerator(r)) and IsSquare(Z!Denominator(r));
end function;

function CurvePolynomial(vals)
    return X*&*[X+z^2:z in vals];
end function;

function Usable(vals)
    if #vals ne 4 or &or[z eq 0:z in vals] then return false; end if;
    sq := [z^2:z in vals];
    return #Set(sq) eq 4;
end function;

prime_list := [13,11,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,
               79,83,89,97,101,103,107,109,113,127,131,137,139,149,
               151,157,163,167,173,179,181,191,193,197,199];

function ThreeBound(f)
    C := HyperellipticCurve(f); g := 0; used := [];
    for p in prime_list do
        try
            fp := ChangeRing(f,GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
            n := Z!Evaluate(LPolynomial(HyperellipticCurve(fp)),1);
            g := (g eq 0) select n else GCD(g,n);
            Append(~used,<p,n,g>);
            if g mod 3 ne 0 then return false,g,used; end if;
        catch e
            continue;
        end try;
    end for;
    return g ne 0 and g mod 3 eq 0,g,used;
end function;

function SimpleWitness(f)
    C := HyperellipticCurve(f);
    for p in prime_list do
        try
            fp := ChangeRing(f,GF(p));
            if Discriminant(fp) eq 0 then continue; end if;
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

params := HeightRationals(height);
seen := {}; tested := 0; usable := 0; cover := 0; survivors := 0; hits := 0;
elim := AssociativeArray(); for p in prime_list do elim[p] := 0; end for;

print "TARGET_22224_A2228_CURVES_PLUS3_START","height",height,"parameters",#params;
for family in [1..3] do
    fam_tested := 0; fam_usable := 0; fam_cover := 0; fam_survivors := 0;
    for t in params do
        ok := true; vals := [];
        if family eq 1 then
            if t in {Q!0,Q!-1} then continue; end if;
            vals := Filip1(t);
        elif family eq 2 then
            ok,vals := Filip2(t);
        else
            ok,vals := Adam(t);
        end if;
        if not ok then continue; end if;
        tested +:= 1; fam_tested +:= 1;
        if not Usable(vals) then continue; end if;
        usable +:= 1; fam_usable +:= 1;
        rads := CoverRadicands(vals);
        if not &and[r ne 0 and RationalSquare(r):r in rads] then continue; end if;
        cover +:= 1; fam_cover +:= 1;
        f := CurvePolynomial(vals);
        // Deduplicate isomorphic common-scale repetitions by the monic
        // coefficient string after translating no branch: family symmetries
        // may still survive, but never create false hits.
        key := Sprint(f);
        if key in seen then continue; end if;
        Include(~seen,key);
        survives,g,used := ThreeBound(f);
        if not survives then
            elim[used[#used][1]] +:= 1;
            continue;
        end if;
        survivors +:= 1; fam_survivors +:= 1;
        print "REDUCTION_SURVIVOR","family",family,"t",t,"vals",vals,
              "gcd_bound",g,"used",used;
        G,mp := TorsionSubgroup(Jacobian(HyperellipticCurve(f)));
        inv := Invariants(G);
        print "EXACT","family",family,"t",t,"torsion",inv,"curve",f;
        if #G mod 3 eq 0 then
            hits +:= 1;
            simple,p,chi := SimpleWitness(f);
            print "TARGET_22224_HIT","family",family,"t",t,"vals",vals,
                  "torsion",inv,"simple",simple,"prime",p,"chi",chi;
        end if;
        if progress_interval gt 0 and tested mod progress_interval eq 0 then
            print "PROGRESS",tested,"cover",cover,"survivors",survivors,"hits",hits;
        end if;
    end for;
    print "FAMILY_DONE",family,"tested",fam_tested,"usable",fam_usable,
          "cover",fam_cover,"reduction_survivors",fam_survivors;
end for;

print "ELIMINATIONS",[<p,elim[p]>:p in prime_list|elim[p] gt 0];
print "TARGET_22224_A2228_CURVES_PLUS3_DONE","tested",tested,"usable",usable,
      "cover",cover,"unique",#seen,"reduction_survivors",survivors,"hits",hits;
UnsetLogFile();
quit;
