//////////////////////////////////////////////////////////////////////
// opus_abc_classify.m  (2026-07-19)
//
// Takes (A,B,C) points on the (2,2,2,12) moduli surface, maps them
// through (s,m,n) to the genus-2 curve, and reports exact torsion +
// absolute G2 invariants, so that points can be grouped into
// ISOMORPHISM classes of curves (the real "up to automorphism" test).
//
// Also runs a two-prime root-power geometric simplicity certificate
// on each distinct curve.
//
// Run: magma -b pts:="60,109,133;120,218,143;..." \
//        code/opus_abc_classify.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(2*10^9);

QQ := Rationals();
Z  := Integers();
PP<x> := PolynomialRing(QQ);

if not assigned pts then
    pts := "60,109,133;120,218,143;120,241,143;120,241,266;143,218,266;143,241,266;143,437,408";
end if;

function smnFromABC(a, b, c)
    den1 := 2*a^2*b - 2*b^3 + 2*b*c^2;
    den2 := 2*a^3*b - 2*a*b^3 + 2*a*b*c^2;
    if den1 eq 0 or den2 eq 0 or b eq 0 then return false, 0, 0, 0; end if;
    s := (a^3 - a*b^2 + 2*a*c^2) / den1;
    m := (-2*a^2*b^2 + 2*b^4 + 2*a^2*c^2 - 6*b^2*c^2 + 4*c^4) / den2;
    n := a / b;
    return true, s, m, n;
end function;

function CurveFromSMN(s, m, n)
    Bc := [ 2*s^2 - s*n,
            2*s^2 + s*m - 2*s*n - m*n,
            2*s^2 + s*m - s*n - m*n,
            -m*n,
            4*s^2 - 4*s*n - m*n ];
    Ac := [1,1,1,2,2];
    return PP!(&*[ Ac[i] + Bc[i]*x : i in [1..5] ]);
end function;

// integral square-scaled model: y^2 = f  ->  (den*y)^2 = den^2*f
function IntModel(f)
    den := LCM([Denominator(Coefficient(f,i)) : i in [0..Degree(f)]]);
    fI := PP!(den^2*f);
    return PP![Z!co : co in Coefficients(fI)];
end function;

// two-prime root-power geometric simplicity certificate
function SimplicityCert(f)
    fI := IntModel(f);
    C := HyperellipticCurve(fI);
    D := Z!Discriminant(C);
    ncert := 0; witness := [];
    for p in [11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,137,149] do
        if ncert ge 2 then break; end if;
        if D mod p eq 0 then continue; end if;
        ok := true;
        try
            chi := PP!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(C, GF(p))))));
        catch e ok := false; end try;
        if not ok or not IsIrreducible(chi) then continue; end if;
        K<pi> := NumberField(chi);
        good := true;
        for nn in [2..12] do
            if Degree(MinimalPolynomial(pi^nn)) ne 4 then good := false; break; end if;
        end for;
        if good then ncert +:= 1; Append(~witness, <p, chi>); end if;
    end for;
    return ncert, witness;
end function;

seen := AssociativeArray();
order := [];
for chunk in Split(pts, ";") do
    cs := [StringToInteger(t) : t in Split(chunk, ",")];
    a := cs[1]; b := cs[2]; c := cs[3];
    ok, s, m, n := smnFromABC(a, b, c);
    if not ok then printf "POINT %o,%o,%o : degenerate map\n", a, b, c; continue; end if;
    f := CurveFromSMN(s, m, n);
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        printf "POINT %o,%o,%o : singular curve\n", a, b, c; continue;
    end if;
    Cv := HyperellipticCurve(f);
    g2 := G2Invariants(Cv);
    key := Sprintf("%o", g2);
    T := TorsionSubgroup(Jacobian(HyperellipticCurve(IntModel(f))));
    printf "POINT %o,%o,%o (s,m,n)=(%o,%o,%o) torsion=%o\n", a, b, c, s, m, n, Invariants(T);
    if IsDefined(seen, key) then
        printf "   SAME CURVE as %o\n", seen[key];
    else
        seen[key] := Sprintf("%o,%o,%o", a, b, c);
        Append(~order, <key, a, b, c, f, Invariants(T)>);
        printf "   NEW CURVE class; G2=%o\n", g2;
    end if;
end for;

printf "\nDISTINCT CURVES: %o\n", #order;
for t in order do
    printf "CLASS from (%o,%o,%o) torsion=%o\n", t[2], t[3], t[4], t[6];
    nc, wit := SimplicityCert(t[5]);
    printf "   simplicity certificates: %o", nc;
    if nc ge 2 then printf "  (p=%o, p=%o)", wit[1][1], wit[2][1]; end if;
    printf "\n";
    printf "   f = %o\n", t[5];
end for;
printf "CLASSIFY_DONE\n";
quit;
