// verify_isoglue_hits.m — independent fresh-session verification of glued
// witnesses (2026-08-13 session).  For each entry: rebuild from raw sextic
// coefficients, recompute exact torsion on a minimal/integral model, check
// L_p(C) = L_p(E) L_p(F) at all good p < 200 against the STORED elliptic
// pair (labels or a-invariant lists), and report.  Populate HITS from
// lane_isoglue.log ISOGLUE_RESULT lines before running.
// Usage: cd product/code && magma -b verify_isoglue_hits.m > ../logs/verify_isoglue_hits.log
SetColumns(0);
SetMemoryLimit(8*10^9);
load "split_lab.m";

// < invs-claimed, f-coeffs (ascending), aInvariants E, aInvariants F, tag >
HITS := [*
*];

QZ<T> := PolynomialRing(Integers());
for h in HITS do
    g := RQx!h[2];
    E := EllipticCurve([Rationals()| a : a in h[3] ]);
    F := EllipticCurve([Rationals()| a : a in h[4] ]);
    printf "==== %o: claimed %o ====\n", h[5], h[1];
    einv := ExactTorsion(g);
    printf "exact torsion: %o  (match=%o)\n", einv, einv eq h[1];
    dz := Integers()!Discriminant(g);
    lc := Integers()!LeadingCoefficient(g);
    bad := 2*dz*lc*Conductor(E)*Conductor(F);
    nok := 0; nbad := 0;
    for p in PrimesUpTo(200) do
        if bad mod p eq 0 then continue; end if;
        gp := PolynomialRing(GF(p))!g;
        Lp := Numerator(ZetaFunction(HyperellipticCurve(gp)));
        PZ<t> := Parent(Lp);
        LE := 1 - TraceOfFrobenius(E,p)*t + p*t^2;
        LF := 1 - TraceOfFrobenius(F,p)*t + p*t^2;
        if Lp eq LE*LF then nok +:= 1; else nbad +:= 1; printf "  Lp MISMATCH p=%o\n", p; end if;
    end for;
    printf "Lp-split: %o ok, %o bad\n", nok, nbad;
    printf "VERDICT %o: %o\n", h[5],
        (einv eq h[1] and nbad eq 0 and nok ge 10) select "VERIFIED" else "PROBLEM";
end for;
printf "VERIFY_HITS_DONE\n";
quit;
