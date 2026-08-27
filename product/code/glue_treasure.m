// glue_treasure.m — process (N,N)-congruent treasure pairs: analytic gluing
// (analytic_glue.m) followed by certification: L_p(C) = L_p(E) L_p(F) at good
// p < 120, and exact TorsionSubgroup of the glued Jacobian (expected to
// contain T1 x T2).  Edit TREASURES below (source: lane_congtors.log /
// lane_cong_sieve.log), then run:
//   cd product/code && magma -b glue_treasure.m > ../logs/glue_treasure.log
SetColumns(0);
SetMemoryLimit(8*10^9);
load "split_lab.m";
load "analytic_glue.m";

// <label-or-ainvs-string E, same F, N, comment>
TREASURES := [*
    // <"11a3", "66c1", 5, "smoke test: [5] x [1] -- nothing to gain, pipeline check">
*];

for tr in TREASURES do
    E := Type(tr[1]) eq MonStgElt select EllipticCurve(tr[1]) else EllipticCurve(tr[1]);
    F := Type(tr[2]) eq MonStgElt select EllipticCurve(tr[2]) else EllipticCurve(tr[2]);
    N := tr[3];
    T1 := Invariants(TorsionSubgroup(E)); T2 := Invariants(TorsionSubgroup(F));
    expected := Invariants(AbelianGroup(T1 cat T2));
    printf "==== TREASURE %o x %o (N=%o): T1=%o T2=%o expect >= %o ====\n",
        tr[1], tr[2], N, T1, T2, expected;
    cands := [* *];
    for prec in [160, 320] do
        cands := GlueCandidates(E, F, N : prec := prec);
        if #cands gt 0 then break; end if;
        printf "no candidates at prec %o, retrying higher\n", prec;
    end for;
    for c in cands do
        C := c[1];
        ok, g := IntegralSextic(C);
        if not ok then printf "TREASURE: bad model\n"; continue; end if;
        // certification 1: L-poly split at good p
        dz := Integers()!Discriminant(g);
        lc := Integers()!LeadingCoefficient(g);
        bad := 2*dz*lc*Conductor(E)*Conductor(F);
        nok := 0; nbad := 0;
        for p in PrimesUpTo(120) do
            if bad mod p eq 0 then continue; end if;
            gp := PolynomialRing(GF(p))!g;
            Lp := Numerator(ZetaFunction(HyperellipticCurve(gp)));
            PZ<T> := Parent(Lp);
            LE := 1 - TraceOfFrobenius(E,p)*T + p*T^2;
            LF := 1 - TraceOfFrobenius(F,p)*T + p*T^2;
            if Lp eq LE*LF then nok +:= 1; else nbad +:= 1; end if;
        end for;
        printf "TREASURE M=%o twist=%o Lp-split %o ok / %o bad\n", c[2], c[3], nok, nbad;
        if nbad gt 0 or nok lt 5 then printf "TREASURE: NOT the right surface, skip\n"; continue; end if;
        einv := [Integers()|-1];
        try
            einv := ExactTorsion(g);
        catch e
            printf "TREASURE: ExactTorsion failed\n"; continue;
        end try;
        isnew := not einv in KNOWN;
        printf "TREASURE_RESULT invs=%o %o (expected >= %o) g=%o\n",
            einv, isnew select "*** NEW GROUP ***" else "known", expected, g;
    end for;
end for;
printf "TREASURE_DONE\n";
quit;
