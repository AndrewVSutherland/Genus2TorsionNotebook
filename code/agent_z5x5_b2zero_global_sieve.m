//////////////////////////////////////////////////////////////////////
//  Z/5 x Z/5 full Mumford norm, b2=0 branch:
//  global/local sieve from smooth F7/F11 p-adic charts.
//
//  This file is intentionally standalone and writes no files.  Typical
//  bounded run from the repository root:
//
//      magma \
//          code/agent_z5x5_b2zero_global_sieve.m \
//          > results/z5x5_b2zero_global_default.log 2>&1
//////////////////////////////////////////////////////////////////////

SetColumns(0);

function ToIntegerParameter(x)
    if Type(x) eq MonStgElt then
        return StringToInteger(x);
    end if;
    return x;
end function;

function ToBoolParameter(x)
    if Type(x) eq MonStgElt then
        return (x eq "true") or (x eq "True") or (x eq "1");
    end if;
    return x;
end function;

if not assigned low_7_power then
    low_7_power := 4;
else
    low_7_power := ToIntegerParameter(low_7_power);
end if;

if not assigned low_11_power then
    low_11_power := 3;
else
    low_11_power := ToIntegerParameter(low_11_power);
end if;

if not assigned high_7_power then
    high_7_power := 7;
else
    high_7_power := ToIntegerParameter(high_7_power);
end if;

if not assigned high_11_power then
    high_11_power := 6;
else
    high_11_power := ToIntegerParameter(high_11_power);
end if;

if not assigned frontier_height then
    frontier_height := 500000;
else
    frontier_height := ToIntegerParameter(frontier_height);
end if;

if not assigned search_num_bound then
    search_num_bound := 5000;
else
    search_num_bound := ToIntegerParameter(search_num_bound);
end if;

if not assigned search_den_bound then
    search_den_bound := 5000;
else
    search_den_bound := ToIntegerParameter(search_den_bound);
end if;

if not assigned tall_num_bound then
    tall_num_bound := 20000;
else
    tall_num_bound := ToIntegerParameter(tall_num_bound);
end if;

if not assigned tall_den_bound then
    tall_den_bound := 1000;
else
    tall_den_bound := ToIntegerParameter(tall_den_bound);
end if;

if not assigned wide_num_bound then
    wide_num_bound := 1000;
else
    wide_num_bound := ToIntegerParameter(wide_num_bound);
end if;

if not assigned wide_den_bound then
    wide_den_bound := 20000;
else
    wide_den_bound := ToIntegerParameter(wide_den_bound);
end if;

if not assigned max_search_tuples then
    max_search_tuples := 20000000;
else
    max_search_tuples := ToIntegerParameter(max_search_tuples);
end if;

if not assigned lift_count_cap then
    lift_count_cap := 10000;
else
    lift_count_cap := ToIntegerParameter(lift_count_cap);
end if;

if not assigned do_low_crt_search then
    do_low_crt_search := true;
else
    do_low_crt_search := ToBoolParameter(do_low_crt_search);
end if;

if not assigned do_frontier then
    do_frontier := true;
else
    do_frontier := ToBoolParameter(do_frontier);
end if;

if not assigned do_adjacent_probe then
    do_adjacent_probe := true;
else
    do_adjacent_probe := ToBoolParameter(do_adjacent_probe);
end if;

if not assigned adjacent_power then
    adjacent_power := 3;
else
    adjacent_power := ToIntegerParameter(adjacent_power);
end if;

Names5 := ["K", "s", "t", "b0", "b1"];

function ForcedAData(h, kap, U, B)
    P := Parent(h);
    x := P.1;
    f := h^2 - kap*x^5;
    A := x^5;
    E := A^2 - B^2*f - U^5;
    avec := [];

    for deg in [9,8,7,6,5] do
        a := -Coefficient(E, deg)/2;
        Append(~avec, a);
        A +:= a*x^(deg - 5);
        E := A^2 - B^2*f - U^5;
    end for;

    return A, E, f, avec;
end function;

function SliceResidualsQ(h1fix, h2fix)
    Q := Rationals();
    R<kap,s,t,b0,b1> := PolynomialRing(Q, 5, "grevlex");
    Fr := FieldOfFractions(R);
    P<x> := PolynomialRing(Fr);

    h := 1 + (Fr!h1fix)*x + (Fr!h2fix)*x^2;
    U := x^2 + s*x + t;
    B := b0 + b1*x;
    A, E, f, avec := ForcedAData(h, kap, U, B);
    res := [R!Numerator(Coefficient(E, i)) : i in [0..4]];
    return R, res, A, f, U, B;
end function;

function IntegralMultiple(f)
    RQ := Parent(f);
    nvars := Rank(RQ);
    Z := Integers();
    RZ := PolynomialRing(Z, nvars, "grevlex");
    coeffs := Coefficients(f);
    mons := Monomials(f);
    den := 1;
    for c in coeffs do
        den := LCM(den, Denominator(c));
    end for;

    g := RZ!0;
    for i in [1..#coeffs] do
        exps := Exponents(mons[i]);
        term := RZ!(Z!(den*coeffs[i]));
        for j in [1..nvars] do
            term *:= RZ.j^exps[j];
        end for;
        g +:= term;
    end for;
    return g, den;
end function;

function IntegralResiduals(resQ)
    resZ := [];
    dens := [];
    for e in resQ do
        g, den := IntegralMultiple(e);
        Append(~resZ, g);
        Append(~dens, den);
    end for;
    return resZ, dens;
end function;

function ResidualsZeroMod(resZ, vals, m)
    Z := Integers();
    for e in resZ do
        if (Z!Evaluate(e, vals)) mod m ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function ResidualsZeroExact(resZ, vals)
    Q := Rationals();
    for e in resZ do
        if Q!Evaluate(e, vals) ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function BoundaryValuesQ(vals, h1fix, h2fix)
    Q := Rationals();
    P<x> := PolynomialRing(Q);
    kap := Q!vals[1]; s := Q!vals[2]; t := Q!vals[3];
    b0 := Q!vals[4]; b1 := Q!vals[5];
    h := 1 + (Q!h1fix)*x + (Q!h2fix)*x^2;
    f := h^2 - kap*x^5;
    U := x^2 + s*x + t;
    B := b0 + b1*x;
    discU := Discriminant(U);
    resBU := Resultant(B, U);
    discF := Discriminant(f);
    return [kap, b1, b0, discU, resBU, discF];
end function;

function IsOpenQ(vals, h1fix, h2fix)
    return &and[b ne 0 : b in BoundaryValuesQ(vals, h1fix, h2fix)];
end function;

function ResidualJacobianRankModP(resZ, vals, p)
    F := GF(p);
    Z := Integers();
    entries := [];
    for i in [1..#resZ] do
        for j in [1..5] do
            Append(~entries, F!(Z!Evaluate(Derivative(resZ[i], j), vals) mod p));
        end for;
    end for;
    M := Matrix(F, #resZ, 5, entries);
    return Rank(M), Determinant(M);
end function;

function BalancedResidue(a, m)
    Z := Integers();
    r := (Z!a) mod m;
    if 2*r gt m then
        return r - m;
    end if;
    return r;
end function;

function BalancedVector(vals, m)
    return [BalancedResidue(Integers()!v, m) : v in vals];
end function;

function VectorMod(vals, m)
    return [(Integers()!v) mod m : v in vals];
end function;

function RelationPairs(J, D0, D2)
    rels := [];
    for a in [0..4] do
        for b in [0..4] do
            if a eq 0 and b eq 0 then
                continue;
            end if;
            if a*D0 + b*D2 eq J!0 then
                Append(~rels, <a,b>);
            end if;
        end for;
    end for;
    return rels;
end function;

function CheckClassesQ(f, h, U, V)
    P := Parent(f);
    x := P.1;
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    D0 := J![x, Evaluate(h, 0)];
    D2 := J![U, V];

    if D0 eq J!0 or D2 eq J!0 then
        return false, [], "zero class";
    end if;
    if 5*D0 ne J!0 then
        return false, [], "5*D0 != 0";
    end if;
    if 5*D2 ne J!0 then
        return false, [], "5*D2 != 0";
    end if;

    rels := RelationPairs(J, D0, D2);
    return (#rels eq 0), rels, "checked";
end function;

function B2ZeroModelQ(vals, h1fix, h2fix)
    Q := Rationals();
    P<x> := PolynomialRing(Q);
    h := 1 + (Q!h1fix)*x + (Q!h2fix)*x^2;
    kap := Q!vals[1]; s := Q!vals[2]; t := Q!vals[3];
    b0 := Q!vals[4]; b1 := Q!vals[5];
    U := x^2 + s*x + t;
    B := b0 + b1*x;
    A, E, f, avec := ForcedAData(h, kap, U, B);

    if &or[Coefficient(E, i) ne 0 : i in [0..4]] then
        return false, A, B, f, U, P!0, [], "residual nonzero";
    end if;
    if not IsOpenQ(vals, h1fix, h2fix) then
        return false, A, B, f, U, P!0, [], "outside open chart";
    end if;
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false, A, B, f, U, P!0, [], "bad contact curve";
    end if;
    if Discriminant(U) eq 0 or GCD(U, f) ne 1 then
        return false, A, B, f, U, P!0, [], "bad U";
    end if;
    if GCD(B, U) ne 1 then
        return false, A, B, f, U, P!0, [], "B not invertible mod U";
    end if;

    Binv := InverseMod(B, U);
    V := (-A*Binv) mod U;
    if (V^2 - f) mod U ne 0 then
        return false, A, B, f, U, V, [], "V recovery failed";
    end if;

    independent, rels, reason := CheckClassesQ(f, h, U, V);
    if not independent then
        return false, A, B, f, U, V, rels, reason;
    end if;
    return true, A, B, f, U, V, rels, "independent";
end function;

function LiftSliceChart(label, p, h1fix, h2fix, base5, max_power, resZ)
    printf "\n# lift_chart %o p=%o h1=%o h2=%o max_power=%o\n",
        label, p, h1fix, h2fix, max_power;
    printf "  base_mod_%o=<K=%o,s=%o,t=%o,b0=%o,b1=%o>\n",
        p, base5[1], base5[2], base5[3], base5[4], base5[5];
    printf "  base_boundary_values=%o\n",
        [Integers()!b mod p : b in BoundaryValuesQ(base5, h1fix, h2fix)];

    rank, det := ResidualJacobianRankModP(resZ, base5, p);
    printf "  slice_jacobian_rank_mod_%o=%o determinant=%o\n", p, rank, det;

    sols := [VectorMod(base5, p)];
    printf "  level=1 modulus=%o lift_count=%o sample=%o balanced=%o\n",
        p, #sols, sols[1], BalancedVector(sols[1], p);

    for e in [2..max_power] do
        oldmod := p^(e - 1);
        newmod := p^e;
        newsols := [];
        capped := false;
        for sol in sols do
            for d1 in [0..p-1] do
            for d2 in [0..p-1] do
            for d3 in [0..p-1] do
            for d4 in [0..p-1] do
            for d5 in [0..p-1] do
                cand := [
                    (sol[1] + oldmod*d1) mod newmod,
                    (sol[2] + oldmod*d2) mod newmod,
                    (sol[3] + oldmod*d3) mod newmod,
                    (sol[4] + oldmod*d4) mod newmod,
                    (sol[5] + oldmod*d5) mod newmod
                ];
                if ResidualsZeroMod(resZ, cand, newmod) then
                    Append(~newsols, cand);
                    if #newsols ge lift_count_cap then
                        capped := true;
                        break d5;
                    end if;
                end if;
            end for;
                if capped then break d4; end if;
            end for;
                if capped then break d3; end if;
            end for;
                if capped then break d2; end if;
            end for;
                if capped then break d1; end if;
            end for;
            if capped then
                break;
            end if;
        end for;

        sols := newsols;
        if #sols gt 0 then
            rank_e, det_e := ResidualJacobianRankModP(resZ, sols[1], p);
            printf "  level=%o modulus=%o lift_count=%o capped=%o sample=%o balanced=%o jac_rank_mod_%o=%o determinant=%o\n",
                e, newmod, #sols, capped, sols[1], BalancedVector(sols[1], newmod),
                p, rank_e, det_e;
        else
            printf "  level=%o modulus=%o lift_count=0 capped=%o\n", e, newmod, capped;
            break;
        end if;
    end for;

    return sols, p^max_power;
end function;

function CRTPair(a, m, b, n)
    Zn := Integers(n);
    delta := Integers()!((Zn!(b - a))/(Zn!m));
    return (a + m*delta) mod (m*n);
end function;

function CRTVector(a, m, b, n)
    return [CRTPair(a[i], m, b[i], n) : i in [1..#a]];
end function;

function RationalCandidatesFast(residue, modulus, num_bound, den_bound)
    Q := Rationals();
    Z := Integers();
    r := (Z!residue) mod modulus;
    cands := [];
    seen := {};
    for den in [1..den_bound] do
        if GCD(den, modulus) ne 1 then
            continue;
        end if;
        n0 := BalancedResidue((r*(den mod modulus)) mod modulus, modulus);
        kmin := Ceiling(((-num_bound) - n0)/modulus);
        kmax := Floor((num_bound - n0)/modulus);
        for k in [kmin..kmax] do
            num := n0 + k*modulus;
            if Abs(num) le num_bound and GCD(Abs(num), den) eq 1 then
                q := (Q!num)/(Q!den);
                if not q in seen then
                    Include(~seen, q);
                    Append(~cands, q);
                end if;
            end if;
        end for;
    end for;
    Sort(~cands);
    return cands;
end function;

function ModResiduals(resZ, q)
    RZ := Parent(resZ[1]);
    S := PolynomialRing(GF(q), Rank(RZ), "grevlex");
    phi := hom<RZ -> S | [S.i : i in [1..Rank(RZ)]]>;
    return [phi(e) : e in resZ];
end function;

function RationalMod(qrat, q)
    F := GF(q);
    return (F!Integers()!Numerator(qrat))/(F!Integers()!Denominator(qrat));
end function;

function PassesModularResiduals(resMod, vals)
    for e in resMod do
        if Evaluate(e, vals) ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

procedure CRTSearch(label, modulus, residues, resZ, h1fix, h2fix, num_bound, den_bound)
    printf "\n# crt_rational_search %o modulus=%o num_bound=%o den_bound=%o\n",
        label, modulus, num_bound, den_bound;
    cand_lists := [
        RationalCandidatesFast(residues[i], modulus, num_bound, den_bound)
            : i in [1..5]
    ];
    counts := [#L : L in cand_lists];
    prod := 1;
    for c in counts do
        prod *:= c;
    end for;
    printf "  candidate_counts=<K=%o,s=%o,t=%o,b0=%o,b1=%o> tuple_count=%o\n",
        counts[1], counts[2], counts[3], counts[4], counts[5], prod;
    for i in [1..5] do
        preview := cand_lists[i];
        if #preview gt 8 then
            preview := preview[1..8];
        end if;
        printf "  candidates_%o_%o_preview=%o\n", i, Names5[i], preview;
    end for;

    if prod eq 0 then
        print "  search_result=no_coordinate_candidates";
        return;
    end if;
    if prod gt max_search_tuples then
        printf "  search_result=skipped_tuple_cap max_search_tuples=%o\n",
            max_search_tuples;
        return;
    end if;

    q1 := NextPrime(den_bound + 100);
    while q1 in {2,5,7,11} do
        q1 := NextPrime(q1 + 1);
    end while;
    q2 := NextPrime(q1 + 200);
    while q2 in {2,5,7,11} do
        q2 := NextPrime(q2 + 1);
    end while;
    sieve_primes := [q1, q2];
    printf "  sieve_primes=%o\n", sieve_primes;

    resMod1 := ModResiduals(resZ, q1);
    resMod2 := ModResiduals(resZ, q2);
    candMods1 := [[RationalMod(c, q1) : c in cand_lists[j]] : j in [1..5]];
    candMods2 := [[RationalMod(c, q2) : c in cand_lists[j]] : j in [1..5]];

    survivors := [0, 0];
    exact_tested := 0;
    exact_hits := 0;

    for i1 in [1..#cand_lists[1]] do
    for i2 in [1..#cand_lists[2]] do
    for i3 in [1..#cand_lists[3]] do
    for i4 in [1..#cand_lists[4]] do
    for i5 in [1..#cand_lists[5]] do
        idxs := [i1,i2,i3,i4,i5];
        passed_all := true;
        valsq1 := [candMods1[j][idxs[j]] : j in [1..5]];
        if PassesModularResiduals(resMod1, valsq1) then
            survivors[1] +:= 1;
        else
            continue;
        end if;

        valsq2 := [candMods2[j][idxs[j]] : j in [1..5]];
        if PassesModularResiduals(resMod2, valsq2) then
            survivors[2] +:= 1;
        else
            continue;
        end if;

        exact_tested +:= 1;
        vals := [
            cand_lists[1][i1], cand_lists[2][i2], cand_lists[3][i3],
            cand_lists[4][i4], cand_lists[5][i5]
        ];
        if ResidualsZeroExact(resZ, vals) then
            exact_hits +:= 1;
            ok, A, B, f, U, V, rels, reason := B2ZeroModelQ(vals, h1fix, h2fix);
            printf "  exact_hit_%o vals=<K=%o,s=%o,t=%o,b0=%o,b1=%o> class=%o\n",
                exact_hits, vals[1], vals[2], vals[3], vals[4], vals[5], reason;
            printf "    f=%o\n", f;
            printf "    U=%o\n", U;
            printf "    V=%o\n", V;
            printf "    A=%o\n", A;
            printf "    B=%o\n", B;
            if #rels gt 0 then
                printf "    rels=%o\n", rels;
            end if;
        end if;
    end for;
    end for;
    end for;
    end for;
    end for;

    printf "  modular_survivors_after_prefix_primes=%o exact_tested=%o exact_hits=%o\n",
        survivors, exact_tested, exact_hits;
    if exact_hits eq 0 then
        print "  search_result=no_exact_rational_candidate";
    end if;
end procedure;

procedure CRTFrontier(label, modulus, residues, max_height)
    Z := Integers();
    printf "\n# crt_height_frontier %o modulus=%o max_height=%o\n",
        label, modulus, max_height;
    all_found := true;
    min_heights := [];
    for i in [1..5] do
        r := (Z!residues[i]) mod modulus;
        bestH := max_height + 1;
        bestNum := 0;
        bestDen := 0;
        for den in [1..max_height] do
            if GCD(den, modulus) ne 1 then
                continue;
            end if;
            num := BalancedResidue((r*(den mod modulus)) mod modulus, modulus);
            if GCD(Abs(num), den) ne 1 then
                continue;
            end if;
            H := Max(Abs(num), den);
            if H lt bestH then
                bestH := H;
                bestNum := num;
                bestDen := den;
            end if;
        end for;
        if bestH gt max_height then
            all_found := false;
            Append(~min_heights, max_height + 1);
            printf "  coord_%o_%o min_height_le_%o=none\n",
                i, Names5[i], max_height;
        else
            Append(~min_heights, bestH);
            printf "  coord_%o_%o min_height_le_%o=%o best=%o/%o\n",
                i, Names5[i], max_height, bestH, bestNum, bestDen;
        end if;
    end for;
    if all_found then
        printf "  tuple_height_lower_bound_from_coordinates=%o\n", Max(min_heights);
    else
        printf "  tuple_height_lower_bound_from_coordinates=>%o\n", max_height;
    end if;
end procedure;

Q := Rationals();
R, resQ, Atemplate, ftemplate, Utemplate, Btemplate := SliceResidualsQ(Q!1, Q!0);
resZ, dens := IntegralResiduals(resQ);

print "# Z/5 x Z/5 b2=0 global/local sieve driver";
printf "parameters low_7_power=%o low_11_power=%o high_7_power=%o high_11_power=%o frontier_height=%o\n",
    low_7_power, low_11_power, high_7_power, high_11_power, frontier_height;
printf "parameters search_box=(%o,%o) tall_box=(%o,%o) wide_box=(%o,%o) max_search_tuples=%o\n",
    search_num_bound, search_den_bound, tall_num_bound, tall_den_bound,
    wide_num_bound, wide_den_bound, max_search_tuples;
printf "integral_residual_denominators=%o\n", dens;
for i in [1..#resZ] do
    printf "  E%o degree=%o total_degree=%o terms=%o\n",
        i - 1, Degree(resZ[i]), TotalDegree(resZ[i]), #Terms(resZ[i]);
end for;

// Open h1=1,h2=0 slice points from agent_z5x5_b2zero_elim/lift.
F7A := [2,0,6,3,1];
F7B := [2,0,6,4,6];
F11A := [7,5,10,2,4];
F11B := [7,5,10,9,7];

sols7A, mod7High := LiftSliceChart("F7_slice_A", 7, Q!1, Q!0, F7A, high_7_power, resZ);
sols7B, mod7HighB := LiftSliceChart("F7_slice_B", 7, Q!1, Q!0, F7B, high_7_power, resZ);
sols11A, mod11High := LiftSliceChart("F11_slice_A", 11, Q!1, Q!0, F11A, high_11_power, resZ);
sols11B, mod11HighB := LiftSliceChart("F11_slice_B", 11, Q!1, Q!0, F11B, high_11_power, resZ);

mod7Low := 7^low_7_power;
mod11Low := 11^low_11_power;
modLow := mod7Low*mod11Low;
modHigh := mod7High*mod11High;

if #sols7A gt 0 and #sols7B gt 0 and #sols11A gt 0 and #sols11B gt 0 then
    C7 := [
        <"F7A", VectorMod(sols7A[1], mod7Low), VectorMod(sols7A[1], mod7High)>,
        <"F7B", VectorMod(sols7B[1], mod7Low), VectorMod(sols7B[1], mod7High)>
    ];
    C11 := [
        <"F11A", VectorMod(sols11A[1], mod11Low), VectorMod(sols11A[1], mod11High)>,
        <"F11B", VectorMod(sols11B[1], mod11Low), VectorMod(sols11B[1], mod11High)>
    ];

    printf "\n# crt_summary low_modulus=%o high_modulus=%o\n", modLow, modHigh;
    for c7 in C7 do
        for c11 in C11 do
            label := "CRT_" cat c7[1] cat "_" cat c11[1];
            lowResidues := CRTVector(c7[2], mod7Low, c11[2], mod11Low);
            highResidues := CRTVector(c7[3], mod7High, c11[3], mod11High);
            printf "  %o low_balanced=%o high_balanced=%o\n",
                label, BalancedVector(lowResidues, modLow),
                BalancedVector(highResidues, modHigh);

            if do_low_crt_search then
                CRTSearch(label cat "_balanced", modLow, lowResidues, resZ, Q!1, Q!0,
                    search_num_bound, search_den_bound);
                CRTSearch(label cat "_tall", modLow, lowResidues, resZ, Q!1, Q!0,
                    tall_num_bound, tall_den_bound);
                CRTSearch(label cat "_wide", modLow, lowResidues, resZ, Q!1, Q!0,
                    wide_num_bound, wide_den_bound);
            end if;

            if do_frontier then
                CRTFrontier(label cat "_high", modHigh, highResidues, frontier_height);
            end if;
        end for;
    end for;
else
    print "# crt_summary skipped because at least one selected chart failed to lift";
end if;

if do_adjacent_probe then
    print "\n# adjacent_local_chart_probe";

    Radj7, resQadj7, Aadj7, fadj7, Uadj7, Badj7 := SliceResidualsQ(Q!2, Q!0);
    resZadj7, densAdj7 := IntegralResiduals(resQadj7);
    adj7 := [1,0,5,1,3];
    solsAdj7, modAdj7 := LiftSliceChart("F7_adjacent_h1_2_h2_0", 7, Q!2, Q!0,
        adj7, adjacent_power, resZadj7);
    printf "  adjacent_F7_integral_residual_denominators=%o final_count=%o final_modulus=%o\n",
        densAdj7, #solsAdj7, modAdj7;

    Radj11, resQadj11, Aadj11, fadj11, Uadj11, Badj11 := SliceResidualsQ(Q!1, Q!1);
    resZadj11, densAdj11 := IntegralResiduals(resQadj11);
    adj11 := [4,3,2,2,10];
    solsAdj11, modAdj11 := LiftSliceChart("F11_adjacent_h1_1_h2_1", 11, Q!1, Q!1,
        adj11, adjacent_power, resZadj11);
    printf "  adjacent_F11_integral_residual_denominators=%o final_count=%o final_modulus=%o\n",
        densAdj11, #solsAdj11, modAdj11;
end if;

print "\n# verdict_marker";
print "global_sieve_completed";
print "check crt_rational_search exact_hits and crt_height_frontier lower bounds above";

quit;
