//////////////////////////////////////////////////////////////////////
//  Z/5 x Z/5 full Mumford norm, b2=0 branch: p-adic lifts and
//  bounded rational reconstruction on the h1=1, h2=0 slice.
//
//  This script is intentionally disjoint from
//  agent_z5x5_b2zero_elim.m: it reads no files and writes no files.
//  Typical run from the repository root:
//
//      magma \
//          code/agent_z5x5_b2zero_lift.m \
//          > results/z5x5_b2zero_lift_default.log 2>&1
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned max_7_power then
    max_7_power := 3;
elif Type(max_7_power) eq MonStgElt then
    max_7_power := StringToInteger(max_7_power);
end if;

if not assigned max_11_power then
    max_11_power := 3;
elif Type(max_11_power) eq MonStgElt then
    max_11_power := StringToInteger(max_11_power);
end if;

if not assigned rational_num_bound then
    rational_num_bound := 80;
elif Type(rational_num_bound) eq MonStgElt then
    rational_num_bound := StringToInteger(rational_num_bound);
end if;

if not assigned rational_den_bound then
    rational_den_bound := 40;
elif Type(rational_den_bound) eq MonStgElt then
    rational_den_bound := StringToInteger(rational_den_bound);
end if;

if not assigned max_search_tuples then
    max_search_tuples := 2000000;
elif Type(max_search_tuples) eq MonStgElt then
    max_search_tuples := StringToInteger(max_search_tuples);
end if;

if not assigned lift_count_cap then
    lift_count_cap := 200000;
elif Type(lift_count_cap) eq MonStgElt then
    lift_count_cap := StringToInteger(lift_count_cap);
end if;

if not assigned do_crt_search then
    do_crt_search := true;
elif Type(do_crt_search) eq MonStgElt then
    do_crt_search := (do_crt_search eq "true") or (do_crt_search eq "True");
end if;

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

function B2ZeroSliceResidualsQ()
    Q := Rationals();
    R<kap,s,t,b0,b1> := PolynomialRing(Q, 5, "grevlex");
    Fr := FieldOfFractions(R);
    P<x> := PolynomialRing(Fr);

    h := 1 + x;
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

function BoundaryValuesQ(vals)
    Q := Rationals();
    kap := Q!vals[1]; s := Q!vals[2]; t := Q!vals[3];
    b0 := Q!vals[4]; b1 := Q!vals[5];
    discU := s^2 - 4*t;
    resBU := -s*b0*b1 + t*b1^2 + b0^2;
    // On h1=1, h2=0, disc(f) = 108*K^3 + 3125*K^4.
    discF := 108*kap^3 + 3125*kap^4;
    return [kap, b1, b0, discU, resBU, discF];
end function;

function IsOpenQ(vals)
    return &and[b ne 0 : b in BoundaryValuesQ(vals)];
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
    r := a mod m;
    if 2*r gt m then
        return r - m;
    end if;
    return r;
end function;

function BalancedVector(vals, m)
    return [BalancedResidue(Integers()!v, m) : v in vals];
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

function B2ZeroModelQ(vals)
    Q := Rationals();
    P<x> := PolynomialRing(Q);
    h := 1 + x;
    kap := Q!vals[1]; s := Q!vals[2]; t := Q!vals[3];
    b0 := Q!vals[4]; b1 := Q!vals[5];
    U := x^2 + s*x + t;
    B := b0 + b1*x;
    A, E, f, avec := ForcedAData(h, kap, U, B);

    if &or[Coefficient(E, i) ne 0 : i in [0..4]] then
        return false, A, B, f, U, P!0, [], "residual nonzero";
    end if;
    if not IsOpenQ(vals) then
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

function LiftChart(label, p, base5, max_power, resZ)
    printf "\n# lift_chart %o p=%o max_power=%o\n", label, p, max_power;
    printf "  base_mod_%o=<K=%o,s=%o,t=%o,b0=%o,b1=%o>\n",
        p, base5[1], base5[2], base5[3], base5[4], base5[5];
    printf "  base_boundary_values=%o\n", [Integers()!b mod p : b in BoundaryValuesQ(base5)];

    rank, det := ResidualJacobianRankModP(resZ, base5, p);
    printf "  slice_jacobian_rank_mod_%o=%o determinant=%o\n", p, rank, det;

    sols := [base5];
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

function RationalCandidatesForResidue(residue, modulus, num_bound, den_bound)
    Q := Rationals();
    Z := Integers();
    Zm := Integers(modulus);
    want := Zm!residue;
    cands := [];
    seen := {};
    for den in [1..den_bound] do
        if GCD(den, modulus) ne 1 then
            continue;
        end if;
        den_mod := Zm!den;
        for num in [-num_bound..num_bound] do
            if GCD(Abs(num), den) ne 1 then
                continue;
            end if;
            if (Zm!num)/den_mod eq want then
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

procedure RationalSearchChart(label, modulus, residues, resZ)
    printf "\n# rational_search %o modulus=%o num_bound=%o den_bound=%o\n",
        label, modulus, rational_num_bound, rational_den_bound;
    cand_lists := [
        RationalCandidatesForResidue(residues[i], modulus, rational_num_bound,
            rational_den_bound) : i in [1..5]
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
        printf "  candidates_%o_preview=%o\n", i, preview;
    end for;

    if prod eq 0 then
        print "  search_result=no_candidates";
        return;
    end if;
    if prod gt max_search_tuples then
        printf "  search_result=skipped_tuple_cap max_search_tuples=%o\n",
            max_search_tuples;
        return;
    end if;

    hits := [];
    tested := 0;
    for kapv in cand_lists[1] do
    for sv in cand_lists[2] do
    for tv in cand_lists[3] do
    for b0v in cand_lists[4] do
    for b1v in cand_lists[5] do
        vals := [kapv, sv, tv, b0v, b1v];
        tested +:= 1;
        if ResidualsZeroExact(resZ, vals) then
            Append(~hits, vals);
            ok, A, B, f, U, V, rels, reason := B2ZeroModelQ(vals);
            printf "  exact_hit_%o vals=<K=%o,s=%o,t=%o,b0=%o,b1=%o> class=%o\n",
                #hits, kapv, sv, tv, b0v, b1v, reason;
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
    printf "  tested=%o exact_hits=%o\n", tested, #hits;
    if #hits eq 0 then
        print "  search_result=no_exact_rational_candidate";
    end if;
end procedure;

function CRTPair(a, m, b, n)
    Zn := Integers(n);
    delta := Integers()!((Zn!(b - a))/(Zn!m));
    return (a + m*delta) mod (m*n);
end function;

function CRTVector(a, m, b, n)
    return [CRTPair(a[i], m, b[i], n) : i in [1..#a]];
end function;

R, resQ, Atemplate, ftemplate, Utemplate, Btemplate := B2ZeroSliceResidualsQ();
resZ, dens := IntegralResiduals(resQ);

print "# Z/5 x Z/5 b2=0 h1=1 h2=0 lift/reconstruction driver";
printf "parameters max_7_power=%o max_11_power=%o rational_num_bound=%o rational_den_bound=%o max_search_tuples=%o\n",
    max_7_power, max_11_power, rational_num_bound, rational_den_bound,
    max_search_tuples;
printf "integral_residual_denominators=%o\n", dens;
for i in [1..#resZ] do
    printf "  E%o degree=%o total_degree=%o terms=%o\n",
        i - 1, Degree(resZ[i]), TotalDegree(resZ[i]), #Terms(resZ[i]);
end for;

// Open slice points supplied by agent_z5x5_b2zero_elim.
F7A := [2,0,6,3,1];
F7B := [2,0,6,4,6];
F11A := [7,5,10,2,4];

sols7A, mod7A := LiftChart("F7_slice_A", 7, F7A, max_7_power, resZ);
sols7B, mod7B := LiftChart("F7_slice_B", 7, F7B, max_7_power, resZ);
sols11A, mod11A := LiftChart("F11_slice_A", 11, F11A, max_11_power, resZ);

if #sols7A gt 0 then
    RationalSearchChart("F7_slice_A", mod7A, sols7A[1], resZ);
end if;
if #sols7B gt 0 then
    RationalSearchChart("F7_slice_B", mod7B, sols7B[1], resZ);
end if;
if #sols11A gt 0 then
    RationalSearchChart("F11_slice_A", mod11A, sols11A[1], resZ);
end if;

if do_crt_search and #sols7A gt 0 and #sols11A gt 0 then
    crtA := CRTVector(sols7A[1], mod7A, sols11A[1], mod11A);
    RationalSearchChart("CRT_F7A_F11A", mod7A*mod11A, crtA, resZ);
end if;

if do_crt_search and #sols7B gt 0 and #sols11A gt 0 then
    crtB := CRTVector(sols7B[1], mod7B, sols11A[1], mod11A);
    RationalSearchChart("CRT_F7B_F11A", mod7B*mod11A, crtB, resZ);
end if;

print "\n# verdict_marker";
print "smooth_lift_targets_have_unique_lifts_when_slice_jacobian_rank_is_5";
print "bounded_rational_search_completed";

quit;
