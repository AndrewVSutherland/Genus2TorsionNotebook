/*
Close the non-quartic pieces in the four best A(2,12) -> A(2,24) fibers.

This extends agent_A2_24_quartic_extract.m by inspecting every translated
order-12 class in the common split fiber, including the residual degree-8
extra factors and the O/TR degree-16 factors.

Run:
    magma code/agent_A2_24_branch_closure.m

Optional:
    magma -b SearchHeight:=8 PrintPolynomials:=false \
        code/agent_A2_24_branch_closure.m
*/

SetColumns(0);

QQ := Rationals();
ZZ := Integers();
P<x> := PolynomialRing(QQ);

if assigned SearchHeight and Type(SearchHeight) eq MonStgElt then
    SearchHeight := StringToInteger(SearchHeight);
end if;
if not assigned SearchHeight then
    SearchHeight := 8;
end if;

if assigned PrintPolynomials and Type(PrintPolynomials) eq MonStgElt then
    PrintPolynomials := PrintPolynomials in {"true", "True", "1", "yes", "Yes"};
end if;
if not assigned PrintPolynomials then
    PrintPolynomials := true;
end if;

if assigned BoundaryOverQuadratic and Type(BoundaryOverQuadratic) eq MonStgElt then
    BoundaryOverQuadratic := BoundaryOverQuadratic in {"true", "True", "1", "yes", "Yes"};
end if;
if not assigned BoundaryOverQuadratic then
    BoundaryOverQuadratic := true;
end if;

T1u := x^2 + x + QQ!1/QQ!3;
T2u := x^2 + QQ!8/QQ!3*x + QQ!16/QQ!3;
TRu := x^2 + x + 7;

DminusU := x^2 - QQ!2/QQ!3*x - QQ!4/QQ!3;
DminusV := -QQ!7/QQ!3*x - QQ!8/QQ!3;
DplusU := x^2 + QQ!13/QQ!3*x + QQ!26/QQ!3;
DplusV := QQ!5/QQ!3*x + QQ!16/QQ!3;

function SumInts(xs)
    if #xs eq 0 then
        return 0;
    end if;
    return &+xs;
end function;

function MakeMonic(g)
    if g eq 0 then
        return g;
    end if;
    return g/LeadingCoefficient(g);
end function;

function ContainsPoint(seq, T)
    for S in seq do
        if S eq T then
            return true;
        end if;
    end for;
    return false;
end function;

function HeightRationals(H)
    vals := [QQ!0];
    if H le 0 then
        return vals;
    end if;
    for den in [1..H] do
        for num in [-H..H] do
            if GCD(num, den) eq 1 then
                Append(~vals, QQ!num/QQ!den);
            end if;
        end for;
    end for;
    return Sort(Setseq(Seqset(vals)));
end function;

function FactorDegreeMults(g)
    if g eq 0 or Degree(g) lt 1 then
        return [];
    end if;
    return [<Degree(fe[1]), fe[2]> : fe in Factorization(g)];
end function;

function FactorDegreeMultsInVariable(facs, var)
    return [<Degree(fe[1], var), fe[2]> : fe in facs];
end function;

function A12Data(p_in, z_in, r_in)
    p := QQ!p_in;
    z := QQ!z_in;
    r := QQ!r_in;

    if p eq 0 or z eq 0 then
        return false, P!0, P!0, P!0, P!0, P!0;
    end if;

    s := (z^2 - 4*p^2 + 1)/(2*z);
    if s^2 eq 1 then
        return false, P!0, P!0, P!0, P!0, P!0;
    end if;

    t := (z^2 + 4*p^2 - 1)^2/(8*p^2*z);
    mu := ((s^2 - 1)*(2*p*r + 1) - p^2*(2*s*t - 4))/(4*p^3);
    lambda := (4 - mu^2)*p^2/(s^2 - 1);
    if lambda eq 0 then
        return false, P!0, P!0, P!0, P!0, P!0;
    end if;

    T := p*x + r;
    Rpol := (T^2 + x - 1)/lambda;
    ell := s*x + t;
    Qpol := 2*T + mu*Rpol;
    F := Rpol*x^2 + 4*(Rpol + x - 1)*(Rpol - 1);
    if F ne Qpol^2 + Rpol*ell^2 then
        return false, P!0, P!0, P!0, P!0, P!0;
    end if;

    return true, Rpol, Qpol, ell, F, Rpol*F;
end function;

function A12JacobianData(Rpol, Qpol, ell, f)
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    O := J!0;

    TR := J![MakeMonic(Rpol), P!0];
    u4 := MakeMonic(Qpol);
    v4 := (Rpol*ell) mod u4;
    P4 := J![u4, v4];

    u6 := MakeMonic(Rpol + x - 1);
    v6 := (x*Rpol) mod u6;
    P6 := J![u6, v6];

    return J, O, TR, P4 + P6;
end function;

function RationalTwoTorsionClasses(J, f)
    O := J!0;
    fac := Factorization(f);
    atoms := [];
    for fe in fac do
        for i in [1..fe[2]] do
            Append(~atoms, MakeMonic(fe[1]));
        end for;
    end for;

    gens := [];
    if #atoms gt 0 then
        for mask in [1..2^#atoms - 1] do
            g := Parent(f)!1;
            for i in [1..#atoms] do
                if ((mask div 2^(i - 1)) mod 2) eq 1 then
                    g *:= atoms[i];
                end if;
            end for;
            if Degree(g) gt 2 then
                continue;
            end if;
            g := MakeMonic(g);
            try
                T := J![g, P!0];
                if T ne O and 2*T eq O and not ContainsPoint(gens, T) then
                    Append(~gens, T);
                end if;
            catch e
                continue;
            end try;
        end for;
    end if;

    twos := [O];
    changed := true;
    while changed do
        changed := false;
        snapshot := twos;
        for A in snapshot do
            for G in gens do
                S := A + G;
                if not ContainsPoint(twos, S) then
                    Append(~twos, S);
                    changed := true;
                end if;
            end for;
        end for;
    end while;
    return twos;
end function;

function TwoLabel(T, O, TR)
    if T eq O then
        return "O";
    end if;
    if T eq TR then
        return "TR";
    end if;
    uv := Eltseq(T);
    if #uv ge 2 and uv[2] eq 0 then
        if uv[1] eq T1u then
            return "T1_extra";
        end if;
        if uv[1] eq T2u then
            return "T2_extra";
        end if;
    end if;
    return "extra_unknown";
end function;

function DLabel(D)
    uv := Eltseq(D);
    if #uv lt 2 then
        return "non_mumford";
    end if;
    if uv[1] eq DminusU and uv[2] eq DminusV then
        return "D_minus";
    end if;
    if uv[1] eq DplusU and uv[2] eq DplusV then
        return "D_plus";
    end if;
    return "O_or_TR_class";
end function;

function HalvingEquations(u, v, f)
    A<M,N> := PolynomialRing(QQ, 2);
    RX<X> := PolynomialRing(A);
    phi := hom<P -> RX | X>;

    uX := phi(u);
    vX := phi(v);
    fX := phi(f);
    ell := vX + uX*(M*X + N);

    if (ell^2 - fX) mod uX ne 0 then
        return false, A!0, A!0, RX!0;
    end if;

    S := ExactQuotient(ell^2 - fX, uX);
    if Degree(S) ne 4 then
        return false, A!0, A!0, S;
    end if;

    s4 := Coefficient(S, 4);
    s3 := Coefficient(S, 3);
    s2 := Coefficient(S, 2);
    s1 := Coefficient(S, 1);
    s0 := Coefficient(S, 0);

    E1 := 8*s4^2*s1 - s3*(4*s4*s2 - s3^2);
    E0 := 64*s4^3*s0 - (4*s4*s2 - s3^2)^2;

    return true, E1, E0, S;
end function;

function ToUnivariateInM(qA)
    A := Parent(qA);
    PM<m> := PolynomialRing(QQ);
    phi := hom<A -> PM | [m, PM!0]>;
    return phi(qA);
end function;

function RationalSolutionsAtM(E1, E0, S, Mv)
    A := Parent(E1);
    PN<n> := PolynomialRing(QQ);
    psi := hom<A -> PN | [PN!Mv, n]>;
    g1 := psi(E1);
    g0 := psi(E0);
    g := GCD(g1, g0);
    pts := [];

    if g eq PN!0 then
        return pts, "positive_dimensional_N_gcd";
    end if;
    if Degree(g) lt 1 then
        return pts, "no_N_root";
    end if;

    g := MakeMonic(g);
    s4 := Coefficient(S, 4);
    for rt in Roots(g) do
        Nv := rt[1];
        if Evaluate(E1, [Mv, Nv]) eq 0 and Evaluate(E0, [Mv, Nv]) eq 0 and
                Evaluate(s4, [Mv, Nv]) ne 0 then
            Append(~pts, <Mv, Nv>);
        end if;
    end for;
    return pts, "finite_N_gcd";
end function;

function DirectRationalSolutions(E1, E0, S, H)
    if H le 0 then
        return [];
    end if;
    s4 := Coefficient(S, 4);
    vals := HeightRationals(H);
    pts := [];
    for Mv in vals do
        for Nv in vals do
            if Evaluate(E1, [Mv, Nv]) eq 0 and Evaluate(E0, [Mv, Nv]) eq 0 and
                    Evaluate(s4, [Mv, Nv]) ne 0 then
                Append(~pts, <Mv, Nv>);
            end if;
        end for;
    end for;
    return pts;
end function;

function SquareIntegralPolynomial(g)
    den := 1;
    for c in Coefficients(g) do
        den := LCM(den, Denominator(c));
    end for;
    return Parent(g)!(den^2*g);
end function;

function SquareModelScale(g)
    den := 1;
    for c in Coefficients(g) do
        den := LCM(den, Denominator(c));
    end for;
    return den;
end function;

function IntegralSquareModelDivisor(f, D)
    uv := Eltseq(D);
    if #uv lt 2 then
        return false, Jacobian(HyperellipticCurve(SquareIntegralPolynomial(f)))!0,
            SquareIntegralPolynomial(f), 1;
    end if;

    den := SquareModelScale(f);
    fI := Parent(f)!(den^2*f);
    JI := Jacobian(HyperellipticCurve(fI));
    vI := (den*uv[2]) mod uv[1];
    return true, JI![uv[1], vI], fI, den;
end function;

function TorsionOrderFromInvariants(invs)
    if #invs eq 0 then
        return 1;
    end if;
    return &*invs;
end function;

function ContainsSubgroup(invG, invH)
    all := invG cat invH;
    primes := Sort(Setseq(Seqset(&cat[PrimeDivisors(n) : n in all | n ne 0])));
    for p in primes do
        emax := Max([Valuation(n, p) : n in invH]);
        for e in [1..emax] do
            if #[n : n in invG | Valuation(n, p) ge e] lt
               #[n : n in invH | Valuation(n, p) ge e] then
                return false;
            end if;
        end for;
    end for;
    return true;
end function;

procedure CertifyCandidate(fiberIndex, p0, z0, r0, twoLabel, dLabel, D12, J, f, pts)
    if #pts eq 0 then
        return;
    end if;

    printf "    RATIONAL_HALVING_CANDIDATES fiber=%o p=%o z=%o r=%o two_label=%o D_label=%o points=%o\n",
        fiberIndex, p0, z0, r0, twoLabel, dLabel, pts;

    okI, DI, fI, denI := IntegralSquareModelDivisor(f, D12);
    if okI then
        divisible := IsDivisibleBy(DI, 2);
        if divisible then
            _, half := IsDivisibleBy(DI, 2);
            printf "    EXACT_HALVING_CERT integral_den=%o half_order=%o twice_half_matches=%o\n",
                denI, Order(half), 2*half eq DI;
        else
            printf "    EXACT_HALVING_CERT integral_den=%o IsDivisibleBy(D,2)=false\n", denI;
        end if;
    else
        printf "    EXACT_HALVING_CERT skipped_integral_model_conversion_failed\n";
    end if;

    try
        TI, mp := TorsionSubgroup(Jacobian(HyperellipticCurve(fI)));
        invs := Invariants(TI);
        printf "    EXACT_TORSION integral_square_model_invariants=%o order=%o contains_Z2xZ24=%o fI=%o\n",
            invs, TorsionOrderFromInvariants(invs), ContainsSubgroup(invs, [2, 24]), fI;
    catch e
        printf "    EXACT_TORSION_FAILED error=%o\n", e`Object;
    end try;
end procedure;

procedure AnalyzeBoundary(fiberIndex, p0, z0, r0, twoLabel, dLabel, E1, E0, S,
        ~boundaryRationalClasses)
    A := Parent(E1);
    s4 := Coefficient(S, 4);
    qS4 := MakeMonic(ToUnivariateInM(s4));
    rootsQ := Roots(qS4);
    rationalBoundaryPts := [];

    for rt in rootsQ do
        Mv := rt[1];
        ptsM, status := RationalSolutionsAtM(E1, E0, S, Mv);
        rationalBoundaryPts cat:= ptsM;
    end for;

    if #rationalBoundaryPts gt 0 then
        boundaryRationalClasses +:= 1;
    end if;

    printf "    BOUNDARY_S4 fiber=%o p=%o z=%o r=%o two_label=%o D_label=%o s4=%o s4_roots_Q=%o rational_boundary_points=%o legitimate_halving=false reason=%o\n",
        fiberIndex, p0, z0, r0, twoLabel, dLabel, s4, rootsQ,
        rationalBoundaryPts, (#rootsQ eq 0 select "no_rational_M_on_s4=0" else "s4_zero_degenerate_square_quartic");

    if BoundaryOverQuadratic and Degree(qS4) eq 2 and #rootsQ eq 0 then
        try
            K<a> := NumberField(qS4);
            PN<n> := PolynomialRing(K);
            psi := hom<A -> PN | [PN!a, n]>;
            gN := GCD(psi(E1), psi(E0));
            gNdeg := (gN eq PN!0) select -1 else Degree(gN);
            specDeg := -1;
            for j in [0..Degree(S)] do
                if psi(Coefficient(S, j)) ne PN!0 then
                    specDeg := j;
                end if;
            end for;
            printf "      boundary_over_Qsqrt5 gcd_N_degree=%o specialized_S_degree=%o classification=discarded_boundary_fake\n",
                gNdeg, specDeg;
        catch e
            printf "      boundary_over_quadratic_failed error=%o\n", e`Object;
        end try;
    end if;
end procedure;

procedure AnalyzeAffineFactor(fiberIndex, p0, z0, r0, twoLabel, dLabel,
        factorIndex, qA, mult, E1, E0, S, D12, J, f,
        ~totalAffineFactors, ~affineFactorsWithQRoot, ~rationalPointClasses)
    totalAffineFactors +:= 1;

    qM := MakeMonic(ToUnivariateInM(qA));
    qFac := Factorization(qM);
    roots := Roots(qM);
    qFacDegrees := [<Degree(fe[1]), fe[2]> : fe in qFac];
    points := [];
    statuses := [];

    for rt in roots do
        ptsM, status := RationalSolutionsAtM(E1, E0, S, rt[1]);
        points cat:= ptsM;
        Append(~statuses, <rt[1], status, ptsM>);
    end for;

    if #roots gt 0 then
        affineFactorsWithQRoot +:= 1;
    end if;
    if #points gt 0 then
        rationalPointClasses +:= 1;
    end if;

    printf "    AFFINE_FACTOR fiber=%o p=%o z=%o r=%o two_label=%o D_label=%o factor=%o degree_M=%o mult=%o factorization_degrees_Q=%o rational_M_roots=%o rational_points=%o\n",
        fiberIndex, p0, z0, r0, twoLabel, dLabel, factorIndex,
        Degree(qM), mult, qFacDegrees, roots, points;
    if #roots gt 0 then
        printf "      rational_M_root_lifts=%o\n", statuses;
    end if;
    if PrintPolynomials then
        printf "      qM_monic=%o\n", qM;
        printf "      qM_factorization_Q=%o\n", qFac;
    end if;

    CertifyCandidate(fiberIndex, p0, z0, r0, twoLabel, dLabel, D12, J, f, points);
end procedure;

procedure AnalyzeTranslatedClass(fiberIndex, p0, z0, r0, twoLabel, twoIndex,
        twoCount, T, D12, J, O, f,
        ~classRows, ~affineClasses, ~totalAffineFactors,
        ~affineFactorsWithQRoot, ~rationalPointClasses,
        ~boundaryRationalClasses, ~exactDivisibleClasses)
    uvT := Eltseq(T);
    uvD := Eltseq(D12);
    dLabel := DLabel(D12);

    classRows +:= 1;
    printf "  TRANSLATED_CLASS fiber=%o p=%o z=%o r=%o two_label=%o two_index=%o two_count=%o D_label=%o\n",
        fiberIndex, p0, z0, r0, twoLabel, twoIndex, twoCount, dLabel;
    if #uvT ge 2 then
        printf "    T_u=%o\n", uvT[1];
        printf "    T_v=%o\n", uvT[2];
    else
        printf "    T_is_origin=true\n";
    end if;
    if #uvD lt 2 or Degree(uvD[1]) ne 2 then
        printf "    CLASS_SKIP reason=non_degree2_translated_class D_data=%o\n", uvD;
        return;
    end if;

    printf "    D_u=%o\n", uvD[1];
    printf "    D_v=%o\n", uvD[2];
    printf "    D_exact_order12=%o\n",
        12*D12 eq O and &and[n*D12 ne O : n in [1..11]];

    try
        okI, DI, fI, denI := IntegralSquareModelDivisor(f, D12);
        if okI then
            divisible := IsDivisibleBy(DI, 2);
            printf "    exact_integral_square_model_den=%o IsDivisibleBy_D_by_2=%o\n",
                denI, divisible;
            if divisible then
                exactDivisibleClasses +:= 1;
                _, half := IsDivisibleBy(DI, 2);
                printf "    exact_half_order=%o twice_half_matches=%o\n", Order(half), 2*half eq DI;
            end if;
        else
            printf "    exact_integral_square_model_conversion_failed\n";
        end if;
    catch e
        printf "    exact_integral_IsDivisibleBy_failed error=%o\n", e`Object;
    end try;

    ok, E1, E0, S := HalvingEquations(uvD[1], uvD[2], f);
    if not ok then
        printf "    CLASS_ERROR reason=halving_equations_failed\n";
        return;
    end if;

    A := Parent(E1);
    M := A.1;
    N := A.2;
    s4 := Coefficient(S, 4);
    RN := Resultant(E1, E0, N);
    RNfac := Factorization(RN);
    affineFacs := [<fe[1], fe[2]> : fe in RNfac | GCD(fe[1], s4) eq A!1];
    boundaryFacs := [<fe[1], fe[2]> : fe in RNfac | GCD(fe[1], s4) ne A!1];
    affineDegree := SumInts([Degree(fe[1], M)*fe[2] : fe in affineFacs]);
    boundaryDegree := SumInts([Degree(fe[1], M)*fe[2] : fe in boundaryFacs]);
    affineClasses +:= 1;

    printf "    E_degrees E1=(total %o,M %o,N %o) E0=(total %o,M %o,N %o) gcd_E_degree=%o\n",
        TotalDegree(E1), Degree(E1, M), Degree(E1, N),
        TotalDegree(E0), Degree(E0, M), Degree(E0, N), TotalDegree(GCD(E1, E0));
    printf "    RESULTANT_SUMMARY raw_degree_M=%o boundary_degree=%o boundary_factor_degrees=%o saturated_degree=%o affine_factor_degrees=%o\n",
        Degree(RN, M), boundaryDegree, FactorDegreeMultsInVariable(boundaryFacs, M),
        affineDegree, FactorDegreeMultsInVariable(affineFacs, M);

    AnalyzeBoundary(fiberIndex, p0, z0, r0, twoLabel, dLabel, E1, E0, S,
        ~boundaryRationalClasses);

    pts := DirectRationalSolutions(E1, E0, S, SearchHeight);
    printf "    direct_Q_search_height=%o saturated_points=%o points=%o\n",
        SearchHeight, #pts, pts;
    if #pts gt 0 then
        rationalPointClasses +:= 1;
        CertifyCandidate(fiberIndex, p0, z0, r0, twoLabel, dLabel, D12, J, f, pts);
    end if;

    factorIndex := 0;
    for fe in affineFacs do
        factorIndex +:= 1;
        AnalyzeAffineFactor(fiberIndex, p0, z0, r0, twoLabel, dLabel,
            factorIndex, fe[1], fe[2], E1, E0, S, D12, J, f,
            ~totalAffineFactors, ~affineFactorsWithQRoot, ~rationalPointClasses);
    end for;

    if #affineFacs eq 0 then
        printf "    PROJECTION_OBSTRUCTION status=no_affine_factor_after_saturation\n";
    elif &and[#Roots(MakeMonic(ToUnivariateInM(fe[1]))) eq 0 : fe in affineFacs] then
        printf "    PROJECTION_OBSTRUCTION status=no_rational_M_in_any_saturated_factor\n";
    else
        printf "    PROJECTION_OBSTRUCTION status=rational_M_seen_needs_lift_check\n";
    end if;
end procedure;

procedure AnalyzeFiber(fiberIndex, p0, z0, r0,
        ~totalFibers, ~classRows, ~affineClasses, ~totalAffineFactors,
        ~affineFactorsWithQRoot, ~rationalPointClasses,
        ~boundaryRationalClasses, ~exactDivisibleClasses)
    ok, Rpol, Qpol, ell, F, f := A12Data(p0, z0, r0);
    if not ok or Degree(f) ne 6 or Discriminant(f) eq 0 then
        printf "FIBER_ERROR index=%o p=%o z=%o r=%o reason=bad_A12_or_singular\n",
            fiberIndex, p0, z0, r0;
        return;
    end if;

    totalFibers +:= 1;
    printf "\nFIBER index=%o p=%o z=%o r=%o\n", fiberIndex, p0, z0, r0;
    printf "  R=%o\n", Rpol;
    printf "  F=%o\n", F;
    printf "  f=%o\n", f;
    printf "  F_factorization=%o F_factor_degrees=%o f_factor_degrees=%o\n",
        Factorization(F), FactorDegreeMults(F), FactorDegreeMults(f);

    J, O, TR, P12 := A12JacobianData(Rpol, Qpol, ell, f);
    uvP12 := Eltseq(P12);
    uvTR := Eltseq(TR);
    printf "  P12_exact_order12=%o\n",
        12*P12 eq O and &and[n*P12 ne O : n in [1..11]];
    printf "  P12_u=%o\n", uvP12[1];
    printf "  P12_v=%o\n", uvP12[2];
    printf "  cyclic_two_torsion_TR_u=%o TR_v=%o\n", uvTR[1], uvTR[2];

    twos := RationalTwoTorsionClasses(J, f);
    printf "  rational_two_torsion_count=%o\n", #twos;
    for i in [1..#twos] do
        D12 := P12 + twos[i];
        if not (12*D12 eq O and &and[n*D12 ne O : n in [1..11]]) then
            printf "  TWO_CLASS_SKIP index=%o label=%o reason=translated_not_exact_order12\n",
                i, TwoLabel(twos[i], O, TR);
            continue;
        end if;
        AnalyzeTranslatedClass(fiberIndex, p0, z0, r0, TwoLabel(twos[i], O, TR),
            i, #twos, twos[i], D12, J, O, f,
            ~classRows, ~affineClasses, ~totalAffineFactors,
            ~affineFactorsWithQRoot, ~rationalPointClasses,
            ~boundaryRationalClasses, ~exactDivisibleClasses);
    end for;
end procedure;

print "A2_24_BRANCH_CLOSURE_START";
printf "parameters SearchHeight=%o PrintPolynomials=%o BoundaryOverQuadratic=%o\n",
    SearchHeight, PrintPolynomials, BoundaryOverQuadratic;

targets := [
    <QQ!(-1)/QQ!3, QQ!(-1), QQ!4/QQ!3>,
    <QQ!(-1)/QQ!3, QQ!1, QQ!4/QQ!3>,
    <QQ!1/QQ!3, QQ!(-1), QQ!(-4)/QQ!3>,
    <QQ!1/QQ!3, QQ!1, QQ!(-4)/QQ!3>
];

totalFibers := 0;
classRows := 0;
affineClasses := 0;
totalAffineFactors := 0;
affineFactorsWithQRoot := 0;
rationalPointClasses := 0;
boundaryRationalClasses := 0;
exactDivisibleClasses := 0;

for i in [1..#targets] do
    AnalyzeFiber(i, targets[i][1], targets[i][2], targets[i][3],
        ~totalFibers, ~classRows, ~affineClasses, ~totalAffineFactors,
        ~affineFactorsWithQRoot, ~rationalPointClasses,
        ~boundaryRationalClasses, ~exactDivisibleClasses);
end for;

printf "\nA2_24_BRANCH_CLOSURE_DONE fibers=%o translated_classes=%o affine_classes=%o affine_factors=%o affine_factors_with_Q_root=%o rational_point_classes=%o boundary_rational_classes=%o exact_divisible_classes=%o cputime=%o\n",
    totalFibers, classRows, affineClasses, totalAffineFactors,
    affineFactorsWithQRoot, rationalPointClasses, boundaryRationalClasses,
    exactDivisibleClasses, Cputime();
