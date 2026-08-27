/*
Extract and analyze the quartic saturated halving branches in the best
A(2,12) -> A(2,24) split fibers from the third pass.

Run:
    magma code/agent_A2_24_quartic_extract.m

Optional:
    magma -b SearchHeight:=10 \
        code/agent_A2_24_quartic_extract.m
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

function SumInts(xs)
    if #xs eq 0 then
        return 0;
    end if;
    return &+xs;
end function;

function MakeMonic(g)
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
    return "extra";
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

function QuarticResolvent(q)
    PM := Parent(q);
    K := BaseRing(PM);
    R<eta> := PolynomialRing(K);

    qm := MakeMonic(q);
    B := Coefficient(qm, 3);
    C := Coefficient(qm, 2);
    D := Coefficient(qm, 1);
    E := Coefficient(qm, 0);

    return eta^3 - C*eta^2 + (B*D - 4*E)*eta
        + (4*C*E - B^2*E - D^2);
end function;

function DirectRationalSolutions(E1, E0, S, H)
    A := Parent(E1);
    s4 := Coefficient(S, 4);
    vals := HeightRationals(H);
    pts := [];
    for Mv in vals do
        for Nv in vals do
            ev := [Mv, Nv];
            if Evaluate(E1, ev) eq 0 and Evaluate(E0, ev) eq 0
                    and Evaluate(s4, ev) ne 0 then
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

procedure CertifyCandidate(fiberIndex, p0, z0, r0, twoIndex, D12, J, f, pts)
    if #pts eq 0 then
        return;
    end if;

    printf "  RATIONAL_HALVING_CANDIDATES fiber=%o p=%o z=%o r=%o two_index=%o points=%o\n",
        fiberIndex, p0, z0, r0, twoIndex, pts;

    divisible := IsDivisibleBy(D12, 2);
    if divisible then
        _, half := IsDivisibleBy(D12, 2);
        printf "  EXACT_HALVING_CERT half_order=%o twice_half_matches=%o\n",
            Order(half), 2*half eq D12;
    else
        printf "  EXACT_HALVING_CERT failed IsDivisibleBy(D12,2)=false\n";
    end if;

    try
        fI := SquareIntegralPolynomial(f);
        TI, mp := TorsionSubgroup(Jacobian(HyperellipticCurve(fI)));
        invs := Invariants(TI);
        printf "  EXACT_TORSION integral_square_model_invariants=%o order=%o contains_Z2xZ24=%o fI=%o\n",
            invs, TorsionOrderFromInvariants(invs), ContainsSubgroup(invs, [2, 24]), fI;
    catch e
        printf "  EXACT_TORSION_FAILED error=%o\n", e`Object;
    end try;
end procedure;

procedure AnalyzeQuarticFactor(fiberIndex, p0, z0, r0, twoIndex, quarticIndex,
        E1, E0, S, qA, mult, ~totalQuartics, ~rationalQuartics)
    totalQuartics +:= 1;

    qM := MakeMonic(ToUnivariateInM(qA));
    qFac := Factorization(qM);
    roots := Roots(qM);
    if #roots gt 0 then
        rationalQuartics +:= 1;
    end if;

    disc := Discriminant(qM);
    discNum := ZZ!Numerator(disc);
    discDen := ZZ!Denominator(disc);
    discSquare := IsSquare(disc);
    resolvent := QuarticResolvent(qM);

    printf "    QUARTIC_BRANCH fiber=%o p=%o z=%o r=%o two_index=%o branch=%o mult=%o\n",
        fiberIndex, p0, z0, r0, twoIndex, quarticIndex, mult;
    printf "      qM_monic=%o\n", qM;
    printf "      q_factorization_Q=%o rational_roots=%o\n", qFac, roots;
    printf "      discriminant=%o disc_is_square=%o disc_num_factorization=%o disc_den_factorization=%o\n",
        disc, discSquare, Factorization(Abs(discNum)), Factorization(discDen);
    printf "      cubic_resolvent=%o resolvent_factorization=%o resolvent_roots_Q=%o\n",
        resolvent, Factorization(resolvent), Roots(resolvent);

    if Degree(qM) eq 4 and #qFac eq 1 and qFac[1][2] eq 1 then
        try
            K<a> := NumberField(qM);
            PN<n> := PolynomialRing(K);
            A := Parent(E1);
            psi := hom<A -> PN | [PN!a, n]>;
            gN := GCD(psi(E1), psi(E0));
            if Degree(gN) ge 1 then
                gN := gN/LeadingCoefficient(gN);
            end if;
            s4N := psi(Coefficient(S, 4));
            printf "      branch_field=Q<a>/(qM) gcd_N=%o degree_N=%o\n",
                gN, Degree(gN);
            if Degree(gN) eq 1 then
                Nexpr := -Coefficient(gN, 0);
                printf "      N_expression=%o\n", Nexpr;
                printf "      s4_at_branch=%o nonzero=%o\n",
                    Evaluate(s4N, Nexpr), Evaluate(s4N, Nexpr) ne 0;
                printf "      component_type=zero_dimensional_linear_N no_low_genus_model\n";
            else
                printf "      component_type=nonlinear_N_gcd needs_plane_curve_followup\n";
            end if;
        catch e
            printf "      branch_field_analysis_failed=%o\n", e`Object;
        end try;
    else
        printf "      branch_field_analysis_skipped reducible_or_nonquartic\n";
    end if;
end procedure;

procedure AnalyzeExtraClass(fiberIndex, p0, z0, r0, twoIndex, twoCount,
        T, D12, J, O, f, ~totalQuartics, ~rationalQuartics, ~rationalPts)
    uvT := Eltseq(T);
    uvD := Eltseq(D12);
    if #uvD lt 2 or Degree(uvD[1]) ne 2 then
        printf "  EXTRA_SKIP two_index=%o reason=non_degree2_D u_data=%o\n",
            twoIndex, uvD;
        return;
    end if;

    ok, E1, E0, S := HalvingEquations(uvD[1], uvD[2], f);
    if not ok then
        printf "  EXTRA_ERROR two_index=%o reason=halving_equations_failed\n", twoIndex;
        return;
    end if;

    A := Parent(E1);
    M := A.1;
    N := A.2;
    s4 := Coefficient(S, 4);
    RN := Resultant(E1, E0, N);
    RNfac := Factorization(RN);
    affineFacs := [<fe[1], fe[2]> : fe in RNfac | GCD(fe[1], s4) eq 1];
    boundaryFacs := [<fe[1], fe[2]> : fe in RNfac | GCD(fe[1], s4) ne 1];
    affineDegrees := [<Degree(fe[1], M), fe[2]> : fe in affineFacs];
    boundaryDegrees := [<Degree(fe[1], M), fe[2]> : fe in boundaryFacs];
    affineDegree := SumInts([Degree(fe[1], M)*fe[2] : fe in affineFacs]);

    printf "  EXTRA_TRANSLATION fiber=%o p=%o z=%o r=%o two_index=%o two_count=%o\n",
        fiberIndex, p0, z0, r0, twoIndex, twoCount;
    printf "    T_u=%o\n", uvT[1];
    printf "    T_v=%o\n", uvT[2];
    printf "    D=P12+T u=%o\n", uvD[1];
    printf "    D=P12+T v=%o\n", uvD[2];
    printf "    D_exact_order12=%o\n",
        12*D12 eq O and &and[n*D12 ne O : n in [1..11]];
    printf "    E_degrees E1=(total %o,M %o,N %o) E0=(total %o,M %o,N %o) gcd_E_degree=%o\n",
        TotalDegree(E1), Degree(E1, M), Degree(E1, N),
        TotalDegree(E0), Degree(E0, M), Degree(E0, N), TotalDegree(GCD(E1, E0));
    printf "    s4=%o\n", s4;
    printf "    s4_factor_degrees=%o raw_resultant_degree_M=%o boundary_factor_degrees=%o saturated_degree=%o affine_factor_degrees=%o\n",
        FactorDegreeMults(s4), Degree(RN, M), boundaryDegrees, affineDegree, affineDegrees;

    pts := DirectRationalSolutions(E1, E0, S, SearchHeight);
    rationalPts +:= #pts;
    printf "    direct_Q_search_height=%o saturated_points=%o points=%o\n",
        SearchHeight, #pts, pts;
    CertifyCandidate(fiberIndex, p0, z0, r0, twoIndex, D12, J, f, pts);

    quarticIndex := 0;
    for fe in affineFacs do
        qA := fe[1];
        if Degree(qA, M) eq 4 then
            quarticIndex +:= 1;
            AnalyzeQuarticFactor(fiberIndex, p0, z0, r0, twoIndex,
                quarticIndex, E1, E0, S, qA, fe[2],
                ~totalQuartics, ~rationalQuartics);
        else
            printf "    NONQUARTIC_AFFINE_FACTOR degree_M=%o mult=%o rational_roots=%o factor=%o\n",
                Degree(qA, M), fe[2], Roots(MakeMonic(ToUnivariateInM(qA))),
                MakeMonic(ToUnivariateInM(qA));
        end if;
    end for;
end procedure;

procedure AnalyzeFiber(fiberIndex, p0, z0, r0,
        ~totalFibers, ~extraClasses, ~totalQuartics, ~rationalQuartics,
        ~rationalPts)
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
    printf "  P12_exact_order12=%o\n",
        12*P12 eq O and &and[n*P12 ne O : n in [1..11]];
    uvTR := Eltseq(TR);
    printf "  cyclic_two_torsion_TR_u=%o TR_v=%o\n", uvTR[1], uvTR[2];

    twos := RationalTwoTorsionClasses(J, f);
    printf "  rational_two_torsion_count=%o\n", #twos;
    for i in [1..#twos] do
        label := TwoLabel(twos[i], O, TR);
        D12 := P12 + twos[i];
        if not (12*D12 eq O and &and[n*D12 ne O : n in [1..11]]) then
            printf "  TWO_CLASS_SKIP index=%o label=%o reason=translated_not_exact_order12\n",
                i, label;
            continue;
        end if;
        if label eq "extra" then
            extraClasses +:= 1;
            AnalyzeExtraClass(fiberIndex, p0, z0, r0, i, #twos,
                twos[i], D12, J, O, f,
                ~totalQuartics, ~rationalQuartics, ~rationalPts);
        else
            okH, E1, E0, S := HalvingEquations(Eltseq(D12)[1], Eltseq(D12)[2], f);
            if okH then
                A := Parent(E1);
                M := A.1;
                N := A.2;
                s4 := Coefficient(S, 4);
                RN := Resultant(E1, E0, N);
                aff := [<Degree(fe[1], M), fe[2]> : fe in Factorization(RN)
                    | GCD(fe[1], s4) eq 1];
                printf "  TWO_CLASS label=%o index=%o affine_factor_degrees=%o\n",
                    label, i, aff;
            end if;
        end if;
    end for;
end procedure;

print "A2_24_QUARTIC_EXTRACT_START";
printf "parameters SearchHeight=%o\n", SearchHeight;

targets := [
    <QQ!(-1)/QQ!3, QQ!(-1), QQ!4/QQ!3>,
    <QQ!(-1)/QQ!3, QQ!1, QQ!4/QQ!3>,
    <QQ!1/QQ!3, QQ!(-1), QQ!(-4)/QQ!3>,
    <QQ!1/QQ!3, QQ!1, QQ!(-4)/QQ!3>
];

totalFibers := 0;
extraClasses := 0;
totalQuartics := 0;
rationalQuartics := 0;
rationalPts := 0;

for i in [1..#targets] do
    AnalyzeFiber(i, targets[i][1], targets[i][2], targets[i][3],
        ~totalFibers, ~extraClasses, ~totalQuartics, ~rationalQuartics,
        ~rationalPts);
end for;

printf "\nA2_24_QUARTIC_EXTRACT_DONE fibers=%o extra_classes=%o quartic_branches=%o quartics_with_Q_root=%o direct_Q_points_height_%o=%o cputime=%o\n",
    totalFibers, extraClasses, totalQuartics, rationalQuartics,
    SearchHeight, rationalPts, Cputime();
