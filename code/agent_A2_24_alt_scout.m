/*
Alternate low-priority scout for saturated A(2,12) -> A(2,24) fibers.

This is a bounded height-6 shell scout beyond the closed four fibers and the
cold height-5 split scan.  It enumerates rational A(12) parameters of bounded
height, optionally keeps only the exact height shell, keeps split A(2,12)
fibers, translates the visible order-12 class by rational 2-torsion, removes
the s4=0 square-quartic boundary, and records the affine M-factor degrees.
Rational M-roots are lifted to N and certified exactly.

Run examples:
    magma code/agent_A2_24_alt_scout.m
    magma -b Height:=6 ShellOnly:=true Progress:=5000 \
        code/agent_A2_24_alt_scout.m

Partition by p-index if needed:
    magma -b Height:=6 ShellOnly:=true PStart:=1 PStop:=12 \
        code/agent_A2_24_alt_scout.m
*/

SetColumns(0);

QQ := Rationals();
ZZ := Integers();
P<x> := PolynomialRing(QQ);

if assigned Height and Type(Height) eq MonStgElt then
    Height := StringToInteger(Height);
end if;
if not assigned Height then
    Height := 6;
end if;

if assigned ShellOnly and Type(ShellOnly) eq MonStgElt then
    ShellOnly := ShellOnly in {"true", "True", "1", "yes", "Yes"};
end if;
if not assigned ShellOnly then
    ShellOnly := true;
end if;

if assigned LowDegree and Type(LowDegree) eq MonStgElt then
    LowDegree := StringToInteger(LowDegree);
end if;
if not assigned LowDegree then
    LowDegree := 4;
end if;

if assigned MaxRows and Type(MaxRows) eq MonStgElt then
    MaxRows := StringToInteger(MaxRows);
end if;
if not assigned MaxRows then
    MaxRows := 0;
end if;

if assigned MaxFibers and Type(MaxFibers) eq MonStgElt then
    MaxFibers := StringToInteger(MaxFibers);
end if;
if not assigned MaxFibers then
    MaxFibers := 0;
end if;

if assigned MaxChecked and Type(MaxChecked) eq MonStgElt then
    MaxChecked := StringToInteger(MaxChecked);
end if;
if not assigned MaxChecked then
    MaxChecked := 0;
end if;

if assigned PStart and Type(PStart) eq MonStgElt then
    PStart := StringToInteger(PStart);
end if;
if not assigned PStart then
    PStart := 1;
end if;

if assigned PStop and Type(PStop) eq MonStgElt then
    PStop := StringToInteger(PStop);
end if;
if not assigned PStop then
    PStop := 0;
end if;

if assigned Progress and Type(Progress) eq MonStgElt then
    Progress := StringToInteger(Progress);
end if;
if not assigned Progress then
    Progress := 5000;
end if;

if assigned SkipClosedBest and Type(SkipClosedBest) eq MonStgElt then
    SkipClosedBest := SkipClosedBest in {"true", "True", "1", "yes", "Yes"};
end if;
if not assigned SkipClosedBest then
    SkipClosedBest := true;
end if;

if assigned VerboseRows and Type(VerboseRows) eq MonStgElt then
    VerboseRows := VerboseRows in {"true", "True", "1", "yes", "Yes"};
end if;
if not assigned VerboseRows then
    VerboseRows := false;
end if;

if assigned DirectSearchHeight and Type(DirectSearchHeight) eq MonStgElt then
    DirectSearchHeight := StringToInteger(DirectSearchHeight);
end if;
if not assigned DirectSearchHeight then
    DirectSearchHeight := 0;
end if;

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

procedure BumpKey(key, ~keys, ~counts)
    for i in [1..#keys] do
        if keys[i] eq key then
            counts[i] +:= 1;
            return;
        end if;
    end for;
    Append(~keys, key);
    Append(~counts, 1);
end procedure;

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

function RatHeight(q)
    q := QQ!q;
    if q eq 0 then
        return 0;
    end if;
    return Max(Abs(Numerator(q)), Denominator(q));
end function;

function TripleHeight(p0, z0, r0)
    return Max([RatHeight(p0), RatHeight(z0), RatHeight(r0)]);
end function;

function CountActiveTriples(vals, pvals, pStart, pStop, H, shellOnly)
    total := 0;
    active := 0;
    for pIndex in [pStart..pStop] do
        p0 := pvals[pIndex];
        for z0 in vals do
            if z0 eq 0 then
                continue;
            end if;
            for r0 in vals do
                total +:= 1;
                if shellOnly and TripleHeight(p0, z0, r0) lt H then
                    continue;
                end if;
                active +:= 1;
            end for;
        end for;
    end for;
    return total, active;
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

function ExpandedDegreesInVariable(facs, var)
    degs := [];
    for fe in facs do
        for j in [1..fe[2]] do
            Append(~degs, Degree(fe[1], var));
        end for;
    end for;
    return Sort(degs);
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

function QuadraticSubfactors(g)
    gm := MakeMonic(g);
    fac := Factorization(gm);
    atoms := [];
    for fe in fac do
        for i in [1..fe[2]] do
            Append(~atoms, MakeMonic(fe[1]));
        end for;
    end for;

    out := [];
    if #atoms eq 0 then
        return out;
    end if;

    for mask in [1..2^#atoms - 1] do
        h := Parent(gm)!1;
        for i in [1..#atoms] do
            if ((mask div 2^(i - 1)) mod 2) eq 1 then
                h *:= atoms[i];
            end if;
        end for;
        if Degree(h) eq 2 then
            h := MakeMonic(h);
            if &and[h ne old : old in out] then
                Append(~out, h);
            end if;
        end if;
    end for;
    return out;
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

function IsClosedBest(p0, z0, r0)
    closed := [
        <QQ!(-1)/QQ!3, QQ!(-1), QQ!4/QQ!3>,
        <QQ!(-1)/QQ!3, QQ!1, QQ!4/QQ!3>,
        <QQ!1/QQ!3, QQ!(-1), QQ!(-4)/QQ!3>,
        <QQ!1/QQ!3, QQ!1, QQ!(-4)/QQ!3>
    ];
    for T in closed do
        if p0 eq T[1] and z0 eq T[2] and r0 eq T[3] then
            return true;
        end if;
    end for;
    return false;
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

function SpecializeS(S, Mv, Nv)
    out := P!0;
    for j in [0..Degree(S)] do
        out +:= Evaluate(Coefficient(S, j), [Mv, Nv])*x^j;
    end for;
    return out;
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

procedure CertifyCandidate(row, p0, z0, r0, twoLabel, D12, J, f, S, pts,
        ~rationalPointRows, ~exactDivisibleRows, ~torsionCertRows)
    if #pts eq 0 then
        return;
    end if;

    rationalPointRows +:= 1;
    uv := Eltseq(D12);
    printf "RATIONAL_HALVING_CANDIDATES row=%o p=%o z=%o r=%o two_label=%o points=%o\n",
        row, p0, z0, r0, twoLabel, pts;

    for pt in pts do
        Mv := pt[1];
        Nv := pt[2];
        candidateLine := uv[2] + uv[1]*(Mv*x + Nv);
        Spt := SpecializeS(S, Mv, Nv);
        printf "  CANDIDATE_POINT M=%o N=%o line_y=%o S_special=%o S_factorization=%o\n",
            Mv, Nv, candidateLine, Spt, Factorization(Spt);
    end for;

    okI, DI, fI, denI := IntegralSquareModelDivisor(f, D12);
    if okI then
        divisible := IsDivisibleBy(DI, 2);
        if divisible then
            exactDivisibleRows +:= 1;
            _, half := IsDivisibleBy(DI, 2);
            printf "  EXACT_HALVING_CERT integral_den=%o half_order=%o twice_half_matches=%o\n",
                denI, Order(half), 2*half eq DI;
        else
            printf "  EXACT_HALVING_CERT integral_den=%o IsDivisibleBy(D,2)=false\n",
                denI;
        end if;
    else
        printf "  EXACT_HALVING_CERT skipped_integral_model_conversion_failed\n";
    end if;

    try
        TI, mp := TorsionSubgroup(Jacobian(HyperellipticCurve(fI)));
        invs := Invariants(TI);
        torsionCertRows +:= 1;
        printf "  EXACT_TORSION integral_square_model_invariants=%o order=%o contains_Z2xZ24=%o fI=%o\n",
            invs, TorsionOrderFromInvariants(invs), ContainsSubgroup(invs, [2, 24]), fI;
    catch e
        printf "  EXACT_TORSION_FAILED error=%o\n", e`Object;
    end try;
end procedure;

procedure AnalyzeTranslatedClass(row, p0, z0, r0, source, tlabel, twoIndex,
        twoCount, F, f, D12, J, ~lowRows, ~degree4Rows, ~rationalMRootRows,
        ~rationalMRootFactors, ~rationalPointRows, ~exactDivisibleRows,
        ~torsionCertRows, ~positiveDimRows, ~errorRows, ~degreeKeys,
        ~degreeCounts, ~minDegreeKeys, ~minDegreeCounts)
    uv := Eltseq(D12);
    if #uv lt 2 or Degree(uv[1]) ne 2 then
        errorRows +:= 1;
        printf "ROW_ERROR row=%o source=%o p=%o z=%o r=%o two_label=%o reason=non_degree2_mumford u_degree=%o\n",
            row, source, p0, z0, r0, tlabel, (#uv ge 1 select Degree(uv[1]) else -1);
        return;
    end if;

    ok, E1, E0, S := HalvingEquations(uv[1], uv[2], f);
    if not ok then
        errorRows +:= 1;
        printf "ROW_ERROR row=%o source=%o p=%o z=%o r=%o two_label=%o reason=halving_equation_setup_failed\n",
            row, source, p0, z0, r0, tlabel;
        return;
    end if;

    A := Parent(E1);
    M := A.1;
    N := A.2;
    s4 := Coefficient(S, 4);
    gcdDeg := TotalDegree(GCD(E1, E0));

    try
        RN := Resultant(E1, E0, N);
        if RN eq 0 then
            positiveDimRows +:= 1;
            lowRows +:= 1;
            BumpKey("positive_dimensional", ~degreeKeys, ~degreeCounts);
            BumpKey("positive_dimensional", ~minDegreeKeys, ~minDegreeCounts);
            printf "LOW_BRANCH row=%o source=%o p=%o z=%o r=%o two_label=%o reason=resultant_zero_or_positive_dim gcd_E_degree=%o\n",
                row, source, p0, z0, r0, tlabel, gcdDeg;
            return;
        end if;

        RNfac := Factorization(RN);
        affineFacs := [<fe[1], fe[2]> : fe in RNfac | GCD(fe[1], s4) eq A!1];
        boundaryFacs := [<fe[1], fe[2]> : fe in RNfac | GCD(fe[1], s4) ne A!1];
        degs := ExpandedDegreesInVariable(affineFacs, M);
        degKey := Sprint(degs);
        minAffine := (#degs eq 0) select -1 else Min(degs);
        minKey := Sprint(minAffine);
        affineDegree := SumInts(degs);
        hasLow := minAffine gt 0 and minAffine le LowDegree;
        hasDegree4 := minAffine gt 0 and minAffine le 4;
        rootsAll := [];
        rootLiftRecords := [];
        factorIndex := 0;

        BumpKey(degKey, ~degreeKeys, ~degreeCounts);
        BumpKey(minKey, ~minDegreeKeys, ~minDegreeCounts);

        for fe in affineFacs do
            factorIndex +:= 1;
            qM := MakeMonic(ToUnivariateInM(fe[1]));
            roots := Roots(qM);
            if #roots gt 0 then
                rationalMRootFactors +:= 1;
                for rt in roots do
                    Append(~rootsAll, rt[1]);
                    ptsM, status := RationalSolutionsAtM(E1, E0, S, rt[1]);
                    Append(~rootLiftRecords, <factorIndex, Degree(qM), rt[1], status, ptsM>);
                    CertifyCandidate(row, p0, z0, r0, tlabel, D12, J, f, S, ptsM,
                        ~rationalPointRows, ~exactDivisibleRows, ~torsionCertRows);
                end for;
            end if;
        end for;

        if hasLow then
            lowRows +:= 1;
        end if;
        if hasDegree4 then
            degree4Rows +:= 1;
        end if;
        if #rootsAll gt 0 then
            rationalMRootRows +:= 1;
        end if;

        if VerboseRows or hasLow or #rootsAll gt 0 or gcdDeg gt 0 then
            printf "ROW_SUMMARY row=%o source=%o p=%o z=%o r=%o two_label=%o two_index=%o two_count=%o F_factor_degrees=%o s4_factor_degrees=%o raw_degree_M=%o boundary_factor_degrees=%o saturated_degree=%o affine_factor_degrees=%o min_affine_factor_degree=%o rational_M_roots=%o gcd_E_degree=%o\n",
                row, source, p0, z0, r0, tlabel, twoIndex, twoCount,
                FactorDegreeMults(F), FactorDegreeMults(s4), Degree(RN, M),
                FactorDegreeMultsInVariable(boundaryFacs, M), affineDegree,
                FactorDegreeMultsInVariable(affineFacs, M), minAffine,
                rootsAll, gcdDeg;
        end if;

        if hasLow then
            printf "LOW_BRANCH row=%o source=%o p=%o z=%o r=%o two_label=%o affine_factor_degrees=%o min_affine_factor_degree=%o rational_M_roots=%o\n",
                row, source, p0, z0, r0, tlabel,
                FactorDegreeMultsInVariable(affineFacs, M), minAffine,
                rootsAll;
        end if;

        if #rootsAll gt 0 then
            printf "RATIONAL_M_ROOT_BRANCH row=%o source=%o p=%o z=%o r=%o two_label=%o affine_factor_degrees=%o lift_records=%o\n",
                row, source, p0, z0, r0, tlabel,
                FactorDegreeMultsInVariable(affineFacs, M), rootLiftRecords;
        end if;

        if DirectSearchHeight gt 0 and (hasLow or #rootsAll gt 0) then
            pts := DirectRationalSolutions(E1, E0, S, DirectSearchHeight);
            printf "DIRECT_SEARCH row=%o height=%o saturated_points=%o points=%o\n",
                row, DirectSearchHeight, #pts, pts;
            CertifyCandidate(row, p0, z0, r0, tlabel, D12, J, f, S, pts,
                ~rationalPointRows, ~exactDivisibleRows, ~torsionCertRows);
        end if;
    catch e
        errorRows +:= 1;
        printf "ROW_ERROR row=%o source=%o p=%o z=%o r=%o two_label=%o reason=resultant_or_factorization_failed error=%o\n",
            row, source, p0, z0, r0, tlabel, e`Object;
    end try;
end procedure;

procedure TryFiber(p0, z0, r0, source, ~rows, ~fibers, ~splitFibers,
        ~lowRows, ~degree4Rows, ~rationalMRootRows, ~rationalMRootFactors,
        ~rationalPointRows, ~exactDivisibleRows, ~torsionCertRows,
        ~positiveDimRows, ~errorRows, ~degreeKeys, ~degreeCounts,
        ~minDegreeKeys, ~minDegreeCounts, ~quarticKeys, ~quarticCounts,
        ~twoCountKeys, ~twoCountCounts)
    ok, Rpol, Qpol, ell, F, f := A12Data(p0, z0, r0);
    if not ok or Degree(f) ne 6 or Discriminant(f) eq 0 then
        return;
    end if;

    qsubs := QuadraticSubfactors(F);
    if #qsubs eq 0 then
        return;
    end if;

    splitFibers +:= 1;
    BumpKey(Sprint(FactorDegreeMults(F)), ~quarticKeys, ~quarticCounts);
    try
        J, O, TR, P12 := A12JacobianData(Rpol, Qpol, ell, f);
        if not (12*P12 eq O and &and[n*P12 ne O : n in [1..11]]) then
            return;
        end if;
        fibers +:= 1;

        twos := RationalTwoTorsionClasses(J, f);
        BumpKey(Sprint(#twos), ~twoCountKeys, ~twoCountCounts);
        printf "FIBER source=%o p=%o z=%o r=%o F_factor_degrees=%o quadratic_subfactors=%o two_count=%o\n",
            source, p0, z0, r0, FactorDegreeMults(F), #qsubs, #twos;

        for i in [1..#twos] do
            if MaxRows gt 0 and rows ge MaxRows then
                return;
            end if;
            D12 := P12 + twos[i];
            if not (12*D12 eq O and &and[n*D12 ne O : n in [1..11]]) then
                continue;
            end if;
            rows +:= 1;
            AnalyzeTranslatedClass(rows, p0, z0, r0, source, TwoLabel(twos[i], O, TR),
                i, #twos, F, f, D12, J, ~lowRows, ~degree4Rows,
                ~rationalMRootRows, ~rationalMRootFactors, ~rationalPointRows,
                ~exactDivisibleRows, ~torsionCertRows, ~positiveDimRows,
                ~errorRows, ~degreeKeys, ~degreeCounts, ~minDegreeKeys,
                ~minDegreeCounts);
        end for;
    catch e
        errorRows +:= 1;
        printf "FIBER_ERROR source=%o p=%o z=%o r=%o error=%o\n",
            source, p0, z0, r0, e`Object;
    end try;
end procedure;

print "A2_24_ALT_SCOUT_START";
printf "parameters Height=%o ShellOnly=%o LowDegree=%o SkipClosedBest=%o PStart=%o PStop=%o MaxRows=%o MaxFibers=%o MaxChecked=%o Progress=%o VerboseRows=%o DirectSearchHeight=%o\n",
    Height, ShellOnly, LowDegree, SkipClosedBest, PStart, PStop, MaxRows, MaxFibers,
    MaxChecked, Progress, VerboseRows, DirectSearchHeight;

vals := HeightRationals(Height);
pvals := [v : v in vals | v ne 0];
if PStop eq 0 or PStop gt #pvals then
    PStop := #pvals;
end if;
if PStart lt 1 then
    PStart := 1;
end if;
fullEnumerated, fullActive := CountActiveTriples(vals, pvals, 1, #pvals, Height, ShellOnly);
partitionEnumerated, partitionActive := CountActiveTriples(vals, pvals, PStart, PStop, Height, ShellOnly);
printf "height_values_count=%o nonzero_p_count=%o active_p_index_range=[%o,%o]\n",
    #vals, #pvals, PStart, PStop;
printf "partition_estimate full_enumerated_triples=%o full_active_triples=%o partition_enumerated_triples=%o partition_active_triples=%o active_fraction=%o\n",
    fullEnumerated, fullActive, partitionEnumerated, partitionActive,
    (fullActive eq 0 select RealField(12)!0 else RealField(12)!partitionActive/RealField(12)!fullActive);

rows := 0;
fibers := 0;
splitFibers := 0;
lowRows := 0;
degree4Rows := 0;
rationalMRootRows := 0;
rationalMRootFactors := 0;
rationalPointRows := 0;
exactDivisibleRows := 0;
torsionCertRows := 0;
positiveDimRows := 0;
errorRows := 0;
checked := 0;
closedSkipped := 0;
degreeKeys := [];
degreeCounts := [];
minDegreeKeys := [];
minDegreeCounts := [];
quarticKeys := [];
quarticCounts := [];
twoCountKeys := [];
twoCountCounts := [];
enumerated := 0;
shellSkipped := 0;
stopNow := false;

for pIndex in [PStart..PStop] do
    if stopNow then
        break;
    end if;
    p0 := pvals[pIndex];
    for z0 in vals do
        if stopNow then
            break;
        end if;
        if z0 eq 0 then
            continue;
        end if;
        for r0 in vals do
            enumerated +:= 1;
            if ShellOnly and TripleHeight(p0, z0, r0) lt Height then
                shellSkipped +:= 1;
                continue;
            end if;
            if MaxChecked gt 0 and checked ge MaxChecked then
                stopNow := true;
                break;
            end if;
            checked +:= 1;
            if Progress gt 0 and checked mod Progress eq 0 then
                printf "SCAN_PROGRESS active_checked=%o enumerated=%o shell_skipped=%o split_fibers=%o order12_split_fibers=%o rows=%o low_rows=%o degree4_rows=%o rational_M_root_rows=%o errors=%o cputime=%o\n",
                    checked, enumerated, shellSkipped, splitFibers, fibers, rows, lowRows, degree4Rows,
                    rationalMRootRows, errorRows, Cputime();
            end if;

            if MaxRows gt 0 and rows ge MaxRows then
                stopNow := true;
                break;
            end if;
            if MaxFibers gt 0 and fibers ge MaxFibers then
                stopNow := true;
                break;
            end if;

            if SkipClosedBest and IsClosedBest(p0, z0, r0) then
                closedSkipped +:= 1;
                printf "CLOSED_BEST_SKIPPED p=%o z=%o r=%o known_affine_factor_degrees_extra=[4,4,8] known_rational_M_roots=[]\n",
                    p0, z0, r0;
                continue;
            end if;

            TryFiber(p0, z0, r0, (ShellOnly select "height_shell" else "height_box"),
                ~rows, ~fibers, ~splitFibers,
                ~lowRows, ~degree4Rows, ~rationalMRootRows,
                ~rationalMRootFactors, ~rationalPointRows, ~exactDivisibleRows,
                ~torsionCertRows, ~positiveDimRows, ~errorRows, ~degreeKeys,
                ~degreeCounts, ~minDegreeKeys, ~minDegreeCounts, ~quarticKeys,
                ~quarticCounts, ~twoCountKeys, ~twoCountCounts);
        end for;
    end for;
end for;

print "AFFINE_DEGREE_DISTRIBUTION_START";
for i in [1..#degreeKeys] do
    printf "AFFINE_DEGREE_DISTRIBUTION degrees=%o count=%o\n",
        degreeKeys[i], degreeCounts[i];
end for;
print "AFFINE_DEGREE_DISTRIBUTION_END";

print "MIN_AFFINE_DEGREE_DISTRIBUTION_START";
for i in [1..#minDegreeKeys] do
    printf "MIN_AFFINE_DEGREE_DISTRIBUTION min_degree=%o count=%o\n",
        minDegreeKeys[i], minDegreeCounts[i];
end for;
print "MIN_AFFINE_DEGREE_DISTRIBUTION_END";

print "RESIDUAL_QUARTIC_FACTOR_DISTRIBUTION_START";
for i in [1..#quarticKeys] do
    printf "RESIDUAL_QUARTIC_FACTOR_DISTRIBUTION degrees=%o count=%o\n",
        quarticKeys[i], quarticCounts[i];
end for;
print "RESIDUAL_QUARTIC_FACTOR_DISTRIBUTION_END";

print "RATIONAL_TWO_TORSION_COUNT_DISTRIBUTION_START";
for i in [1..#twoCountKeys] do
    printf "RATIONAL_TWO_TORSION_COUNT_DISTRIBUTION two_count=%o count=%o\n",
        twoCountKeys[i], twoCountCounts[i];
end for;
print "RATIONAL_TWO_TORSION_COUNT_DISTRIBUTION_END";

printf "KNOWN_CLOSED_BEST_SKIPPED closed_fibers=%o known_degree4_translated_rows_if_full_four=%o known_quartic_branches_if_full_four=%o\n",
    closedSkipped, 2*closedSkipped, 4*closedSkipped;
printf "A2_24_ALT_SCOUT_DONE enumerated=%o shell_skipped=%o active_checked=%o estimated_partition_active=%o estimated_full_active=%o split_fibers=%o order12_split_fibers=%o translated_order12_rows=%o low_rows_le_%o=%o degree4_or_less_rows=%o rational_M_root_rows=%o rational_M_root_factors=%o rational_point_rows=%o exact_divisible_rows=%o torsion_cert_rows=%o positive_dim_rows=%o errors=%o closed_best_skipped=%o stopped_early=%o cputime=%o\n",
    enumerated, shellSkipped, checked, partitionActive, fullActive,
    splitFibers, fibers, rows, LowDegree, lowRows, degree4Rows,
    rationalMRootRows, rationalMRootFactors, rationalPointRows,
    exactDivisibleRows, torsionCertRows, positiveDimRows, errorRows,
    closedSkipped, stopNow, Cputime();
