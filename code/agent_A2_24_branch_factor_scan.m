/*
Bounded scanner for saturated A(2,12) -> A(2,24) branch fibers.

This is a small exact invariant scan, not a broad height search.  It enumerates
small rational A(12) parameters, keeps fibers where the residual quartic F has
a rational quadratic pairing, translates the visible order-12 class by all
rational 2-torsion classes found from the factorization, and records the
saturated square-quartic halving fiber after removing s4=0.

Run examples:
    magma code/agent_A2_24_branch_factor_scan.m
    magma -b Height:=5 MaxRows:=12 code/agent_A2_24_branch_factor_scan.m
*/

SetColumns(0);

QQ := Rationals();
P<x> := PolynomialRing(QQ);

if assigned Height and Type(Height) eq MonStgElt then
    Height := StringToInteger(Height);
end if;
if not assigned Height then
    Height := 4;
end if;

if assigned MaxRows and Type(MaxRows) eq MonStgElt then
    MaxRows := StringToInteger(MaxRows);
end if;
if not assigned MaxRows then
    MaxRows := 12;
end if;

if assigned Progress and Type(Progress) eq MonStgElt then
    Progress := StringToInteger(Progress);
end if;
if not assigned Progress then
    Progress := 5000;
end if;

if assigned IncludeKnown and Type(IncludeKnown) eq MonStgElt then
    IncludeKnown := IncludeKnown in {"true", "True", "1", "yes", "Yes"};
end if;
if not assigned IncludeKnown then
    IncludeKnown := true;
end if;

function SumInts(xs)
    if #xs eq 0 then
        return 0;
    end if;
    return &+xs;
end function;

function ContainsPoint(seq, T)
    for S in seq do
        if S eq T then
            return true;
        end if;
    end for;
    return false;
end function;

function MakeMonic(g)
    return g/LeadingCoefficient(g);
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

procedure AnalyzeTranslatedClass(row, p0, z0, r0, tlabel, twoIndex, twoCount,
        F, f, D12, ~lowRows, ~errorRows)
    uv := Eltseq(D12);
    if #uv lt 2 or Degree(uv[1]) ne 2 then
        printf "BRANCH_SKIP row=%o p=%o z=%o r=%o two_index=%o reason=non_degree2_mumford u_degree=%o\n",
            row, p0, z0, r0, twoIndex, (#uv ge 1 select Degree(uv[1]) else -1);
        return;
    end if;

    ok, E1, E0, S := HalvingEquations(uv[1], uv[2], f);
    if not ok then
        errorRows +:= 1;
        printf "BRANCH_ERROR row=%o p=%o z=%o r=%o two_index=%o reason=halving_equation_setup_failed\n",
            row, p0, z0, r0, twoIndex;
        return;
    end if;

    A := Parent(E1);
    M := A.1;
    N := A.2;
    s4 := Coefficient(S, 4);
    gcdE := GCD(E1, E0);
    gcdDeg := TotalDegree(gcdE);

    try
        RN := Resultant(E1, E0, N);
        if RN eq 0 then
            lowRows +:= 1;
            printf "BRANCH_ROW row=%o p=%o z=%o r=%o two_label=%o two_index=%o two_count=%o F_factor_degrees=%o s4=%o raw_degree_N=0 saturated_degree=positive_dim affine_factor_degrees=[] gcd_E_degree=%o quotient_flag=true\n",
                row, p0, z0, r0, tlabel, twoIndex, twoCount,
                FactorDegreeMults(F), s4, gcdDeg;
            printf "LOW_DEGREE_BRANCH row=%o reason=resultant_zero_or_positive_dim\n", row;
            return;
        end if;

        RNfac := Factorization(RN);
        affineFacs := [<Degree(fe[1], M), fe[2]> : fe in RNfac | GCD(fe[1], s4) eq 1];
        boundaryFacs := [<Degree(fe[1], M), fe[2]> : fe in RNfac | GCD(fe[1], s4) ne 1];
        affineDegree := SumInts([d[1]*d[2] : d in affineFacs]);
        minAffine := (#affineFacs eq 0) select -1 else Min([d[1] : d in affineFacs]);

        RMfacDegrees := [];
        try
            RM := Resultant(E1, E0, M);
            if RM ne 0 then
                RMfacDegrees := [<Degree(fe[1], N), fe[2]> : fe in Factorization(RM)];
            end if;
        catch e
            RMfacDegrees := [];
        end try;

        low := (minAffine gt 0 and minAffine le 8) or gcdDeg gt 0;
        if low then
            lowRows +:= 1;
        end if;

        printf "BRANCH_ROW row=%o p=%o z=%o r=%o two_label=%o two_index=%o two_count=%o F_factor_degrees=%o s4_factor_degrees=%o raw_degree_N=%o boundary_factor_degrees=%o saturated_degree=%o affine_factor_degrees=%o min_affine_factor_degree=%o resultant_M_factor_degrees=%o gcd_E_degree=%o low_degree_or_quotient=%o\n",
            row, p0, z0, r0, tlabel, twoIndex, twoCount,
            FactorDegreeMults(F), FactorDegreeMults(s4), Degree(RN, M),
            boundaryFacs, affineDegree, affineFacs, minAffine, RMfacDegrees,
            gcdDeg, low;

        if low then
            printf "LOW_DEGREE_BRANCH row=%o p=%o z=%o r=%o two_label=%o affine_factor_degrees=%o gcd_E_degree=%o\n",
                row, p0, z0, r0, tlabel, affineFacs, gcdDeg;
        end if;
    catch e
        errorRows +:= 1;
        printf "BRANCH_ERROR row=%o p=%o z=%o r=%o two_index=%o reason=resultant_or_factorization_failed error=%o\n",
            row, p0, z0, r0, twoIndex, e`Object;
    end try;
end procedure;

procedure TryFiber(p0, z0, r0, source, ~rows, ~fibers, ~splitFibers,
        ~lowRows, ~errorRows)
    ok, Rpol, Qpol, ell, F, f := A12Data(p0, z0, r0);
    if not ok or Degree(f) ne 6 or Discriminant(f) eq 0 then
        return;
    end if;

    qsubs := QuadraticSubfactors(F);
    if #qsubs eq 0 then
        return;
    end if;

    splitFibers +:= 1;
    try
        J, O, TR, P12 := A12JacobianData(Rpol, Qpol, ell, f);
        if not (12*P12 eq O and &and[n*P12 ne O : n in [1..11]]) then
            return;
        end if;
        fibers +:= 1;

        twos := RationalTwoTorsionClasses(J, f);
        printf "FIBER source=%o p=%o z=%o r=%o F_factor_degrees=%o quadratic_subfactors=%o two_count=%o\n",
            source, p0, z0, r0, FactorDegreeMults(F), #qsubs, #twos;

        for i in [1..#twos] do
            if rows ge MaxRows then
                return;
            end if;
            D12 := P12 + twos[i];
            if not (12*D12 eq O and &and[n*D12 ne O : n in [1..11]]) then
                continue;
            end if;
            rows +:= 1;
            AnalyzeTranslatedClass(rows, p0, z0, r0, TwoLabel(twos[i], O, TR),
                i, #twos, F, f, D12, ~lowRows, ~errorRows);
        end for;
    catch e
        errorRows +:= 1;
        printf "FIBER_ERROR source=%o p=%o z=%o r=%o error=%o\n",
            source, p0, z0, r0, e`Object;
    end try;
end procedure;

print "A2_24_BRANCH_FACTOR_SCAN_START";
printf "parameters Height=%o MaxRows=%o IncludeKnown=%o Progress=%o\n",
    Height, MaxRows, IncludeKnown, Progress;

rows := 0;
fibers := 0;
splitFibers := 0;
lowRows := 0;
errorRows := 0;
checked := 0;

if IncludeKnown then
    TryFiber(QQ!(-5)/QQ!3, QQ!1, QQ!2/QQ!3, "known_p-5/3_z1_r2/3",
        ~rows, ~fibers, ~splitFibers, ~lowRows, ~errorRows);
end if;

vals := HeightRationals(Height);
for p0 in vals do
    if rows ge MaxRows then
        break;
    end if;
    if p0 eq 0 then
        continue;
    end if;
    for z0 in vals do
        if rows ge MaxRows then
            break;
        end if;
        if z0 eq 0 then
            continue;
        end if;
        for r0 in vals do
            if rows ge MaxRows then
                break;
            end if;
            checked +:= 1;
            if Progress gt 0 and checked mod Progress eq 0 then
                printf "SCAN_PROGRESS checked=%o split_fibers=%o order12_split_fibers=%o rows=%o low_rows=%o errors=%o cputime=%o\n",
                    checked, splitFibers, fibers, rows, lowRows, errorRows, Cputime();
            end if;
            if IncludeKnown and p0 eq QQ!(-5)/QQ!3 and z0 eq QQ!1 and r0 eq QQ!2/QQ!3 then
                continue;
            end if;
            TryFiber(p0, z0, r0, "height_scan", ~rows, ~fibers,
                ~splitFibers, ~lowRows, ~errorRows);
        end for;
    end for;
end for;

printf "A2_24_BRANCH_FACTOR_SCAN_DONE checked=%o split_fibers=%o order12_split_fibers=%o rows=%o low_rows=%o errors=%o cputime=%o\n",
    checked, splitFibers, fibers, rows, lowRows, errorRows, Cputime();
