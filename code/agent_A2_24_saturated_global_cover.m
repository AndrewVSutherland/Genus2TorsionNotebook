/*
Saturated A(2,24) halving cover over one A(2,12) slice.

Chosen A(2,12) branch/slice:
    z = 1, p = -5/3
and the fully split quadratic-factor cover of the residual quartic F.
This is the slice through the small split point recorded in
a2_12_resolvent.tex:
    (p,z,r) = (-5/3, 1, 2/3).

The translated order-12 class is
    D = P12 + T_extra,
where T_extra is the 2-torsion class represented by the recovered quadratic
factor U of F.  For D = [u,v], the saturated halving equations are built from
    ell = v + u*(M*x + N),
    S = (ell^2 - f)/u = s4*x^4 + ... + s0,
    E1 = 8*s4^2*s1 - s3*(4*s4*s2 - s3^2),
    E0 = 64*s4^3*s0 - (4*s4*s2 - s3^2)^2,
and the boundary s4=0 is removed before recording invariants.

Run:
    magma code/agent_A2_24_saturated_global_cover.m
*/

SetColumns(0);

QQ := Rationals();

if assigned PrimeList and Type(PrimeList) eq MonStgElt then
    PrimeList := [StringToInteger(s) : s in Split(PrimeList, ",")];
end if;
if not assigned PrimeList then
    PrimeList := [7, 11, 13, 17, 19, 23, 29, 37, 41];
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

function A12Data(K, p_in, z_in, r_in)
    P<X> := PolynomialRing(K);
    p := K!p_in;
    z := K!z_in;
    r := K!r_in;

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

    T := p*X + r;
    Rpol := (T^2 + X - 1)/lambda;
    ell := s*X + t;
    Qpol := 2*T + mu*Rpol;
    F := Rpol*X^2 + 4*(Rpol + X - 1)*(Rpol - 1);
    f := Rpol*F;

    return true, Rpol, Qpol, ell, F, f;
end function;

function QuarticResolvent(F)
    P := Parent(F);
    K := BaseRing(P);
    S<eta> := PolynomialRing(K);

    Fm := F/LeadingCoefficient(F);
    B := Coefficient(Fm, 3);
    C := Coefficient(Fm, 2);
    D := Coefficient(Fm, 1);
    E := Coefficient(Fm, 0);

    return eta^3 - C*eta^2 + (B*D - 4*E)*eta
        + (4*C*E - B^2*E - D^2);
end function;

function QuadraticFactorFromPairing(F, eta0, delta0)
    P := Parent(F);
    X := P.1;
    K := BaseRing(P);
    Fm := F/LeadingCoefficient(F);

    B := Coefficient(Fm, 3);
    D := Coefficient(Fm, 1);

    if delta0 eq 0 then
        return false, P!0, P!0;
    end if;

    u := (eta0 - delta0)/2;
    v := (eta0 + delta0)/2;
    a := (D - B*u)/delta0;
    U := X^2 + a*X + u;
    V := X^2 + (B - a)*X + v;

    return Fm eq U*V, U, V;
end function;

function A12JacobianData(Rpol, Qpol, ell, f)
    P := Parent(f);
    X := P.1;
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    O := J!0;

    TR := J![MakeMonic(Rpol), P!0];
    u4 := MakeMonic(Qpol);
    v4 := (Rpol*ell) mod u4;
    P4 := J![u4, v4];

    u6 := MakeMonic(Rpol + X - 1);
    v6 := (X*Rpol) mod u6;
    P6 := J![u6, v6];

    return J, O, TR, P4 + P6;
end function;

function HalvingEquations(u, v, f)
    K := BaseRing(Parent(f));
    A<M,N> := PolynomialRing(K, 2);
    RX<X> := PolynomialRing(A);
    phi := hom<Parent(f) -> RX | X>;

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

function CountAffineFiber(E1, E0, S)
    K := BaseRing(Parent(E1));
    s4 := Coefficient(S, 4);
    fiber := 0;
    boundary := 0;

    for Mv in K do
        for Nv in K do
            vals := [Mv, Nv];
            if Evaluate(E1, vals) eq 0 and Evaluate(E0, vals) eq 0 then
                if Evaluate(s4, vals) eq 0 then
                    boundary +:= 1;
                else
                    fiber +:= 1;
                end if;
            end if;
        end for;
    end for;

    return fiber, boundary;
end function;

procedure ExactFiberAtSplitPoint()
    ok, R0, Q0, ell0, F0, f0 := A12Data(QQ, QQ!(-5/3), QQ!1, QQ!(2/3));
    assert ok;
    assert Discriminant(f0) ne 0;

    P := Parent(f0);
    Ffac := Factorization(F0);
    assert #Ffac eq 2 and &and[Degree(fe[1]) eq 2 : fe in Ffac];

    J, O, TR, P12 := A12JacobianData(R0, Q0, ell0, f0);
    Textra := J![MakeMonic(Ffac[1][1]), P!0];
    D12 := P12 + Textra;
    uv := Eltseq(D12);

    okHalve, E1, E0, S := HalvingEquations(uv[1], uv[2], f0);
    assert okHalve;

    A := Parent(E1);
    M := A.1;
    N := A.2;
    s4 := Coefficient(S, 4);
    RN := Resultant(E1, E0, N);
    RNfac := Factorization(RN);
    affineFacs := [<Degree(fe[1], M), fe[2]> : fe in RNfac | GCD(fe[1], s4) eq 1];
    boundaryFacs := [<Degree(fe[1], M), fe[2]> : fe in RNfac | GCD(fe[1], s4) ne 1];
    affineDegree := SumInts([d[1]*d[2] : d in affineFacs]);

    printf "EXACT_SPLIT_FIBER p=-5/3 z=1 r=2/3\n";
    printf "  residual_factor_U=%o\n", MakeMonic(Ffac[1][1]);
    printf "  D=P12+T_extra u=%o\n", uv[1];
    printf "  D=P12+T_extra v=%o\n", uv[2];
    printf "  order_checks 12D=0:%o exact_order12:%o independent_extra2:%o\n",
        12*D12 eq O, (12*D12 eq O and &and[n*D12 ne O : n in [1..11]]),
        (2*Textra eq O and Textra ne O and Textra ne TR);
    printf "  E_degrees E1=(total %o, M %o, N %o), E0=(total %o, M %o, N %o)\n",
        TotalDegree(E1), Degree(E1, M), Degree(E1, N),
        TotalDegree(E0), Degree(E0, M), Degree(E0, N);
    printf "  s4=%o factorization=%o\n", s4, Factorization(s4);
    printf "  resultant_N_degree=%o boundary_factor_degrees=%o\n",
        Degree(RN, M), boundaryFacs;
    printf "  saturated_affine_degree_after_s4=%o affine_factor_degrees=%o\n",
        affineDegree, affineFacs;
    try
        Ibd := ideal<A | E1, E0, s4>;
        printf "  boundary_s4_rational_points=%o\n", #Variety(Ibd);
    catch e
        printf "  boundary_s4_rational_points_failed=%o\n", e`Object;
    end try;
end procedure;

procedure CountFiniteFieldSlice(q)
    Fq := GF(q);
    p0 := Fq!(-5)/Fq!3;
    z0 := Fq!1;

    basePoints := 0;
    smoothBasePoints := 0;
    splitPoints := 0;
    translatedOrder12 := 0;
    coverPoints := 0;
    boundarySolutions := 0;
    errors := 0;

    for r0 in Fq do
        ok, R0, Q0, ell0, F0, f0 := A12Data(Fq, p0, z0, r0);
        if not ok then
            continue;
        end if;

        phi := QuarticResolvent(F0);
        Fm := F0/LeadingCoefficient(F0);
        Econst := Coefficient(Fm, 0);

        for eta0 in Fq do
            if Evaluate(phi, eta0) ne 0 then
                continue;
            end if;
            basePoints +:= 1;

            for delta0 in Fq do
                if delta0^2 ne eta0^2 - 4*Econst then
                    continue;
                end if;
                if delta0 eq 0 then
                    continue;
                end if;

                okfac, U, V := QuadraticFactorFromPairing(F0, eta0, delta0);
                if not okfac then
                    continue;
                end if;
                splitPoints +:= 1;

                if Discriminant(f0) eq 0 then
                    continue;
                end if;
                smoothBasePoints +:= 1;

                try
                    P := Parent(f0);
                    J, O, TR, P12 := A12JacobianData(R0, Q0, ell0, f0);
                    Textra := J![U, P!0];
                    D12 := P12 + Textra;
                    if not (12*D12 eq O and &and[n*D12 ne O : n in [1..11]]) then
                        continue;
                    end if;
                    translatedOrder12 +:= 1;

                    uv := Eltseq(D12);
                    if Degree(uv[1]) ne 2 then
                        continue;
                    end if;

                    okHalve, E1, E0, S := HalvingEquations(uv[1], uv[2], f0);
                    if not okHalve then
                        continue;
                    end if;
                    fiber, boundary := CountAffineFiber(E1, E0, S);
                    coverPoints +:= fiber;
                    boundarySolutions +:= boundary;
                catch e
                    errors +:= 1;
                end try;
            end for;
        end for;
    end for;

    printf "FINITE_FIELD_SLICE q=%o base_resolvent_Fq_points=%o split_oriented_delta_ne0=%o smooth_split=%o translated_order12=%o saturated_cover_Fq_points=%o removed_s4_boundary_points=%o errors=%o\n",
        q, basePoints, splitPoints, smoothBasePoints, translatedOrder12,
        coverPoints, boundarySolutions, errors;
end procedure;

print "A2_24_SATURATED_GLOBAL_COVER_START";
print "chosen_slice p=-5/3 z=1; branch=fully split quadratic-factor A(2,12) cover";
print "translated_class D=P12+T_extra, where T_extra=[U,0] for the recovered quadratic factor U";
ExactFiberAtSplitPoint();
for q in PrimeList do
    if q in {2, 3, 5} then
        continue;
    end if;
    CountFiniteFieldSlice(q);
end for;
print "A2_24_SATURATED_GLOBAL_COVER_DONE";
