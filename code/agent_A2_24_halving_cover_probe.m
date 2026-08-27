/*
Structural probe for the A(12) -> A(2,24) halving cover.

This is not another height search.  It analyzes the one A12 local survivor
already printed in the height-10 logs, and separates the boundary s4=0 factors
from the true affine square-quartic halving cover for the four printed
translated order-12 classes.

Run:
    magma code/agent_A2_24_halving_cover_probe.m

Optional:
    magma -b PrintLargeFactors:=true code/agent_A2_24_halving_cover_probe.m
*/

SetColumns(0);

Q := Rationals();
P<x> := PolynomialRing(Q);

if assigned PrintLargeFactors and Type(PrintLargeFactors) eq MonStgElt then
    PrintLargeFactors := PrintLargeFactors in {"true", "True", "1", "yes", "Yes"};
end if;
if not assigned PrintLargeFactors then
    PrintLargeFactors := false;
end if;

f := 81/15625*x^6 + 369/5000*x^5 + 2499/10000*x^4
    - 1/18*x^3 - 151/648*x^2 - 625/648*x + 6875/11664;

classes := [
    <"A", x^2 - 25/6*x + 625/162, 32/9*x - 400/81>,
    <"B", x^2 + 1150/693*x + 3125/6237, -2949/5929*x + 13775/17787>,
    <"C", x^2 + 50/9*x - 625/27, 4/9*x + 100/27>,
    <"D", x^2 + 1100/189*x - 625/1701, 1601/7938*x - 52825/71442>
];

function SumInts(xs)
    if #xs eq 0 then
        return 0;
    end if;
    return &+xs;
end function;

function SquareIntegralPolynomial(g)
    den := 1;
    for c in Coefficients(g) do
        den := LCM(den, Denominator(c));
    end for;
    return Parent(g)!(den^2*g);
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

function HalvingEquations(u, v, f)
    A<M,N> := PolynomialRing(Q, 2);
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

function SpecializeS(S, Mv, Nv)
    out := P!0;
    for i in [0..Degree(S)] do
        out +:= Evaluate(Coefficient(S, i), [Mv, Nv])*x^i;
    end for;
    return out;
end function;

function SquareQuarticDataQ(S)
    RY := Parent(S);
    Y := RY.1;
    if Degree(S) ne 4 then
        return false, Q!0, RY!0;
    end if;
    s4 := Coefficient(S, 4);
    if s4 eq 0 then
        return false, Q!0, RY!0;
    end if;
    s3 := Coefficient(S, 3);
    s2 := Coefficient(S, 2);
    G := Y^2 + (s3/(2*s4))*Y + (4*s4*s2 - s3^2)/(8*s4^2);
    return S eq s4*G^2, s4, G;
end function;

procedure AnalyzeClass(label, u, v, f, J, O)
    ok, E1, E0, S := HalvingEquations(u, v, f);
    if not ok then
        printf "\nCLASS %o setup_failed\n", label;
        return;
    end if;

    A := Parent(E1);
    M := A.1;
    N := A.2;
    s4 := Coefficient(S, 4);

    printf "\nCLASS %o\n", label;
    printf "  u=%o\n  v=%o\n", u, v;

    D := J![u, v];
    printf "  exact_order12=%o\n", 12*D eq O and &and[n*D ne O : n in [1..11]];
    printf "  E_degrees E1=(total %o, M %o, N %o), E0=(total %o, M %o, N %o)\n",
        TotalDegree(E1), Degree(E1, M), Degree(E1, N),
        TotalDegree(E0), Degree(E0, M), Degree(E0, N);
    printf "  s4=%o factorization=%o\n", s4, Factorization(s4);
    printf "  gcd(E1,E0)=%o\n", GCD(E1, E0);

    I := ideal<A | E1, E0>;
    try
        printf "  ideal_dimension=%o\n", Dimension(I);
    catch e
        printf "  ideal_dimension_failed=%o\n", e`Object;
    end try;

    Ib := ideal<A | E1, E0, s4>;
    try
        bpts := Variety(Ib);
        printf "  boundary_s4_rational_points=%o\n", #bpts;
        for pt in bpts do
            Mv := Q!pt[1];
            Nv := Q!pt[2];
            Sq := SpecializeS(S, Mv, Nv);
            sq_ok, lambda, G := SquareQuarticDataQ(Sq);
            printf "    boundary_point M=%o N=%o deg(S)=%o square_quartic=%o\n",
                Mv, Nv, Degree(Sq), sq_ok;
        end for;
    catch e
        printf "  boundary_s4_variety_failed=%o\n", e`Object;
    end try;

    try
        RN := Resultant(E1, E0, N);
        RNfac := Factorization(RN);
        affDegM := SumInts([Degree(fe[1], M)*fe[2] : fe in RNfac | GCD(fe[1], s4) eq 1]);
        boundaryFacs := [<fe[1], fe[2]> : fe in RNfac | GCD(fe[1], s4) ne 1];
        affineSummary := [<Degree(fe[1], M), fe[2]> : fe in RNfac | GCD(fe[1], s4) eq 1];
        printf "  resultant_eliminate_N degree_M=%o\n", Degree(RN, M);
        printf "  resultant_N_boundary_factors=%o\n", boundaryFacs;
        printf "  resultant_N_affine_degree_after_s4=%o affine_factor_degrees=%o\n",
            affDegM, affineSummary;
        if PrintLargeFactors then
            printf "  resultant_N_factorization=%o\n", RNfac;
        end if;
    catch e
        printf "  resultant_eliminate_N_failed=%o\n", e`Object;
    end try;

    try
        RM := Resultant(E1, E0, M);
        RMfac := Factorization(RM);
        printf "  resultant_eliminate_M degree_N=%o factor_degrees=%o\n",
            Degree(RM, N), [<Degree(fe[1], N), fe[2]> : fe in RMfac];
        if PrintLargeFactors then
            printf "  resultant_M_factorization=%o\n", RMfac;
        end if;
    catch e
        printf "  resultant_eliminate_M_failed=%o\n", e`Object;
    end try;
end procedure;

print "A2_24_HALVING_COVER_PROBE_START";
printf "curve f=%o\n", f;
printf "factor_degrees=%o\n", [Degree(fac[1]) : fac in Factorization(f)];
okLead, sqrtLead := IsSquare(LeadingCoefficient(f));
printf "leading_coefficient=%o is_square=%o sqrt=%o\n",
    LeadingCoefficient(f), okLead, okLead select sqrtLead else 0;

C := HyperellipticCurve(f);
J := Jacobian(C);
O := J!0;

try
    A, mp := TorsionSubgroup(Jacobian(HyperellipticCurve(SquareIntegralPolynomial(f))));
    inv := Invariants(A);
    printf "torsion_invariants_integral_square_model=%o contains_Z2_Z24=%o\n",
        inv, ContainsSubgroup(inv, [2, 24]);
catch e
    printf "torsion_subgroup_failed=%o\n", e`Object;
end try;

for rec in classes do
    AnalyzeClass(rec[1], rec[2], rec[3], f, J, O);
end for;

print "\nA2_24_HALVING_COVER_PROBE_DONE";
