// probe: sigma-congruence feasibility on quadratic points of X1(11).
// X1(11) raw model (derived): r^2 - r(s^3-3s^2+4s) + s = 0, with
// c = s(r-1), b = rc, P = (0,0) of order 11 on E(b,c).
// For rational s0 the r-fiber is quadratic: E lives over K = Q(sqrt(D(s0))),
// D(s) = (s^3-3s^2+4s)^2 - 4s.  Tests per point:
//   (2) E[2] ~ E^sigma[2] over K   (<=> 2-division cubic fields match)
//   (3) x-quartic of E[3] matches sigma-conjugate (necessary for 3-glue)
//   (Q) E isogenous to E^sigma of degree 2,3,5,7 (Q-curve; mu-graph route)
SetColumns(0);
SetMemoryLimit(8*10^9);
Q := Rationals();
Px<x> := PolynomialRing(Q);

function HeightRats(H)
    S := [];
    for a in [1..H], b in [1..H] do
        if GCD(a,b) eq 1 then Append(~S, Q!a/b); Append(~S, Q!-a/b); end if;
    end for;
    return S;
end function;

MP := AssociativeArray();
for d in [2,3,5,7] do
    MP[d] := ClassicalModularPolynomial(d);
end for;

n2match := 0; n3match := 0; nqc := 0; ntot := 0;
for s0 in HeightRats(30) do
    D := (s0^3 - 3*s0^2 + 4*s0)^2 - 4*s0;
    if D eq 0 or IsSquare(D) then continue; end if;
    // squarefree kernel of D for the field
    dn := Numerator(D)*Denominator(D);
    K := QuadraticField(Squarefree(dn));
    RK := PolynomialRing(K);
    rr := Roots(RK![s0, -(s0^3-3*s0^2+4*s0), 1]);
    if #rr eq 0 then continue; end if;
    r0 := rr[1][1];
    c0 := s0*(r0-1); b0 := r0*c0;
    if b0 eq 0 then continue; end if;
    E := 0; ok := true;
    try E := EllipticCurve([1-c0, -b0, -b0, 0, 0]); catch e; ok := false; end try;
    if not ok then continue; end if;
    okP := true;
    try okP := Order(E![0,0]) eq 11; catch e; okP := false; end try;
    if not okP then continue; end if;
    ntot +:= 1;
    sig := hom< K -> K | -K.1 >;
    // 2-division cubic (monic integral-ish): x^3 + b2 x^2 + 8 b4 x + 16 b6
    a1 := 1-c0;
    b2 := a1^2 - 4*b0; b4 := -a1*b0; b6 := b0^2;
    cub := RK![16*b6, 8*b4, b2, 1];
    cubs := RK![sig(Coefficient(cub,i)) : i in [0..3]];
    // module iso test: sigma-cubic has a root in K[x]/(cub)?
    m2 := false;
    if IsIrreducible(cub) then
        L := ext< K | cub >;
        m2 := #Roots(PolynomialRing(L)!cubs) gt 0;
    else
        // rational 2-torsion over K would force Z/22: flag loudly
        printf "UNEXPECTED reducible 2-cubic at s0=%o\n", s0;
    end if;
    // 3-division x-quartic match (necessary condition for any 3-glue)
    psi3 := DivisionPolynomial(E, 3);
    psi3s := Parent(psi3)![sig(Coefficient(psi3,i)) : i in [0..Degree(psi3)]];
    m3 := false;
    f3 := Factorization(psi3);
    // compare multisets of factor degrees first; then root-in-field test on
    // the largest factor's field
    d3 := Sort([Degree(f[1]) : f in f3]);
    f3s := Factorization(psi3s);
    d3s := Sort([Degree(f[1]) : f in f3s]);
    if d3 eq d3s then
        big := f3[#f3][1];
        if Degree(big) eq 1 then
            m3 := true;   // all x-coords in K on both sides
        else
            L3 := ext< K | big >;
            m3 := exists{ f : f in f3s | Degree(f[1]) eq Degree(big) and #Roots(PolynomialRing(L3)!f[1]) gt 0 };
        end if;
    end if;
    // Q-curve test: Phi_d(j, j^sigma) = 0
    jv := jInvariant(E); js := sig(jv);
    isqc := jv eq js;   // degree 1 / defined over Q up to twist
    for d in [2,3,5,7] do
        if Evaluate(MP[d], [jv, js]) eq 0 then isqc := true; end if;
    end for;
    if m2 then n2match +:= 1; printf "MATCH2 s0=%o D=%o\n", s0, Squarefree(dn); end if;
    if m3 then n3match +:= 1; printf "MATCH3 s0=%o D=%o\n", s0, Squarefree(dn); end if;
    if isqc then nqc +:= 1; printf "QCURVE s0=%o j=%o\n", s0, jv; end if;
end for;
printf "X111_DONE points=%o match2=%o match3=%o qcurves=%o\n", ntot, n2match, n3match, nqc;
quit;
