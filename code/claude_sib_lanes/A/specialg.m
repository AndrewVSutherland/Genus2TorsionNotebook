// specialg.m — probe: what distinguishes the special fibers g=725/288, 459/23?
// For fixed g the three surface conditions are elliptic:
//   E2(g): y^2 = q(u)((g^2+1)+(6-8u)g)               (cubic in u)
//   E1(g): y^2 = (u-1)q(u)((u-1)(g^2+1)+2ug)         (quartic)
//   E4(g): y^2 = (2u-1)q(u)((2u-1)(g^2+1)+(6-4u)g)   (quartic)
// Compute rank bounds/torsion of E2 (and quartic covers if pointed) at the
// special g's vs warm controls (many NEAR2, 0 hits) vs random g.
P<u> := PolynomialRing(Rationals());
q := 4*u^2 - 6*u + 3;
probe := function(g)
    f2 := q*((g^2+6*g+1) - 8*g*u);
    E2 := EllipticCurve(HyperellipticCurve(f2));
    r2lo, r2hi := RankBounds(E2);
    t2 := Invariants(TorsionSubgroup(E2));
    return <r2lo, r2hi, t2>;
end function;
for g in [725/288, 459/23, 32/9, 116/9, -75/11, 66/5, 2/1, 3/1, 7/5, -32/9, -459/23] do
    res := probe(g);
    printf "g=%o: E2 rank in [%o,%o], tors %o\n", g, res[1], res[2], res[3];
end for;
quit;
