// probe: X1(11) raw model + the sigma-congruence test on quadratic points.
// Kubert raw coordinates: c = s(r-1), b = rc; P = (0,0) on
// E(b,c): y^2 + (1-c)xy - by = x^3 - bx^2 has order 11 iff F11(r,s) = 0,
// where F11 is the relevant factor of psi_11(0).
// For each rational r0 with the s-fiber quadratic over Q, we get a genuine
// quadratic point (E/K, P11).  Gluing E with E^sigma along 2 (resp. 3)
// requires E[2] ~ E^sigma[2] over K (resp. E[3]): test the 2-division cubic
// field match (and the 3-division quartic field match) over K.
SetColumns(0);
SetMemoryLimit(8*10^9);
Q := Rationals();
Qrs<r,s> := FunctionField(Q, 2);

cS := s*(r-1); bS := r*cS;
E := EllipticCurve([1-cS, -bS, -bS, 0, 0]);
psi11 := DivisionPolynomial(E, 11);
F := Evaluate(psi11, 0);   // vanishing <=> 11*P = O (P=(0,0) not 2-torsion)
// clear to polynomial and factor
P2<R,S> := PolynomialRing(Q, 2);
Fn := Numerator(F);
hh := hom< Parent(Fn) -> P2 | [R, S] >;
FP := hh(Fn);
fac := Factorization(FP);
printf "psi11(0) factors (deg in R, deg in S, multiplicity):\n";
for f in fac do
    printf "  degR=%o degS=%o mult=%o : %o\n", Degree(f[1], R), Degree(f[1], S), f[2],
        (Degree(f[1],R)+Degree(f[1],S)) le 6 select f[1] else Parent(f[1])!0;
end for;
quit;
