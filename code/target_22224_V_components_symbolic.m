// Symbolic expansion of the universal L=0 contact component.
// The normal coordinate is delta=2*v-U and M=L^2.
Q<M,U,delta,e1,e2,e3,e4> := PolynomialRing(Rationals(), 7);
v := (U + delta)/2;
PP := 4*M*e1 + 12*(U^2 + v^2) - (M + 3*U)^2;
F1 := (M + 3*U)*PP + 16*v^3 - 8*U^3 - 48*U*v^2 - 8*M*e2;
F2 := PP^2 + 64*(M + 3*U)*v^3
      - 192*(U^2*v^2 + v^4) - 64*M*e3;
F3 := PP*v^3 - 12*U*v^4 - 4*M*e4;
print "F1", F1;
print "F2", F2;
print "F3", F3;

// First weighted blow-up: M=L^2, delta=L*s.  The exceptional equations
// are the coefficients of the smallest powers of L.
S<L,s,U0,E1,E2,E3,E4> := PolynomialRing(Rationals(), 7);
phi := hom<Q -> S | L^2, U0, L*s, E1, E2, E3, E4>;
G1 := phi(F1);
G2 := phi(F2);
G3 := phi(F3);
function LowestLCoefficient(f)
    mons := Monomials(f);
    n := Minimum([ Degree(m, 1) : m in mons ]);
    q := Evaluate(f div L^n, [0,s,U0,E1,E2,E3,E4]);
    return n, q;
end function;
n1,g1 := LowestLCoefficient(G1);
n2,g2 := LowestLCoefficient(G2);
n3,g3 := LowestLCoefficient(G3);
print "WEIGHTED_BLOWUP_F1", n1, Factorization(g1);
print "WEIGHTED_BLOWUP_F2", n2, Factorization(g2);
print "WEIGHTED_BLOWUP_F3", n3, Factorization(g3);
