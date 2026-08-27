//////////////////////////////////////////////////////////////////////
// Exact Magma verification of the packaged contact-6 [6,6] example.
//
// Run with:
//     magma code/contact6_m36_66_package.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
P<x> := PolynomialRing(Q);

F := 1872*x^5 - 3000*x^4 + 6969*x^3 - 1691*x^2 + 4875*x;
C := HyperellipticCurve(F);
J := Jacobian(C);

q := x^2 - 9/4*x + 25/4;
D := J![x - 1, 95];
E := J![q, 25*x - 300];
T := J![x, 0];
W := E + T;

assert Order(D) eq 6;
assert Order(E) eq 3;
assert Order(T) eq 2;
assert Order(W) eq 6;

generated := {i*D + j*W : i in [0..5], j in [0..5]};
assert #generated eq 36;

G, phi := TorsionSubgroup(J);
assert Invariants(G) eq [6, 6];

print "minimal_model", F;
print "factorization", Factorization(F);
print "D", D, "order", Order(D);
print "E", E, "order", Order(E);
print "T", T, "order", Order(T);
print "W=E+T", W, "order", Order(W);
print "generated_subgroup_order", #generated;
print "torsion_invariants", Invariants(G);

quit;
