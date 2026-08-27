//////////////////////////////////////////////////////////////////////
// Self-contained verifier for the Prym rank-1 dependency (review item:
// notes/m612_review_and_top3_plan_2026_07_13.md, Attack 1.2).
//   D    : y^2 = -3x^6+24x^3-75 (minimized Prym curve of E8/E4)
//   rank : RankBounds = 1..1 ; torsion trivial;
//   generator: explicit infinite-order Mumford class;
//   model transformations used by the bigonal transport.
// Usage: magma -b code/contact6_m612_prym_rank_verifier.m
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetMemoryLimit(4*10^9);
Q := Rationals(); Px<x> := PolynomialRing(Q);
fmin := -3*x^6 + 24*x^3 - 75;
Cmin := HyperellipticCurve(fmin);
Jmin := Jacobian(Cmin);
rlo, rhi := RankBounds(Jmin);
printf "RankBounds(J(D)) = %o..%o\n", rlo, rhi;
assert rlo eq 1 and rhi eq 1;
T := TorsionSubgroup(Jmin);
printf "torsion invariants = %o\n", Invariants(T);
assert #Invariants(T) eq 0;
pts := Points(Jmin : Bound := 2000);
gens := [P : P in pts | Order(P) eq 0];
assert #gens ge 1;
g0 := gens[1];
printf "generator (Mumford): (%o, %o)\n", g0[1], g0[2];
assert g0[1] eq x^2 + 2*x + 4 or g0[1] eq x^2 + 2*x + 4;  // support x=-1+-sqrt(-3)
// model transformation to the parametrized dual-tower curve
f2param := -1/192*x^6 + 1/32*x^5 - 5/64*x^4 + 7/16*x^3 - 69/64*x^2 + 33/32*x - 555/64;
Cpar := HyperellipticCurve(f2param);
okI, mpI := IsIsomorphic(Cmin, Cpar);
assert okI;
print "Cmin ~ Cpar isomorphism:", DefiningPolynomials(mpI);
// the un-minimized integral model used earlier, for the record
fint := -192*x^6 + 1152*x^5 - 2880*x^4 + 16128*x^3 - 39744*x^2 + 38016*x - 319680;
okI2 := IsIsomorphic(HyperellipticCurve(fint), Cmin);
assert okI2;
print "integral model ~ minimized: true";
print "PRYM_RANK_VERIFIED";
quit;
