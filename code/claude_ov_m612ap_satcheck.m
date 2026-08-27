//////////////////////////////////////////////////////////////////////
// claude_ov_m612ap_satcheck.m    lane 9 ([6,12])   2026-07-25
//
// Cross-check on the saturation hypothesis of the Abel-Prym MW sieve.
//
// Pr ~ J(D),  D : y^2 = -3x^6+24x^3-75,  J(D)(Q) = <g0> = Z (rank 1,
// trivial torsion, g0 = (x^2+2x+4, 5x+5) in Mumford coordinates:
// code/contact6_m612_prym_rank_verifier.m, re-verified this session).
//
// The sieve assumes W generates Pr(Q) modulo torsion.  Under the isogeny,
// the reduction of a GENERATOR of the rank-1 group has order ord(g0 mod q)
// in J(D)(F_q) = Pr(F_q) (the L-polynomials agree, claude_ov_m612ap_prymid.m).
// If W = n * (generator) then ord(W mod q) = ord(g0 mod q) / gcd(n, ...),
// so ord(g0 mod q) / ord(W mod q) is (up to the degree of the isogeny and to
// the index of the isogeny image) the part of n visible at q.
//
// This prints ord(g0 mod q) beside the ord(W mod q) measured by the sieve.
// Equality at many q is evidence that W is already a generator (index 1);
// a constant ratio would be evidence that it is not.
//
// Usage: code/claude_magma_slot.sh -b code/claude_ov_m612ap_satcheck.m
//////////////////////////////////////////////////////////////////////

SetColumns(0); SetMemoryLimit(4*10^9);
Q := Rationals(); Px<x> := PolynomialRing(Q);
fD := -3*x^6 + 24*x^3 - 75;
CD := HyperellipticCurve(fD);
JD := Jacobian(CD);
pts := Points(JD : Bound := 2000);
gens := [P : P in pts | Order(P) eq 0];
g0 := gens[1];
printf "D : y^2 = %o\n", fD;
printf "generator g0 = (%o, %o)\n", g0[1], g0[2];

// the orders ord(W mod q) measured by the sieve (results/ap3/ap3_q*.txt)
QS := [ 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
        79, 83, 89, 97, 101, 103, 107, 113, 127, 131, 137, 139, 149, 151 ];
print "q   #J(D)(F_q)   ord(g0 mod q)   cofactor";
for q in QS do
  if Discriminant(HyperellipticCurve(PolynomialRing(GF(q))![GF(q)!c : c in Coefficients(fD)])) eq 0 then
    printf "%-6o bad reduction\n", q; continue;
  end if;
  kq := GF(q); Pq := PolynomialRing(kq);
  Cq := HyperellipticCurve(Pq![kq!c : c in Coefficients(fD)]);
  Jq := Jacobian(Cq);
  h := #Jq;
  aq := Pq![kq!c : c in Coefficients(g0[1])];
  bq := Pq![kq!c : c in Coefficients(g0[2])];
  gq := elt< Jq | aq, bq >;
  o := Order(gq);
  printf "%-6o %-12o %-12o %-12o\n", q, h, o, h div o;
end for;
print "SATCHECK_DONE";
quit;
