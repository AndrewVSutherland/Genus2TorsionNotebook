//////////////////////////////////////////////////////////////////////
// claude_ov_m612ap_satindex.m    lane 9 ([6,12])   2026-07-25
//
// Is g0 = (x^2+2x+4, 5x+5) a GENERATOR of J(D)(Q), D : y^2 = -3x^6+24x^3-75?
// (rank 1, trivial torsion).  If g0 = n*G with n >= 2 then
//      hhat(G) = hhat(g0)/n^2 <= hhat(g0)/4,
// so G has naive height at most hhat(g0)/4 + c, where c = HeightConstant(J)
// bounds |h - hhat|.  Enumerating all points of naive height up to that bound
// and finding none of infinite order below g0 proves n = 1.
//
// This matters because the Abel-Prym MW sieve assumes W generates Pr(Q)
// modulo torsion, and claude_ov_m612ap_satcheck.m shows
//      ord(W mod q) = ord(g0 mod q)   at every prime tested.
//
// Usage: code/claude_magma_slot.sh -b code/claude_ov_m612ap_satindex.m
//////////////////////////////////////////////////////////////////////

SetColumns(0); SetMemoryLimit(6*10^9);
Q := Rationals(); Px<x> := PolynomialRing(Q);
fD := -3*x^6 + 24*x^3 - 75;
CD := HyperellipticCurve(fD);
JD := Jacobian(CD);
printf "D : y^2 = %o\n", fD;
rlo, rhi := RankBounds(JD);
printf "RankBounds = %o..%o\n", rlo, rhi;
T := TorsionSubgroup(JD);
printf "torsion invariants = %o\n", Invariants(T);

pts := Points(JD : Bound := 2000);
gens := [P : P in pts | Order(P) eq 0];
g0 := gens[1];
printf "g0 = (%o, %o)\n", g0[1], g0[2];

// ---- validate Magma's search-bound convention BEFORE relying on it -----
// Points(J : Bound := B) must enumerate exactly the points with
// exp(NaiveHeight) <= B.  g0 has exp(NaiveHeight(g0)) = 73 on the nose, so
// it must appear at Bound = 73 and not at Bound = 72.
nhg := NaiveHeight(g0);
printf "exp(NaiveHeight(g0)) = %o\n", Exp(nhg);
for B in [72, 73] do
  S0 := Points(JD : Bound := B);
  printf "  Bound = %-4o : #points = %-3o  g0 found : %o\n", B, #S0, g0 in S0;
end for;

c := HeightConstant(JD : Effort := 2);
printf "HeightConstant (|h - hhat| <= c) : c = %o\n", c;
h0 := Height(g0);
printf "hhat(g0) = %o\n", h0;
nh0 := NaiveHeight(g0);
printf "h(g0)    = %o\n", nh0;

// a candidate G with g0 = n*G, n >= 2, has hhat(G) <= hhat(g0)/4
for n in [2..6] do
  bnd := h0/n^2 + c;
  printf "n = %o : hhat(G) <= %o, so h(G) <= %o, search bound exp = %o\n",
     n, h0/n^2, bnd, Exp(bnd);
end for;

bnd := h0/4 + c;
B := Ceiling(Exp(bnd));
printf "searching J(D)(Q) for ALL points of naive height <= %o (Bound := %o)\n", bnd, B;
S := Points(JD : Bound := B);
printf "points found: %o\n", #S;
inf := [P : P in S | Order(P) eq 0];
printf "of infinite order: %o\n", #inf;
mn := -1; best := 0;
for P in inf do
  hp := Height(P);
  if mn lt 0 or hp lt mn then mn := hp; best := P; end if;
end for;
printf "smallest canonical height among them: %o\n", mn;
if #inf gt 0 then printf "attained at (%o, %o)\n", best[1], best[2]; end if;
printf "ratio hhat(g0)/min = %o\n", h0/mn;
if mn ge h0 then
  print "CONCLUSION: no point of J(D)(Q) has canonical height below hhat(g0),";
  print "  in particular none with hhat <= hhat(g0)/4, so g0 = n*G forces n = 1:";
  print "  g0 GENERATES J(D)(Q) = Z.";
else
  print "CONCLUSION: a smaller point exists -- g0 is NOT a generator.";
end if;
print "SATINDEX_DONE";
quit;
