// Independent verification of the (2,2,2,12) realization (order 96).
// Curve: y^2 = prod_{i=1..5}(A_i + B_i x) on the M(2,2,2,6) chart,
// (s,m,n) = (336396, -689185, -166464).  Run: magma -b claude_prod_02_22212_verify.m
// Expected output: TORSION INVARIANTS [ 2, 2, 2, 12 ] (order 96) and
// passing simplicity certificates at p = 37, 73, 113, 149
// (Frobenius charpoly irreducible of degree 4 AND its 12th-power transform
// irreducible of degree 4).  QM is independently excluded by
// Laga-Schembri-Shnidman-Voight (arXiv:2308.15193): PQM Jacobians have #tors <= 16.
QQ := Rationals();
P<x> := PolynomialRing(QQ);
A := [1,1,1,2,2];
B := [282322361376, -8243383980, -64241207724, -114724491840, 561915878400];
f := &*[A[i] + B[i]*x : i in [1..5]];
assert Degree(f) eq 5 and Discriminant(f) ne 0;
C := HyperellipticCurve(f);
J := Jacobian(C);
T := TorsionSubgroup(J);
printf "TORSION INVARIANTS: %o  (order %o)\n", Invariants(T), #T;
assert Invariants(T) eq [2,2,2,12];
D := Integers()!Discriminant(f);
ncert := 0;
for p in [37, 73, 113, 149] do
  if D mod p eq 0 then printf "p=%o bad, skipped\n", p; continue; end if;
  chi := EulerFactor(Jacobian(ChangeRing(C, GF(p))));
  chirev := P!Reverse(Coefficients(chi));
  if not IsIrreducible(chirev) then printf "p=%o: chi reducible, skipped\n", p; continue; end if;
  K<a> := NumberField(chirev);
  chi12 := MinimalPolynomial(a^12);
  ok := IsIrreducible(chi12) and Degree(chi12) eq 4;
  printf "p=%o: chi irreducible; 12th-power transform deg %o irreducible %o\n",
         p, Degree(chi12), ok;
  if ok then ncert +:= 1; end if;
end for;
printf "geometric-simplicity certificates passed at %o primes\n", ncert;
assert ncert ge 2;
printf "VERIFIED: (2,2,2,12) on a geometrically simple genus-2 Jacobian over Q\n";
quit;
