// Independent verification of (2,2,2,12) curve #2 (2026-07-18 sibling hunt).
// (s,m,n)=(2208,-8303,-7200) on M(2,2,2,6). Run: magma -b claude_sib_curve2_verify.m
// Expect: TORSION [2,2,2,12] order 96; certs pass at p=71,103,127,137.
QQ := Rationals(); P<x> := PolynomialRing(QQ);
A := [1,1,1,2,2];
B := [25648128,-36568896,-52466496,-59781600,23309856];
f := &*[A[i] + B[i]*x : i in [1..5]];
assert Degree(f) eq 5 and Discriminant(f) ne 0;
C := HyperellipticCurve(f); J := Jacobian(C);
T := TorsionSubgroup(J);
printf "TORSION: %o (order %o)\n", Invariants(T), #T;
assert Invariants(T) eq [2,2,2,12];
D := Integers()!Discriminant(f); nc := 0;
for p in [71,103,127,137,149,157] do
  if D mod p eq 0 then continue; end if;
  chi := P!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(C,GF(p))))));
  if not IsIrreducible(chi) then printf "p=%o chi reducible, skip\n", p; continue; end if;
  K<a> := NumberField(chi); c12 := MinimalPolynomial(a^12);
  ok := IsIrreducible(c12) and Degree(c12) eq 4;
  printf "p=%o: cert %o\n", p, ok; if ok then nc +:= 1; end if;
end for;
printf "certificates: %o\n", nc;
assert nc ge 2;
printf "VERIFIED: (2,2,2,12) on a geometrically simple genus-2 Jacobian over Q (curve #2)\n";
printf "G2: %o\n", G2Invariants(C);
quit;
