// Lane 2: the c = 2 slice of the contact-7 three-root surface is a genus-1 curve.
// Curve (from claude_ov_c7z_slice.m):
//   F = -x^3*y^2 + 2/3*x^3*y - x^2*y^3 + 14/3*x^2*y^2 - 8/3*x^2*y + 2/3*x*y^3 - 8/3*x*y^2 + x + y - 2/3
// (x,y) = the two extra rational z's; z1 = 2; residual quadratic w^2 - s w + p,
//   s = 4 - (x+y),  p = 1/(x*y).   Order-112 <=> s^2 - 4p is a rational SQUARE.
SetColumns(0);
Q := Rationals();
A2<x,y> := AffineSpace(Q,2);
F := -x^3*y^2 + 2/3*x^3*y - x^2*y^3 + 14/3*x^2*y^2 - 8/3*x^2*y + 2/3*x*y^3 - 8/3*x*y^2 + x + y - 2/3;
C := Curve(A2,F);
printf "affine curve is irreducible: %o\n", IsIrreducible(F);
PC := ProjectiveClosure(C);
printf "genus = %o\n", Genus(PC);

known := [ [-1/9,-7/3], [-8/7,19/4], [-61/450,625/114], [297/133,361/525], [49/36,50/63] ];
for pt in known do
  printf "  known point (%o,%o) on C: %o\n", pt[1],pt[2], Evaluate(F,[pt[1],pt[2]]) eq 0;
end for;

P0 := PC ! [known[1][1], known[1][2], 1];
E, mp := EllipticCurve(PC, P0);
E2, m2 := SimplifiedModel(E);
printf "EllipticCurve: %o\n", E2;
printf "Conductor = %o\n", Conductor(E2);
printf "j-invariant = %o\n", jInvariant(E2);
printf "TorsionSubgroup = %o\n", Invariants(TorsionSubgroup(E2));
r1, r2 := RankBounds(E2);
printf "RankBounds = %o .. %o\n", r1, r2;
G, mg, fl1, fl2 := MordellWeilGroup(E2);
printf "MordellWeilGroup = %o  (proven: %o %o)\n", Invariants(G), fl1, fl2;
gens := [ mg(g) : g in Generators(G) ];
printf "generators: %o\n", gens;

// images of the known points
printf "--- known points on E2 ---\n";
imgs := [];
for pt in known do
  pp := PC ! [pt[1],pt[2],1];
  im := m2(mp(pp));
  Append(~imgs, im);
  printf "  (%o,%o) -> %o\n", pt[1],pt[2], im;
end for;
// also the swapped points
for pt in known do
  pp := PC ! [pt[2],pt[1],1];
  im := m2(mp(pp));
  Append(~imgs, im);
end for;
printf "heights: %o\n", [ Height(p) : p in imgs ];
printf "ELL2_DONE\n";
quit;
