// Same-c correspondence: c(g') = ±c(g), g' != g: each rational point => two rational
// roots of the same quintic F => 2-rank 2 => [2,30] candidate.
SetColumns(0);
A2<u,v> := AffineSpace(Rationals(), 2);   // u = g, v = g'
cn := func<z | z*(z^2+3)^2>;
cd := func<z | 2*(z^2-1)^2>;
// c(u) = cn(u)/cd(u). same c: cn(u)*cd(v) - cn(v)*cd(u) = 0; contains v=u.
P := cn(u)*cd(v) - cn(v)*cd(u);
P1 := P div (v - u);
print "minus-branch: degree", TotalDegree(P1);
Cm := Curve(A2, P1);
print "irreducible components (minus branch):";
for co in IrreducibleComponents(Cm) do
  Cco := Curve(co);
  gen := Genus(Cco);
  printf "  genus %o\n", gen;
  pts := PointSearch(Cco, 100000);
  n := 0;
  for pt in pts do
    if pt[1] ne pt[2] and Abs(pt[1]) ne 1 and Abs(pt[2]) ne 1 then
      printf "    NONTRIVIAL PT %o  <== [2,30] CANDIDATE\n", pt; n +:= 1;
    end if;
  end for;
  printf "    (%o points total, %o nontrivial)\n", #pts, n;
end for;
// plus branch: c(v) = -c(u): cn(u)*cd(v) + cn(v)*cd(u) = 0; contains v=-u (since cn odd, cd even)
Q := cn(u)*cd(v) + cn(v)*cd(u);
Q1 := Q div (v + u);
print "plus-branch: degree", TotalDegree(Q1);
Cq := Curve(A2, Q1);
print "irreducible components (plus branch):";
for co in IrreducibleComponents(Cq) do
  Cco := Curve(co);
  gen := Genus(Cco);
  printf "  genus %o\n", gen;
  pts := PointSearch(Cco, 100000);
  n := 0;
  for pt in pts do
    if pt[1] ne -pt[2] and Abs(pt[1]) ne 1 and Abs(pt[2]) ne 1 then
      printf "    NONTRIVIAL PT %o  <== [2,30] CANDIDATE\n", pt; n +:= 1;
    end if;
  end for;
  printf "    (%o points total, %o nontrivial)\n", #pts, n;
end for;
quit;
