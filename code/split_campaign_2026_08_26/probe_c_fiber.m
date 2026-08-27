// probe: structure of Z(3,2) fibered by j1 -- is the fiber a parametrizable conic?
SetColumns(0);
SetMemoryLimit(8*10^9);
AttachSpec("/home/claude/Magma/magma.spec");
Attach("/home/claude/torsion_jac/product/N-congruences/ZNr-equations.m");

jj, jj1728 := WNrModuli(3,2);
P := Parent(jj);
printf "parent: %o\n", P;
jsum := -(jj1728 - jj - 1728^2)/1728;
nu := Numerator(jj); de := Denominator(jj);
printf "jj: num deg %o, den deg %o\n", [Degree(nu, i) : i in [1..2]], [Degree(de, i) : i in [1..2]];
ns := Numerator(jsum); ds := Denominator(jsum);
printf "jsum: num deg %o, den deg %o\n", [Degree(ns, i) : i in [1..2]], [Degree(ds, i) : i in [1..2]];

// sample fiber over a specific rational j0: j of X1(8) member t=3
Qt := Rationals();
t := Qt!3;
bv := (2*t-1)*(t-1); cv := (2*t-1)*(t-1)/t;
E8 := EllipticCurve([1-cv, -bv, -bv, 0, 0]);
j0 := jInvariant(E8);
printf "j0 = %o\n", j0;

// fiber curve: j0^2 - jsum*j0 + jj = 0  (as plane curve in (u,v))
A2<uu,vv> := AffineSpace(Rationals(), 2);
R2 := CoordinateRing(A2);
h := hom< Parent(nu) -> R2 | [uu, vv] >;
num := h(Numerator(j0^2 - jsum*j0 + jj));
printf "fiber polynomial degree: %o (u-deg %o, v-deg %o)\n", TotalDegree(num), Degree(num, uu), Degree(num, vv);
F := Curve(A2, num);
printf "fiber is irreducible: %o\n", IsIrreducible(num);
comps := IrreducibleComponents(F);
printf "components: %o\n", #comps;
for i in [1..#comps] do
    Ci := Curve(comps[i]);
    printf "  comp %o: degree %o, geometric genus %o\n", i, Degree(comps[i]), GeometricGenus(Ci);
end for;
quit;
