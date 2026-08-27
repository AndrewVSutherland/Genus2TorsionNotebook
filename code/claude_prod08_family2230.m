// (2,30) gate: clean g-parametrization of the [30] family; the (2,3) locus CQ;
// the (1,2,2) cut C122; genera and rational points.
SetColumns(0);
Qx<x> := PolynomialRing(Rationals());

// ---- clean parametrization over Q(g) ----
Kg<g> := FunctionField(Rationals());
t := g/(1+g^2);
w := (1-g^2)/(1+g^2);          // w^2 = 1 - 4t^2
assert w^2 eq 1 - 4*t^2;
A := -4*t^2*(4*t^2-1)^2;
B := -96*t^4 + 16*t^2 + 2;
C0 := -(4*t^2+3)^2;
sq := 2*(4*t^2-1)^2*w;         // sqrt of B^2-4AC
assert sq^2 eq B^2 - 4*A*C0;
r := (-B + sq)/(2*A);
c := t*r;
Kx<X> := PolynomialRing(Kg);
Fint := 2*X^5 - (4*c^2+9)*X^4 + 16*c^2*X^3 + (32*c^4-24*c^2)*X^2 - 96*c^4*X - 64*c^6 - 16*c^4;
printf "check Fint(r)=0: %o\n", Evaluate(Fint, r) eq 0;
printf "r(g) = %o\n", r;
printf "c(g) = %o\n", c;
// sanity anchor g=3 -> (45/4, 27/8)?
printf "at g=3: r=%o c=%o\n", Evaluate(r,3), Evaluate(c,3);
G4 := Fint div (X - r);
printf "G4 = %o\n", G4;
printf "G4 factorization over Q(g): %o\n", [<Degree(f[1]),f[2]> : f in Factorization(G4)];

// ---- C122: {(a,b,g): X^2+aX+b divides G4} ----
print "=== C122 (the (1,2,2) cut => (2,30)) ===";
A3<aa,bb,gg> := AffineSpace(Rationals(), 3);
FA3 := FieldOfFractions(CoordinateRing(A3));
PY<Y> := PolynomialRing(FA3);
tt := gg/(1+gg^2);
ww := (1-gg^2)/(1+gg^2);
AA := -4*tt^2*(4*tt^2-1)^2;
BB := -96*tt^4 + 16*tt^2 + 2;
CC := -(4*tt^2+3)^2;
sq2 := 2*(4*tt^2-1)^2*ww;
rr := (-BB + sq2)/(2*AA);
ccc := tt*rr;
FiY := 2*Y^5 - (4*ccc^2+9)*Y^4 + 16*ccc^2*Y^3 + (32*ccc^4-24*ccc^2)*Y^2 - 96*ccc^4*Y - 64*ccc^6 - 16*ccc^4;
G4Y := FiY div (Y - rr);
quad := Y^2 + FA3!aa*Y + FA3!bb;
rem := G4Y mod quad;
e1 := Numerator(Coefficient(rem,1));
e0 := Numerator(Coefficient(rem,0));
S122 := Scheme(A3, [e1, e0]);
print "dim S122:", Dimension(S122);
comps := IrreducibleComponents(S122);
printf "%o components\n", #comps;
for co in comps do
  if Dimension(co) ne 1 then printf "  comp dim %o (skip)\n", Dimension(co); continue; end if;
  Cco := Curve(co);
  printf "  comp deg %o", Degree(Cco);
  gen := Genus(Cco);
  printf " genus %o\n", gen;
  if gen le 2 then
    ptsx := PointSearch(Cco, 100000);
    printf "    points: %o\n", ptsx;
  end if;
end for;
quit;
