// claude_z96_vertical_check.m — test the vertical components of the 3-contact
// elimination on the Elkies-32 family: rational u-values from the GCD factors.
// Any NONSINGULAR member with 3 | #torsion is a cyclic Z/96 curve.
SetColumns(0);
SetMemoryLimit(6*10^9);
load "code/claude_z96_family_setup.m";
Ku, Px, ftil, fu, zu, ru := BuildFamily();
Q := Rationals(); Pq<x> := PolynomialRing(Q);
cfs := [ Numerator(Coefficient(ftil, i)) : i in [0..5] ];

// rational candidates: the linear GFACTORs + rational roots of the deg-2/deg-4 ones
cands := [Q| 6120, 41820/7, 41310/7, 39780/7 ];
P1<T> := PolynomialRing(Q);
for pol in [ T^2 - 83232/7*T + 1732266000/49,
             T^2 - 82620/7*T + 1707296400/49,
             T^4 - 1814580/77*T^3 + 112260200400/539*T^2 - 440992412856000/539*T + 31834731407931360000/26411 ] do
  for r in Roots(pol) do Append(~cands, r[1]); end for;
end for;
printf "rational candidate u-values: %o\n", cands;

for u0 in cands do
  fv := Pq![ Evaluate(c, u0) : c in cfs ];
  if Degree(fv) lt 5 then printf "u=%o: DEGENERATE (degree %o)\n", u0, Degree(fv); continue; end if;
  if Discriminant(fv) eq 0 then printf "u=%o: DEGENERATE (disc 0)\n", u0; continue; end if;
  den := LCM([Denominator(c) : c in Coefficients(fv)]);
  fint := Pq![ c*den^2 : c in Coefficients(fv) ];
  Cs := SimplifiedModel(ReducedMinimalWeierstrassModel(HyperellipticCurve(fint)));
  T0 := TorsionSubgroup(Jacobian(Cs));
  inv := Invariants(T0);
  printf "u=%o: NONSINGULAR member, torsion %o%o\n", u0, inv,
    (#T0 mod 3 eq 0) select "  <<< 3-PART PRESENT: Z/96 CANDIDATE" else "";
end for;
print "ALL_DONE";
quit;
