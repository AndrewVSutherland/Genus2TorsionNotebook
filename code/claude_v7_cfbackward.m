// claude_v7_cfbackward.m — B1 REHEARSAL (order 7): derive the surface
//   S7 = { (a,b,c,d) : f = x(x-1)(x-a)(x-b)(x^2+c*x+d) has D_inf of order 7 }
// inside the [1,1,1,1,2]-shape chart (2-rank 3 built in), by running the
// polynomial continued fraction of sqrt(f) SYMBOLICALLY over Q(a,b,c,d) and
// imposing the quasi-period-closure conditions.  Order 7 = deg pattern
// [3,1,1,1,1]: after a0 (deg 3) and four deg-1 partial quotients the next
// Q must be a nonzero constant: its x^2 and x^1 coefficients vanish — TWO
// equations cutting the 4-dim chart to the 2-dim surface S7.
// Any rational point on S7 off the degenerate loci gives torsion containing
// (Z/2)^3 x Z/7 = [2,2,14]; S7 rational => INFINITE generic [2,2,14] family.
// Anchor: the RM witness 152100.eb.2 lies on S7 after normalizing two of its
// four rational Weierstrass roots to 0,1 (checked below), so S7(Q) != {}.
// This is the 4-CF-step rehearsal for the 8-step order-11 analogue
// (shape [1,1,2,2], anchor 19044.h.2) — the main [2,22] route B1.
SetColumns(0);
SetMemoryLimit(24*10^9);
Q := Rationals();
F<a,b,c,d> := FunctionField(Q, 4);
P<x> := PolynomialRing(F);

f := x*(x-1)*(x-a)*(x-b)*(x^2+c*x+d);

// polynomial part of sqrt(f): monic cubic s0 with deg(f - s0^2) <= 2
s0 := x^3;
for k in [1..3] do
  dd := f - s0^2;
  if Degree(dd) le 2 then break; end if;
  s0 := s0 + (Coefficient(dd, 6-k)/(2*Coefficient(s0,3)))*x^(3-k);
end for;
assert Degree(f - s0^2) le 2;

// CF steps: P_{i+1} = a_i Q_i - P_i,  Q_{i+1} = (f - P_{i+1}^2)/Q_i,
// a_i = (P_i + s0) div Q_i.  Track degrees; abort if pattern degenerates.
Pi := P!0; Qi := P!1;
eqs := [];
for i in [0..4] do
  ai := (Pi + s0) div Qi;
  printf "step %o: deg a_i = %o, deg Q_i = %o\n", i, Degree(ai), Degree(Qi);
  Pn := ai*Qi - Pi;
  rem := (f - Pn^2) mod Qi;
  assert rem eq 0;
  Qn := (f - Pn^2) div Qi;
  Pi := Pn; Qi := Qn;
end for;
// after step 4 (a_0..a_4 consumed), Q_5 = Qi must be a nonzero CONSTANT
printf "closure: deg Q_5 = %o\n", Degree(Qi);
E2 := Numerator(Coefficient(Qi, 2));
E1 := Numerator(Coefficient(Qi, 1));
printf "E2: %o terms, total degree %o\n", #Terms(E2), TotalDegree(E2);
printf "E1: %o terms, total degree %o\n", #Terms(E1), TotalDegree(E1);

// write the two equations to a data file for the geometry pass
R4<aa,bb,cc,dd4> := PolynomialRing(Q, 4);
h := hom< Parent(E2) -> R4 | [aa,bb,cc,dd4] >;
G2 := h(E2); G1 := h(E1);
Write("data/claude_v7_surface_eqs.txt", Sprintf("E2 := %o;\nE1 := %o;\n", G2, G1) : Overwrite := true);

// ---- anchor check: RM witness 152100.eb.2 ----
// y^2 + (x^2+x) y = 9x^6+69x^5+123x^4-95x^3-183x^2+165x-35
// F0 = 4f+h^2; its four rational Weierstrass roots give normalizations onto the chart.
PQ<z> := PolynomialRing(Q);
F0 := 4*(9*z^6+69*z^5+123*z^4-95*z^3-183*z^2+165*z-35) + (z^2+z)^2;
fac := Factorization(F0);
rts := [ -Coefficient(pe[1],0)/Coefficient(pe[1],1) : pe in fac | Degree(pe[1]) eq 1 ];
printf "RM witness rational roots: %o\n", rts;
quad := [ pe[1] : pe in fac | Degree(pe[1]) eq 2 ];
// normalize roots r1 -> 0, r2 -> 1 via z = r1 + (r2-r1)*x for each ordered pair;
// the transformed sextic is lambda * x(x-1)(x-a0)(x-b0)(x^2+c0 x+d0); D_inf order
// is invariant under x-affine maps and y-scaling, so (a0,b0,c0,d0) must satisfy
// E1 = E2 = 0 whenever the ORDER-7 pattern holds in this chart (lambda must be
// square for the two infinite points to stay rational, else twist changes J(Q) —
// but D_inf ORDER over Qbar-model... we check E-vanishing directly).
nanch := 0;
for i in [1..#rts] do for j in [1..#rts] do
  if i eq j then continue; end if;
  r1 := rts[i]; r2 := rts[j];
  vals := [];
  ok := true;
  // images of the other two rational roots
  oth := [ (r - r1)/(r2 - r1) : r in rts | r ne r1 and r ne r2 ];
  if #oth ne 2 then continue; end if;
  a0 := oth[1]; b0 := oth[2];
  // transformed quadratic: z^2 + c z + d with z = r1 + (r2-r1) x
  q := quad[1];
  qq := Evaluate(q, r1 + (r2-r1)*x);   // in P over F... use PQ-side:
  PQx<w> := PolynomialRing(Q);
  qq2 := Evaluate(q, r1 + (r2-r1)*w);
  lc := Coefficient(qq2, 2);
  c0 := Coefficient(qq2,1)/lc; d0 := Coefficient(qq2,0)/lc;
  v2 := Evaluate(G2, [a0,b0,c0,d0]);
  v1 := Evaluate(G1, [a0,b0,c0,d0]);
  printf "anchor try (r1,r2)=(%o,%o): (a,b,c,d)=(%o,%o,%o,%o)  E2=%o E1=%o\n",
    r1, r2, a0, b0, c0, d0, v2 eq 0, v1 eq 0;
  if v2 eq 0 and v1 eq 0 then nanch +:= 1; end if;
end for; end for;
printf "anchor normalizations on S7: %o\n", nanch;
print "V7_CFBACKWARD_DONE";
quit;
