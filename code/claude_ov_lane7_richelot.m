// claude_ov_lane7_richelot.m -- Lane 7 (overnight 2026-07-25).
//
// DEPTH-1 RICHELOT AUDIT of all eleven contact-7 three-root [2,2,14] curves.
// The project's standing rule is that every new geometrically simple curve of
// rational torsion order >= 48 gets one; it has never been done for these
// (order 56).
//
// For each curve: enumerate the (2,2)-isogenous surfaces over Q
// (RichelotIsogenousSurfaces), and for each codomain that is a genus-2
// Jacobian compute the EXACT TorsionSubgroup.  Two things are being looked for:
//   (1) a neighbour with LARGER or DIFFERENT torsion -- a free new realization;
//   (2) a neighbour that is a PRODUCT of elliptic curves, which would
//       contradict geometric simplicity and so is a consistency check on the
//       End=Z certificates (a (2,2)-isogeny preserves End^0 up to isogeny).
//
// Factor type is [1,1,1,2] (deg-5 model: three rational finite Weierstrass
// points, one rational point at infinity, one irreducible quadratic), so the
// six Weierstrass points admit exactly three Galois-stable partitions into
// three pairs -- {r_i, r_j}, {r_k, inf}, {quadratic} -- hence three expected
// rational Richelot isogenies per curve.  A sextic model is built by moving one
// finite rational Weierstrass point to infinity is NOT needed: Magma handles
// the degree-5 model directly.
//
// POSITIVE CONTROL: the [2,2,2,8] witness x(x+1)(x+55^2)(x+99^2)(x+125^2)
// (2-rank 4) must produce a nonempty neighbour list.
//
// Markers: RICH / NBR / RICHCONTROL / LANE7_RICHELOT_DONE
SetColumns(0);
SetMemoryLimit(3*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q); Z := Integers();

Gfun := func< v | -(v^5 - v^3 - v^2/2)/(v+1)^2 >;
trips := [
  [-10, -10/7, -1/2], [-5, -15/8, -15/22], [-3, -3/4, -3/5],
  [-15/8, -15/19, -1/2], [-5/18, -10/49, 4/17], [-4/9, -4/25, 4/17],
  [-511/61, -511/625, -1/2], [-165/41, -33/16, -165/289],
  [-164/297, -1/2, 164/361], [-17/50, -34/189, 34/121], [-1/2, -13/49, 13/50]
];

procedure Audit(F, lbl)
  C := HyperellipticCurve(F);
  J := Jacobian(C);
  ok := true; RS := [];
  try
    RS := RichelotIsogenousSurfaces(J);
  catch ee
    ok := false;
    printf "RICH %o : RichelotIsogenousSurfaces FAILED: %o\n", lbl, ee`Object;
  end try;
  if not ok then return; end if;
  printf "RICH %o : %o rational (2,2)-neighbours\n", lbl, #RS;
  for i in [1..#RS] do
    A := RS[i];
    tp := Type(A);
    if tp eq SetCart then
      printf "NBR %o [%o] PRODUCT OF ELLIPTIC CURVES (would contradict geometric simplicity): %o\n", lbl, i, A;
      continue;
    end if;
    // Magma returns the neighbour on a model with RATIONAL coefficients;
    // TorsionSubgroup needs an integral one.  y^2 = g  ->  (Ly)^2 = L^2 g with
    // L the lcm of the coefficient denominators, then reduce to the minimal
    // Weierstrass model so the torsion computation is tractable.
    inv := []; okt := true; Nb := 0; emsg := "";
    try
      Cn := Curve(A);
      gn, hn := HyperellipticPolynomials(Cn);
      gfull := 4*gn + hn^2;                       // y^2 = gfull
      L := LCM([Denominator(co) : co in Coefficients(gfull)]);
      Gi := Parent(gfull) ! [ co*L^2 : co in Coefficients(gfull) ];
      Cint := HyperellipticCurve(Gi);
      Cmin := ReducedMinimalWeierstrassModel(Cint);
      gm, hm := HyperellipticPolynomials(Cmin);
      Cv := HyperellipticCurve(4*gm + hm^2);
      inv := Invariants(TorsionSubgroup(Jacobian(Cv)));
      Nb := Z!AbsoluteValue(Discriminant(Cmin));
    catch ee
      okt := false; emsg := Sprintf("%o", ee`Object);
    end try;
    if okt then
      printf "NBR %o [%o] torsion = %o   order = %o\n", lbl, i, inv,
             (#inv eq 0 select 1 else &*inv);
      printf "     minmodel y^2 + (%o)*y = %o\n", hm, gm;
    else
      printf "NBR %o [%o] torsion = ERROR (%o)\n", lbl, i, emsg;
    end if;
  end for;
end procedure;

print "==== POSITIVE CONTROL ====";
Audit(x*(x+1)*(x+55^2)*(x+99^2)*(x+125^2), "CONTROL2228");

print "==== THE ELEVEN [2,2,14] CURVES ====";
for T in trips do
  s := T[1]; t := T[2];
  c4 := (Gfun(s)-Gfun(t))/(s^2-t^2); c0 := Gfun(s) - c4*s^2;
  b := c4 - 2; a := 9/2 - c0 - c4;
  h := 1 - (7/2)*x + a*x^2 + b*x^3;
  f := (h^2 + (x-1)^7) div x^2;
  den := LCM([Denominator(co) : co in Coefficients(f)]);
  F := P![ co*den^2 : co in Coefficients(f) ];
  Audit(F, Sprintf("(%o,%o,%o)", T[1], T[2], T[3]));
end for;

print "LANE7_RICHELOT_DONE";
quit;
