// claude_ov_lane7_richelot2.m -- Lane 7 (overnight 2026-07-25).
//
// DEPTH-1 RICHELOT AUDIT, corrected.  The first pass
// (code/claude_ov_lane7_richelot.m) found exactly 3 rational (2,2)-neighbours
// for each of the eleven contact-7 three-root curves -- as predicted, since the
// quintic factor type is [1,1,1,2], so the six Weierstrass points
// {r1, r2, r3, inf, the conjugate pair} admit exactly three Galois-stable
// partitions into three pairs -- but TorsionSubgroup then failed on almost every
// neighbour, because Magma hands back the codomain as a NON-INTEGRAL sextic with
// enormous coefficients.
//
// Fix: clear denominators (y -> L*y) and pass through
// ReducedMinimalWeierstrassModel + SimplifiedModel before the exact torsion
// call.  Forked, one child per curve.
//
// What we are looking for: a neighbour with LARGER 2-part -- [2,2,28], [2,4,14],
// [4,28] -- which would be a free order-112 realization.  A (2,2)-isogeny fixes
// the odd part (so the 7 survives) and preserves End^0, so every neighbour of a
// certified End=Z curve is again geometrically simple; a codomain that came back
// as a PRODUCT of elliptic curves would contradict the certificates and is
// flagged.
//
// Markers: RICH2 / NBR2 / LANE7_RICHELOT2_DONE
SetColumns(0);
Q := Rationals(); P<x> := PolynomialRing(Q); Z := Integers();
if not assigned MemGB then MemGB := 6; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;

Gfun := func< v | -(v^5 - v^3 - v^2/2)/(v+1)^2 >;
trips := [
  [-10, -10/7, -1/2], [-5, -15/8, -15/22], [-3, -3/4, -3/5],
  [-15/8, -15/19, -1/2], [-5/18, -10/49, 4/17], [-4/9, -4/25, 4/17],
  [-511/61, -511/625, -1/2], [-165/41, -33/16, -165/289],
  [-164/297, -1/2, 164/361], [-17/50, -34/189, 34/121], [-1/2, -13/49, 13/50]
];

outfile := func< i | Sprintf("results/claude_ov_lane7_richelot2_%o.txt", i) >;

for i in [1..#trips] do
  pid := Fork();
  if pid ne 0 then continue; end if;
  // ---------------- child ----------------
  SetMemoryLimit(MemGB*10^9);
  s := trips[i][1]; t := trips[i][2]; u := trips[i][3];
  c4 := (Gfun(s)-Gfun(t))/(s^2-t^2); c0 := Gfun(s) - c4*s^2;
  b := c4 - 2; a := 9/2 - c0 - c4;
  h := 1 - (7/2)*x + a*x^2 + b*x^3;
  f := (h^2 + (x-1)^7) div x^2;
  den := LCM([Denominator(co) : co in Coefficients(f)]);
  F := P![ co*den^2 : co in Coefficients(f) ];
  J := Jacobian(HyperellipticCurve(F));
  txt := "";
  RS := []; okR := true;
  try
    RS := RichelotIsogenousSurfaces(J);
  catch ee
    okR := false;
    txt cat:= Sprintf("RICH2 %o (%o, %o, %o) RichelotIsogenousSurfaces FAILED: %o\n", i, s, t, u, ee`Object);
  end try;
  if okR then
    txt cat:= Sprintf("RICH2 %o (%o, %o, %o) source torsion [2,2,14] : %o rational (2,2)-neighbours\n", i, s, t, u, #RS);
    for k in [1..#RS] do
      A := RS[k];
      if Type(A) ne JacHyp then
        txt cat:= Sprintf("NBR2 %o.%o PRODUCT-OF-ELLIPTIC-CURVES (type %o) -- CONTRADICTS End=Z\n", i, k, Type(A));
        continue;
      end if;
      inv := []; okt := true; lbl := "";
      try
        CA := Curve(A);
        g6 := HyperellipticPolynomials(SimplifiedModel(CA));
        L := LCM([Denominator(co) : co in Coefficients(g6)]);
        // y -> L*y clears denominators without changing the curve
        Gi := P![ Z!(co*L^2) : co in Coefficients(g6) ];
        Cm := ReducedMinimalWeierstrassModel(HyperellipticCurve(Gi));
        fm, hm := HyperellipticPolynomials(Cm);
        inv := Invariants(TorsionSubgroup(Jacobian(HyperellipticCurve(4*fm + hm^2))));
        lbl := Sprintf("y^2 + (%o)*y = %o", hm, fm);
      catch ee
        okt := false; lbl := Sprintf("%o", ee`Object);
      end try;
      if okt then
        txt cat:= Sprintf("NBR2 %o.%o torsion = %o   order = %o\n", i, k, inv,
                          #inv eq 0 select 1 else &*[Z|e : e in inv]);
        txt cat:= Sprintf("     minmodel %o\n", lbl);
      else
        txt cat:= Sprintf("NBR2 %o.%o torsion FAILED: %o\n", i, k, lbl);
      end if;
    end for;
  end if;
  Write(outfile(i), txt : Overwrite := true);
  quit;
end for;

WaitForAllChildren();
for i in [1..#trips] do
  printf "%o", Read(outfile(i));
end for;
print "LANE7_RICHELOT2_DONE";
quit;
