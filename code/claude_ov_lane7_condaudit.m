// claude_ov_lane7_condaudit.m -- Lane 7 (overnight 2026-07-25).
//
// CONDUCTOR AUDIT of all eleven contact-7 three-root curves.
//
// WHY.  The 2026-07-23 note records the five generic [2,2,14] conductors from
// PARI genus2red and concludes "all conductors are odd non-squares (no GL2-type
// smell)".  This lane found that genus2red is WRONG AT 2 on this family: for
// the sixth (RM) curve it returns N = 38025 = 195^2 with no 2-part, while
//   * its own 2-minimal model has v_2(Delta) = 7 != 0 (bad reduction at 2), and
//   * the curve is provably isogenous to the level-390 newform with Hecke field
//     Q(sqrt2) whose a_2 = 1 != 0, i.e. MULTIPLICATIVE reduction at 2, so
//     v_2(N) = 2 and N = 152100 = 390^2 (= Magma's answer).
// So every "odd conductor" in that note is suspect at 2.  This script measures,
// for each of the eleven curves, the one thing that settles the question
// rigorously: v_2 of the discriminant of the reduced minimal Weierstrass model.
// v_2(Delta_min) > 0  <=>  BAD reduction at 2  <=>  2 | N.
//
// It also records Magma's Conductor and whether Magma's own correctness
// disclaimer applies (it prints the Ogg warning exactly when v_2(D) >= 12), so
// the report can say precisely which 2-parts are certain and which are not.
//
// Forked: one child per curve, results/claude_ov_lane7_condaudit_<i>.txt.
// Markers: CA / LANE7_CONDAUDIT_DONE
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
// PARI genus2red conductors, from results/claude_ov_lane7_conductors.log
gpN := [ 9550095752925, 17934836205825, 38025, 14439206077220925,
         37247776656021702225, 21161909161575,
         1363558168459661985137828367867185325,
         11352674691759169446159049425,
         10200972432773189299959691439475,
         3069529451653112094922275, 69195300331841508825 ];

outfile := func< i | Sprintf("results/claude_ov_lane7_condaudit_%o.txt", i) >;

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
  C := HyperellipticCurve(F);
  Cm := ReducedMinimalWeierstrassModel(C);
  fm, hm := HyperellipticPolynomials(Cm);
  Dm := Z!Discriminant(Cm);
  v2 := Valuation(Dm, 2);
  txt := Sprintf("CA %o (%o, %o, %o) v_2(disc_min) = %o  bad_at_2 = %o  magma_ogg_disclaimer_applies = %o\n",
                 i, s, t, u, v2, v2 gt 0, v2 ge 12);
  NN := 0;
  try
    NN := Conductor(Cm);
  catch ee
    NN := -1;
  end try;
  if NN gt 0 then
    txt cat:= Sprintf("   magma N   = %o = %o\n", NN, Factorization(NN));
    txt cat:= Sprintf("   pari  N   = %o\n", gpN[i]);
    txt cat:= Sprintf("   ratio magma/pari = %o ; odd parts agree = %o\n",
                      NN/gpN[i], (NN div 2^Valuation(NN,2)) eq (gpN[i] div 2^Valuation(gpN[i],2)));
  else
    txt cat:= Sprintf("   magma Conductor FAILED (memory/other); pari N = %o\n", gpN[i]);
  end if;
  Write(outfile(i), txt : Overwrite := true);
  quit;
end for;

WaitForAllChildren();
for i in [1..#trips] do
  ok := true;
  try
    txt := Read(outfile(i));
  catch ee
    ok := false;
  end try;
  if ok then printf "%o", txt; else printf "CA %o : child produced no output (died)\n", i; end if;
end for;
print "LANE7_CONDAUDIT_DONE";
quit;
