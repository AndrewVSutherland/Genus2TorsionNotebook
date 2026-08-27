// Lane 2 (claude_ov_c7z): the ORDER-112 (split-all) locus over the c=2 slice.
//
// c=2 slice pencil:  P_m(w) = w^4 - 4w^3 + m w^2 - (2/3) m w + 1,   m = M(w) = -(w^4-4w^3+1)/(w^2-(2/3)w).
//   C  = {(x,y): M(x)=M(y), x != y}   -- >=2 rational roots -> [2,2,14].  GENUS 1, rank 1.
//   C3 = {(x,y,z) distinct : M(x)=M(y)=M(z)}  -- 3 rational roots forces the 4th (e1=4),
//        so C3(Q) = the ORDER-112 points on this slice.  C3 -> C is the double cover
//        "residual quadratic of P_m splits", i.e. w^2 = s^2-4p with s=4-x-y, p=1/(xy).
// This script computes the genus of C3 (both as a function-field extension of C and
// directly as a space curve), and searches it for rational points.
SetColumns(0);
SetMemoryLimit(20*10^9);
Q := Rationals();

// ---------- 1. the base curve C ----------
A2<x,y> := AffineSpace(Q,2);
F := -x^3*y^2 + 2/3*x^3*y - x^2*y^3 + 14/3*x^2*y^2 - 8/3*x^2*y + 2/3*x*y^3 - 8/3*x*y^2 + x + y - 2/3;
C := Curve(A2,F);
PC := ProjectiveClosure(C);
printf "C: irreducible=%o genus=%o\n", IsIrreducible(F), Genus(PC);

// ---------- 2. the double cover, via the function field of C ----------
// FIXED 2026-07-25: FunctionField(C) is a FldFunFracSch, which has no ext<> constructor.
// Pass to the ALGORITHMIC function field (a FldFun) first; then FunctionField(poly) works.
KC := FunctionField(C);
FF, mp := AlgorithmicFunctionField(KC);
xx := mp(KC!(A2.1)); yy := mp(KC!(A2.2));
s := 4 - xx - yy;
p := 1/(xx*yy);
Dfn := s^2 - 4*p;
printf "D = s^2-4p as a function on C: is it a square in K(C)? %o\n", IsSquare(Dfn);
Pt<t> := PolynomialRing(FF);
L := FunctionField(t^2 - Dfn);
printf "COVER genus(C3) = %o   (degree over K(C) = %o)\n", Genus(L), Degree(L);
// also the quadratic twist structure: which quadratic extension?
printf "COVER ramified places of the extension: %o\n", #[ pl : pl in Support(Divisor(Dfn)) | IsOdd(Valuation(Dfn,pl)) ];

// ---------- 3. the same curve as an explicit space curve in A^3 ----------
A3<u1,u2,u3> := AffineSpace(Q,3);
u4 := 4 - u1 - u2 - u3;
e2 := u1*u2+u1*u3+u1*u4+u2*u3+u2*u4+u3*u4;
e3 := u1*u2*u3+u1*u2*u4+u1*u3*u4+u2*u3*u4;
e4 := u1*u2*u3*u4;
Z := Scheme(A3, [Numerator(e4-1), Numerator(3*e3-2*e2)]);
Zd := [ c : c in IrreducibleComponents(Z) ];
printf "split-all space curve: %o irreducible component(s), dims %o\n", #Zd, [Dimension(c) : c in Zd];
for i in [1..#Zd] do
  cpt := Zd[i];
  if Dimension(cpt) eq 1 then
    ok := true; g := -1;
    try g := Genus(Curve(cpt)); catch e ok := false; end try;
    printf "  component %o: degree=%o genus=%o (ok=%o)\n", i, Degree(cpt), g, ok;
  else
    printf "  component %o: dimension %o degree %o\n", i, Dimension(cpt), Degree(cpt);
  end if;
end for;

// ---------- 4. brute rational-point search on the split-all curve ----------
// parametrize by (x,y) on C with small height, test s^2-4p square
cnt := 0; found := 0;
for dn in [1..60] do
  for nu in [-240..240] do
    if GCD(Abs(nu),dn) ne 1 then continue; end if;
    xv := nu/dn;
    if xv eq 0 or xv eq 2/3 then continue; end if;
    fy := Evaluate(F, [xv, PolynomialRing(Q).1]);
    if fy eq 0 then continue; end if;
    for r in Roots(fy) do
      yv := r[1];
      if yv eq xv or yv eq 0 then continue; end if;
      cnt +:= 1;
      sv := 4 - xv - yv; pv := 1/(xv*yv);
      dv := sv^2 - 4*pv;
      isq, rt := IsSquare(dv);
      if isq and dv ne 0 then
        found +:= 1;
        printf "SPLITALL_C2 x=%o y=%o s=%o p=%o D=%o\n", xv, yv, sv, pv, dv;
      end if;
    end for;
  end for;
end for;
printf "C2 point search: %o rational points of C found, %o with D square\n", cnt, found;
printf "COVER_DONE\n";
quit;
