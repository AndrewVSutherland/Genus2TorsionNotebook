// claude_2214_minimal_models.m — minimal Weierstrass models for the five
// verified generic [2,2,14] curves (three-root route); re-verify exact torsion
// and End=Z on the minimal model; print y^2 + h*y = f forms and discriminants.
SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q);
Afun := func< w | (2*w^5+4*w^4+6*w^3+8*w^2+10*w+5)/(2*(w+1)^2) >;

function EndZ(C, ord)
  D := Integers()!Discriminant(C); p0 := 0; chip := P!0; dP := 0; nf := 0;
  for p in PrimesInInterval(11, 400) do
    if D mod p eq 0 then continue; end if;
    chi := P!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(C, GF(p))))));
    if Degree(chi) ne 4 or not IsIrreducible(chi) then continue; end if;
    K<aa> := NumberField(chi); st := true;
    for n in [2..12] do if Degree(MinimalPolynomial(aa^n)) ne 4 then st := false; break; end if; end for;
    dL := Degree(SplittingField(chi));
    if st and p0 eq 0 then p0 := p; chip := chi; dP := dL; continue; end if;
    if p0 ne 0 and st and Degree(SplittingField(chip*chi)) eq dP*dL then return (ord gt 18), p0, p, chip, chi; end if;
    nf +:= 1; if nf ge 20 then break; end if;
  end for;
  return false, 0, 0, P!0, P!0;
end function;

cands := [ [-10,-1/2,-10/7], [-5,-15/8,-15/22], [-1/2,-15/8,-15/19],
           [-4/9,4/17,-4/25], [4/17,-5/18,-10/49] ];
for c in cands do
  s := c[1]; t := c[2]; u := c[3];
  b := (Afun(t)-Afun(s))/(s^2-t^2);
  a := Afun(s) - b*(1-s^2);
  h := 1 - (7/2)*x + a*x^2 + b*x^3;
  f := (h^2 + (x-1)^7) div x^2;
  den := LCM([Denominator(co) : co in Coefficients(f)]);
  Fi := P![ co*den^2 : co in Coefficients(f) ];
  C := HyperellipticCurve(Fi);
  Cm := ReducedMinimalWeierstrassModel(C);
  fm, hm := HyperellipticPolynomials(Cm);
  printf "STU %o %o %o\n", s, t, u;
  printf "  MINMODEL y^2 + (%o)*y = %o\n", hm, fm;
  D := Integers()!Discriminant(Cm);
  printf "  |disc| = %o = %o\n", AbsoluteValue(D), Factorization(AbsoluteValue(D));
  Fm := 4*fm + hm^2;
  Cv := HyperellipticCurve(Fm);
  J := Jacobian(Cv);
  inv := Invariants(TorsionSubgroup(J));
  ez, p0, q0, cp, cq := EndZ(Cv, 56);
  printf "  torsion=%o EndZ=%o p0=%o q0=%o\n  chi_p0=%o\n  chi_q0=%o\n", inv, ez, p0, q0, cp, cq;
end for;
print "MINMODELS_DONE";
quit;
