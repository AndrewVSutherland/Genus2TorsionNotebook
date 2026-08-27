// verify2214.m — EXACT verification of the four generic-[2,2,14] candidates
// found by the contact-7 three-root sweep (probeA).  For each (s,t,u):
//   rebuild (a,b), h, f from scratch; assert the three roots; exact
//   TorsionSubgroup; marked-class order; End=Z certificate
//   (root-power geom-simple + disjoint splitting fields + LSSV ord>18).
SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q);

Afun := func< w | (2*w^5+4*w^4+6*w^3+8*w^2+10*w+5)/(2*(w+1)^2) >;

function EndZ(C, ord)
  D := Integers()!Discriminant(C); p0 := 0; chip := P!0; dP := 0; nf := 0;
  for p in PrimesInInterval(11, 300) do
    if D mod p eq 0 then continue; end if;
    chi := P!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(C, GF(p))))));
    if Degree(chi) ne 4 or not IsIrreducible(chi) then continue; end if;
    K<aa> := NumberField(chi); st := true;
    for n in [2..12] do if Degree(MinimalPolynomial(aa^n)) ne 4 then st := false; break; end if; end for;
    dL := Degree(SplittingField(chi));
    if st and p0 eq 0 then p0 := p; chip := chi; dP := dL; continue; end if;
    // second prime must ALSO be root-power strict (absolute simplicity of the
    // reduction => End^0(reduction) = Q(pi_q), so a real quadratic center would
    // embed in BOTH splitting fields; even chi (e.g. x^4-2x^2+1681) is excluded)
    if p0 ne 0 and st and Degree(SplittingField(chip*chi)) eq dP*dL then return (ord gt 18), p0, p, chip, chi; end if;
    nf +:= 1; if nf ge 14 then break; end if;
  end for;
  return false, 0, 0, P!0, P!0;
end function;

cands := [ [-3,-3/4,-3/5], [-10,-1/2,-10/7], [-5,-15/8,-15/22], [-1/2,-15/8,-15/19],
           [-4/9,4/17,-4/25], [4/17,-5/18,-10/49] ];

for c in cands do
  s := c[1]; t := c[2]; u := c[3];
  b := (Afun(t)-Afun(s))/(s^2-t^2);
  a := Afun(s) - b*(1-s^2);
  h := 1 - (7/2)*x + a*x^2 + b*x^3;
  F7 := h^2 + (x-1)^7;
  assert F7 mod x^2 eq 0;
  f := F7 div x^2;
  assert Degree(f) eq 5 and Discriminant(f) ne 0;
  for w in [s,t,u] do assert Evaluate(f, 1-w^2) eq 0; end for;
  ftype := Sort([Degree(pe[1]) : pe in Factorization(f)]);
  // integral model y^2 = d^2*f(x)  (same curve via y -> d*y)
  den := LCM([Denominator(co) : co in Coefficients(f)]);
  F := P![ co*den^2 : co in Coefficients(f) ];
  assert forall{co : co in Coefficients(F) | Denominator(co) eq 1};
  C := HyperellipticCurve(F);
  J := Jacobian(C);
  h1 := Evaluate(h,1);
  mark := J![x-1, P!(den*h1)];
  om := Order(mark);
  T := TorsionSubgroup(J);
  inv := Invariants(T);
  printf "CAND s=%o t=%o u=%o type=%o ord(marked)=%o TORSION=%o\n", s, t, u, ftype, om, inv;
  printf "  integral model y^2 = %o\n", F;
  if inv eq [2,2,14] then
    ez, p0, q0, cp, cq := EndZ(C, 56);
    printf "  EndZ=%o p0=%o q0=%o\n  chi_p0=%o\n  chi_q0=%o\n", ez, p0, q0, cp, cq;
    if ez then printf "GENERIC_2214_VERIFIED s=%o t=%o u=%o\n  f=%o\n", s,t,u,f; end if;
  end if;
end for;
print "VERIFY_DONE";
quit;
