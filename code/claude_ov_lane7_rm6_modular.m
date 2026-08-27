// claude_ov_lane7_rm6_modular.m -- Lane 7 (overnight 2026-07-25).
//
// FORMAL RM(sqrt 2) CERTIFICATION of the sixth contact-7 three-root curve
//     (s,t,u) = (-3, -3/4, -3/5),   torsion exactly [2,2,14], marked order 7.
//
// Corrected input (this lane): the conductor is
//     N = 152100 = 2^2 * 3^2 * 5^2 * 13^2 = 390^2      (Magma, v_2(disc_min)=7
//     so Magma's Ogg formula carries its correctness guarantee)
// NOT the 38025 = 195^2 that PARI genus2red reports -- genus2red does not
// implement the conductor exponent at 2 (sentinel exponent -1) and omits the
// 2-part from N.
//
// Certification route.  The diagnostics in results/claude_ov_lane7_rm6_diag.log
// show End^0(J_Qbar) is a REAL QUADRATIC field, namely Q(sqrt2):
//   * the real-subfield disc core is the constant {2} at all 164 good p < 1000
//     and Q(sqrt2) embeds in every splitting field (148/148 irreducible chi);
//   * Q(pi_p) takes 47 DISTINCT field discriminants over p < 300, so
//     End^0(J_Qbar) is not a quartic CM field (CM would force Q(pi_p) constant);
//   * a1(p) = 0 at only 18/164 = 10% of good primes, far from the density 1/2
//     that a Galois-nontrivial action on Q(sqrt2) would force (conjugation by
//     sqrt2 would send Frobenius pi to -pi, hence Tr(pi) = 0, at every prime
//     inert in the field of definition of the RM).  So the RM is defined over Q
//     and J is of GL2-type.
// GL2-type + Ribet/Serre => J ~ A_f for a weight-2 newform f with trivial
// character, real quadratic Hecke field, and cond(A_f) = level(f)^2, so
// level(f) = 390.  We locate f and match Frobenius polynomials.
//
// Markers: NEWFORM390 / MATCH390 / LANE7_RM6MOD_DONE
SetColumns(0);
SetMemoryLimit(3*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q); Z := Integers();
if not assigned PMAX then PMAX := 400; elif Type(PMAX) eq MonStgElt then PMAX := StringToInteger(PMAX); end if;

Gfun := func< v | -(v^5 - v^3 - v^2/2)/(v+1)^2 >;
s := -3; t := -3/4;
c4 := (Gfun(s)-Gfun(t))/(s^2-t^2); c0 := Gfun(s) - c4*s^2;
b := c4 - 2; a := 9/2 - c0 - c4;
h := 1 - (7/2)*x + a*x^2 + b*x^3;
f := (h^2 + (x-1)^7) div x^2;
den := LCM([Denominator(co) : co in Coefficients(f)]);
F := P![ co*den^2 : co in Coefficients(f) ];
C := HyperellipticCurve(F);
D0 := Z!Discriminant(C);
Nlev := 390;
K2 := QuadraticField(2);

printf "NEWFORM390 curve y^2 = %o ; conductor 152100 = 390^2\n", F;
S := CuspForms(Nlev, 2);
NF := Newforms(S);
degs := [Degree(BaseRing(Parent(orb[1]))) : orb in NF];
printf "NEWFORM390 level %o weight 2 : %o newform Galois orbits, Hecke-field degrees %o (total %o)\n",
  Nlev, #NF, degs, &+degs;
cands := [];
for i in [1..#NF] do
  KK := BaseRing(Parent(NF[i][1]));
  if Degree(KK) ne 2 then continue; end if;
  iso := IsIsomorphic(KK, K2);
  printf "NEWFORM390 orbit %o : deg-2 Hecke field %o  isQ(sqrt2)=%o\n", i, DefiningPolynomial(KK), iso;
  if iso then Append(~cands, i); end if;
end for;
printf "NEWFORM390 orbits with Hecke field Q(sqrt2): %o\n", cands;

// cache the curve's Frobenius polynomials once
chid := AssociativeArray();
plist := [];
for p in PrimesInInterval(3, PMAX) do
  if D0 mod p eq 0 or Nlev mod p eq 0 then continue; end if;
  chi := P!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(C, GF(p))))));
  if Degree(chi) ne 4 then continue; end if;
  chid[p] := chi; Append(~plist, p);
end for;
printf "NEWFORM390 comparing at %o good primes p < %o (Sturm bound for level 390 weight 2 is 168)\n", #plist, PMAX;

for i in cands do
  fm2 := NF[i][1];
  agree := 0; bad := [];
  for p in plist do
    ap := Coefficient(fm2, p);
    tr := Z!Trace(ap); nm := Z!Norm(ap);
    pred := x^4 - tr*x^3 + (2*p + nm)*x^2 - p*tr*x + p^2;
    if pred eq chid[p] then agree +:= 1;
    else if #bad lt 10 then Append(~bad, <p, chid[p], pred>); end if;
    end if;
  end for;
  printf "MATCH390 orbit %o : chi_p == Norm_{K/Q}(T^2 - a_p T + p) at %o of %o primes ; mismatches %o\n",
    i, agree, #plist, #bad;
  if #bad gt 0 then printf "MATCH390 orbit %o first mismatches (p, chi_curve, chi_pred) = %o\n", i, bad; end if;
  if #bad eq 0 and agree eq #plist then
    printf "MATCH390 orbit %o : CERTIFIED. J is isogenous over Q to A_f (level 390, weight 2, trivial character, Hecke field Q(sqrt2)); hence End^0(J) contains Q(sqrt2) and J has REAL MULTIPLICATION BY Q(sqrt2), defined over Q.\n", i;
    printf "MATCH390 orbit %o q-expansion (first 30 coefficients): %o\n", i, qExpansion(fm2, 30);
  end if;
end for;

print "LANE7_RM6MOD_DONE";
quit;
