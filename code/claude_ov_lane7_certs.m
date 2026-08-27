// claude_ov_lane7_certs.m -- Lane 7 (overnight 2026-07-25).
//
// INDEPENDENT rebuild + strict End=Z certification of ALL ELEVEN contact-7
// three-root points (6 known + 5 new from the 2026-07-25 harvest), in a fresh
// Magma session, per the validate-and-record-a-hit protocol.
//
// The rebuild uses the *G-recipe* stated in the lane brief
//     G(v) = -(v^5 - v^3 - v^2/2)/(v+1)^2
//     c4 = (G(v1)-G(v2))/(v1^2-v2^2),  c0 = G(v1) - c4*v1^2
//     b  = c4 - 2,                     a  = 9/2 - c0 - c4
// and CROSS-CHECKS it against the independent *A-recipe* used by the 07-23
// verification script (code/claude_2214_threeroot_verify.m)
//     A(w) = (2w^5+4w^4+6w^3+8w^2+10w+5)/(2(w+1)^2)
//     b = (A(t)-A(s))/(s^2-t^2),  a = A(s) - b*(1-s^2)
// The two derivations are algebraically distinct; agreement is the independent
// rebuild.
//
// For each curve:
//   * assert f = (h^2+(x-1)^7)/x^2 monic quintic, three asserted rational roots
//   * quintic factor type, exact TorsionSubgroup, order of the marked 7-class
//   * RM pre-screen: squarefree core of the real-subfield disc c3^2-4(c2-2p)
//     over the good primes in [11,400]  (CONSTANT => RM, SCATTERS => End=Z)
//   * STRICT End=Z certificate: SCAN for p0 (chi irreducible, all root powers
//     pi^n of degree 4 for n=2..12) and then for q0 with the SAME strictness
//     whose splitting field is linearly disjoint from p0's
//     (deg SF(chi_p0 * chi_q0) = deg SF(chi_p0) * deg SF(chi_q0)).
//   * reduced minimal Weierstrass model, discriminant factorisation, conductor
//
// Output markers: TRIPLE / CERT / MINMODEL / COND / LANE7_CERTS_DONE
SetColumns(0);
SetMemoryLimit(3*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q); Z := Integers();

Gfun := func< v | -(v^5 - v^3 - v^2/2)/(v+1)^2 >;
Afun := func< w | (2*w^5+4*w^4+6*w^3+8*w^2+10*w+5)/(2*(w+1)^2) >;

trips := [
  [* "known", [-10, -10/7, -1/2] *],
  [* "known", [-5, -15/8, -15/22] *],
  [* "known-RM", [-3, -3/4, -3/5] *],
  [* "known", [-15/8, -15/19, -1/2] *],
  [* "known", [-5/18, -10/49, 4/17] *],
  [* "known", [-4/9, -4/25, 4/17] *],
  [* "NEW", [-511/61, -511/625, -1/2] *],
  [* "NEW", [-165/41, -33/16, -165/289] *],
  [* "NEW", [-164/297, -1/2, 164/361] *],
  [* "NEW", [-17/50, -34/189, 34/121] *],
  [* "NEW", [-1/2, -13/49, 13/50] *]
];

PMAX := 400;

for item in trips do
  tag := item[1]; T := item[2];
  s := T[1]; t := T[2]; u := T[3];

  // --- G-recipe (lane brief) ---
  c4 := (Gfun(s)-Gfun(t))/(s^2-t^2);
  c0 := Gfun(s) - c4*s^2;
  bG := c4 - 2;  aG := 9/2 - c0 - c4;
  // --- A-recipe (independent, 07-23 script) ---
  bA := (Afun(t)-Afun(s))/(s^2-t^2);
  aA := Afun(s) - bA*(1-s^2);
  agree := (aG eq aA) and (bG eq bA);

  a := aG; b := bG;
  h := 1 - (7/2)*x + a*x^2 + b*x^3;
  num := h^2 + (x-1)^7;
  assert num mod x^2 eq 0;
  f := num div x^2;
  assert Degree(f) eq 5 and LeadingCoefficient(f) eq 1 and Discriminant(f) ne 0;
  for w in [s,t,u] do assert Evaluate(f, 1-w^2) eq 0; end for;

  ftype := Sort([Degree(pe[1]) : pe in Factorization(f)]);

  // integral model y^2 = den^2 * f(x)   (y -> den*y, same curve)
  den := LCM([Denominator(co) : co in Coefficients(f)]);
  F := P![ co*den^2 : co in Coefficients(f) ];
  assert forall{co : co in Coefficients(F) | Denominator(co) eq 1};
  C := HyperellipticCurve(F);
  J := Jacobian(C);
  mark := J![x-1, P!(den*Evaluate(h,1))];
  om := Order(mark);
  TG := TorsionSubgroup(J);
  inv := Invariants(TG); tord := #TG;

  printf "TRIPLE %o (s,t,u)=(%o, %o, %o) recipes_agree=%o ftype=%o torsion=%o markedord=%o\n",
         tag, s, t, u, agree, ftype, inv, om;
  printf "  f = %o\n", f;
  printf "  Fint = %o   (y^2 = Fint, den=%o)\n", F, den;

  // ---- prime data over the good primes ----
  D0 := Z!Discriminant(C);
  goodp := []; chis := []; strictL := []; sfdeg := []; rsdisc := [];
  for p in PrimesInInterval(11, PMAX) do
    if D0 mod p eq 0 then continue; end if;
    chi := P!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(C, GF(p))))));
    if Degree(chi) ne 4 then continue; end if;
    Append(~goodp, p); Append(~chis, chi);
    dd := Coefficient(chi,3)^2 - 4*(Coefficient(chi,2) - 2*p);
    Append(~rsdisc, dd eq 0 select 0 else Squarefree(Z!Abs(dd)));
    if not IsIrreducible(chi) then
      Append(~strictL, false); Append(~sfdeg, 0); continue;
    end if;
    K<aa> := NumberField(chi); st := true;
    for n in [2..12] do
      if Degree(MinimalPolynomial(aa^n)) ne 4 then st := false; break; end if;
    end for;
    Append(~strictL, st);
    Append(~sfdeg, st select Degree(SplittingField(chi)) else 0);
  end for;

  sig := Sort([ss : ss in Seqset([r : r in rsdisc | r ne 0])]);
  nstrict := #[i : i in [1..#goodp] | strictL[i]];
  printf "  RMprescreen: #goodprimes=%o  #distinct real-subfield-disc cores=%o  cores=%o\n",
         #goodp, #sig, sig;
  printf "  strictprimes: %o of %o  -> %o\n", nstrict, #goodp,
         [goodp[i] : i in [1..#goodp] | strictL[i]];

  // ---- STRICT certificate: scan p0 then q0 ----
  p0 := 0; q0 := 0; chip := P!0; chiq := P!0; tried := 0;
  idxs := [i : i in [1..#goodp] | strictL[i]];
  for ii in idxs do
    p0 := goodp[ii]; chip := chis[ii]; dP := sfdeg[ii];
    for jj in idxs do
      if jj eq ii then continue; end if;
      tried +:= 1;
      if tried gt 60 then break; end if;
      dQ := sfdeg[jj];
      if Degree(SplittingField(chip*chis[jj])) eq dP*dQ then
        q0 := goodp[jj]; chiq := chis[jj]; break;
      end if;
    end for;
    if q0 ne 0 then break; end if;
    if tried gt 60 then break; end if;
  end for;

  qm := (tord gt 18);
  cert := (p0 ne 0) and (q0 ne 0) and qm;
  printf "CERT (%o, %o, %o) p0=%o q0=%o disjoint_pairs_tried=%o tors_gt_18=%o ==> EndZ=%o\n",
         s, t, u, p0, q0, tried, qm, cert;
  printf "  chi_p0 = %o\n", chip;
  printf "  chi_q0 = %o\n", chiq;

  // ---- minimal model, discriminant, conductor ----
  Cm := ReducedMinimalWeierstrassModel(C);
  fm, hm := HyperellipticPolynomials(Cm);
  Dm := Z!Discriminant(Cm);
  printf "MINMODEL (%o, %o, %o) y^2 + (%o)*y = %o\n", s, t, u, hm, fm;
  printf "  |disc| = %o = %o\n", AbsoluteValue(Dm), Factorization(AbsoluteValue(Dm));
  // conductor is computed OUTSIDE Magma with PARI genus2red (Liu's algorithm);
  // Magma's Conductor() falls back on Ogg's formula when v_2(disc) >= 12 and
  // explicitly disclaims correctness there.  Emit gp-readable coefficients.
  printf "GPMODEL (%o, %o, %o) | %o | %o\n", s, t, u, fm, hm;
  // re-verify torsion on the minimal model (independent of the chart model)
  Cv := HyperellipticCurve(4*fm + hm^2);
  printf "  minmodel torsion = %o\n", Invariants(TorsionSubgroup(Jacobian(Cv)));
end for;

print "LANE7_CERTS_DONE";
quit;
