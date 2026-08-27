// claude_ov_lane1_verify.m -- Lane 1 (overnight 2026-07-25)
// INDEPENDENT rebuild + exact verification of the u = -1/2 slice family of the
// contact-7 three-root surface R(s,t,u)=0.
//
//   R(s,t,-1/2)=0  <=>  2q^3+3pq^2+2pq-p^3 = 0  (p=s+t, q=st): a NODAL cubic,
//   parametrised by p = m*q:  q = 2m/((m+1)^2(m-2)),  p = 2m^2/((m+1)^2(m-2)).
//   s,t individually rational  <=>  w^2 = -m^4+6m^2+4m,  and then
//   s = (m^2+w)/((m+1)^2(m-2)),  t = (m^2-w)/((m+1)^2(m-2)).
//   The quartic is the elliptic curve E: Y^2 = X^3+6X^2-16 via X = 4/m,
//   Y = 4w/m^2.  E(Q) = Z*<(-4,4)> + Z/2*<(-2,0)>  (rank 1).
//
// For each E-point we rebuild (s,t,-1/2) -> (a,b) -> h -> f from scratch and run
//   * the three-root assertion,
//   * quintic factor type (want [1,1,1,2] => 2-rank 3),
//   * exact TorsionSubgroup(J),
//   * order of the marked class [ (1,h(1)) - infty ] (want 7),
//   * the End=Z certificate: two good primes p0,q0, each with chi irreducible
//     and NO degree drop of pi^n for n in [2..12], and with linearly disjoint
//     splitting fields (Degree(SplittingField(chi_p0*chi_q0)) = dP*dQ).
// Usage:  magma -b NLO:=-8 NHI:=-2 claude_ov_lane1_verify.m
SetColumns(0);
if not assigned NLO then NLO := -8; elif Type(NLO) eq MonStgElt then NLO := StringToInteger(NLO); end if;
if not assigned NHI then NHI := -2; elif Type(NHI) eq MonStgElt then NHI := StringToInteger(NHI); end if;
if not assigned MemGB then MemGB := 6; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals(); P<x> := PolynomialRing(Q); Z := Integers();
Afun := func< w | (2*w^5+4*w^4+6*w^3+8*w^2+10*w+5)/(2*(w+1)^2) >;

// ---- the End=Z certificate (protocol of code/claude_2214_threeroot_verify.m) ----
function EndZ(C, ord)
  D := Z!Discriminant(C); p0 := 0; chip := P!0; dP := 0; nf := 0;
  for p in PrimesInInterval(11, 400) do
    if D mod p eq 0 then continue; end if;
    chi := P!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(C, GF(p))))));
    if Degree(chi) ne 4 or not IsIrreducible(chi) then continue; end if;
    K<aa> := NumberField(chi); st := true;
    for n in [2..12] do
      if Degree(MinimalPolynomial(aa^n)) ne 4 then st := false; break; end if;
    end for;
    dL := Degree(SplittingField(chi));
    if st and p0 eq 0 then p0 := p; chip := chi; dP := dL; continue; end if;
    if p0 ne 0 and st and Degree(SplittingField(chip*chi)) eq dP*dL then
      return (ord gt 18), p0, p, chip, chi;
    end if;
    nf +:= 1; if nf ge 20 then break; end if;
  end for;
  return false, 0, 0, P!0, P!0;
end function;

// ---- cheap RM pre-screen: squarefree core of c3^2-4(c2-2p) over good primes ----
function RMsig(C)
  D := Z!Discriminant(C); sig := {};
  for p in PrimesInInterval(11,200) do
    if D mod p ne 0 then
      chi := P!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(C, GF(p))))));
      if Degree(chi) eq 4 then
        d := Coefficient(chi,3)^2 - 4*(Coefficient(chi,2) - 2*p);
        if d ne 0 then Include(~sig, Squarefree(Z!Abs(d))); end if;
      end if;
    end if;
  end for;
  return Sort([s : s in sig]);
end function;

E := EllipticCurve([0,6,0,0,-16]);
printf "E = %o\n", E;
printf "  Conductor = %o   Torsion = %o\n", Conductor(E), Invariants(TorsionSubgroup(E));
lo,hi := RankBounds(E);
printf "  RankBounds = %o %o\n", lo, hi;
Pgen := E![-4,4]; Ttor := E![-2,0];
printf "  generator P = %o, 2-torsion T = %o\n", Pgen, Ttor;

seen := {};
for n in [NLO..NHI] do
 for tw in [0,1] do
  Pt := n*Pgen; if tw eq 1 then Pt := Pt + Ttor; end if;
  if Pt eq E!0 then continue; end if;
  X := Pt[1]/Pt[3]; Y := Pt[2]/Pt[3];
  if X eq 0 then continue; end if;
  m := 4/X;
  if m eq 0 or m eq -1 or m eq 2 then continue; end if;
  w := Y*m^2/4;
  assert w^2 eq -m^4+6*m^2+4*m;                       // the slice quartic
  den := (m+1)^2*(m-2);
  s := (m^2+w)/den; t := (m^2-w)/den; u := Q!(-1/2);
  vs := [s,t,u];
  bad := false;
  for i in [1..3] do if vs[i] eq 0 or vs[i] eq -1 then bad := true; end if; end for;
  for i in [1..3] do for j in [i+1..3] do
    if vs[i] eq vs[j] or vs[i] eq -vs[j] then bad := true; end if; end for; end for;
  if bad then continue; end if;
  key := Sort([s,t]);
  if key in seen then continue; end if;  Include(~seen, key);

  // --- rebuild the chart from scratch ---
  b := (Afun(t)-Afun(s))/(s^2-t^2);
  a := Afun(s) - b*(1-s^2);
  h := 1 - (7/2)*x + a*x^2 + b*x^3;
  F7 := h^2 + (x-1)^7;
  assert F7 mod x^2 eq 0;
  f := F7 div x^2;
  if Degree(f) ne 5 or Discriminant(f) eq 0 then
    printf "n=%o tw=%o DEGENERATE\n", n, tw; continue;
  end if;
  // the three-root assertion (independent of the surface algebra)
  for v in vs do assert Evaluate(f, 1-v^2) eq 0; end for;
  ftype := Sort([Degree(pe[1]) : pe in Factorization(f)]);

  // integral model y^2 = den2^2 * f(x)
  den2 := LCM([Denominator(co) : co in Coefficients(f)]);
  F := P![ co*den2^2 : co in Coefficients(f) ];
  assert forall{co : co in Coefficients(F) | Denominator(co) eq 1};
  C := HyperellipticCurve(F);
  J := Jacobian(C);
  mark := J![x-1, P!(den2*Evaluate(h,1))];
  om := Order(mark);
  tors := Invariants(TorsionSubgroup(J));
  printf "MEMBER n=%o tw=%o m=%o\n", n, tw, m;
  printf "  (s,t,u) = (%o, %o, %o)\n", s, t, u;
  printf "  ftype=%o  ord(marked)=%o  TORSION=%o\n", ftype, om, tors;
  if tors eq [2,2,14] then
    sig := RMsig(C);
    ez, p0, q0, cp, cq := EndZ(C, 56);
    printf "  RMsig #=%o  first=%o\n", #sig, sig[1..Min(6,#sig)];
    printf "  EndZ=%o  p0=%o  q0=%o\n", ez, p0, q0;
    if ez then
      printf "  chi_p0 = %o\n  chi_q0 = %o\n", cp, cq;
      printf "FAMILY_MEMBER_CERTIFIED n=%o tw=%o  torsion=[2,2,14] EndZ p0=%o q0=%o\n", n, tw, p0, q0;
      printf "  f = %o\n", f;
      printf "  Fint = %o\n", F;
    end if;
  else
    printf "  *** torsion is NOT [2,2,14] ***\n";
  end if;
 end for;
end for;
printf "LANE1_VERIFY_DONE NLO=%o NHI=%o\n", NLO, NHI;
quit;
