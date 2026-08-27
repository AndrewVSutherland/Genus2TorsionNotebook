// claude_ov_lane1_certify.m -- Lane 1, overnight 2026-07-25.
//
// STRICT End(J_Qbar) = Z CERTIFICATES for members of the u = -1/2 contact-7
// [2,2,14] family (code/claude_ov_lane1_family.m).
//
// For each member the script does, on the INTEGRAL model:
//   (a) exact TorsionSubgroup(J)  and the marked-class order (must be 7);
//   (b) the quintic factor type and the rational 2-torsion rank;
//   (c) a SCAN over good primes for the strict Frobenius certificate:
//         chi irreducible over Q  AND  deg MinPoly(pi^n) = 4 for n = 2..12;
//       records EVERY certifying prime, not a fixed one;
//   (d) the two-prime RM exclusion: a pair (p0,q0) of certifying primes whose
//       Frobenius SPLITTING FIELDS are linearly disjoint over Q, i.e.
//         [ SplittingField(chi_p0 * chi_q0) : Q ]
//              = [SplittingField(chi_p0):Q] * [SplittingField(chi_q0):Q].
//       (An even chi such as x^4-2x^2+1681 can never supply this.)
//   (e) the cheap RM pre-screen: the squarefree core of c3^2-4(c2-2p) over
//       ~30 good primes.  CONSTANT => RM;  SCATTERS => End = Z.
//
// Usage: magma -b KLIST:="2,3,4,5,6,7" code/claude_ov_lane1_certify.m

SetColumns(0);
if not assigned MemGB then MemGB := 16; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
if not assigned PMAX then PMAX := 400; elif Type(PMAX) eq MonStgElt then PMAX := StringToInteger(PMAX); end if;
if not assigned KLIST then KLIST := "2,3,4,5,6,7"; end if;
Q := Rationals(); Z := Integers();
P<x> := PolynomialRing(Q);
Gfun := func<u0 | -(u0^5 - u0^3 - u0^2/2)/(u0+1)^2>;

E := EllipticCurve([0,6,0,0,-16]);
Gp := E![4,-12]; Tp := E![-2,0];

klist := [StringToInteger(s) : s in Split(KLIST, ",")];

for k in klist do
 for eps in [0,1] do
  Pt := k*Gp + eps*Tp;
  if Pt eq E!0 then continue; end if;
  xc := Pt[1]/Pt[3];
  if xc eq 0 then continue; end if;
  mm := 4/xc; ww := (Pt[2]/Pt[3])*mm^2/4;
  d := (mm+1)^2*(mm-2);
  if d eq 0 then continue; end if;
  s0 := (mm^2+ww)/d; t0 := (mm^2-ww)/d; u0 := Q!(-1/2);
  vs := [s0,t0,u0];
  if #{v : v in vs} ne 3 or #{v^2 : v in vs} ne 3 then continue; end if;
  if (-1 in vs) or (0 in vs) then continue; end if;
  c4 := (Gfun(s0)-Gfun(t0))/(s0^2-t0^2); c0 := Gfun(s0) - c4*s0^2;
  b := c4 - 2; a := 9/2 - c0 - c4;
  h := 1 - (7/2)*x + a*x^2 + b*x^3;
  num := h^2 + (x-1)^7;
  if num mod x^2 ne 0 then continue; end if;
  f := num div x^2;
  if Degree(f) ne 5 or Discriminant(f) eq 0 then continue; end if;
  mden := LCM([Denominator(cc) : cc in Coefficients(f)]);
  F := P![ Z!(Coefficient(f,i)*mden^(6-i)) : i in [0..5] ];
  C := HyperellipticCurve(F); J := Jacobian(C);
  tors := Invariants(TorsionSubgroup(J));
  ft := Sort([Degree(tt[1]) : tt in Factorization(F)]);
  tr := #Invariants(TwoTorsionSubgroup(J));
  Ppt := C ! [mden*1, mden^3*Evaluate(h,1)];
  D := J ! (Ppt - PointsAtInfinity(C)[1]);
  ordD := Order(D);
  D0 := Z!Discriminant(C);

  printf "\n==== MEMBER k=%o eps=%o  m=%o ====\n", k, eps, mm;
  printf "  (s,t,u) = (%o, %o, -1/2)\n", s0, t0;
  printf "  (a,b) = (%o, %o)\n", a, b;
  printf "  f = %o\n", f;
  printf "  integral model F = %o\n", F;
  printf "  ftype=%o  2rank=%o  markedclassorder=%o  TORSION=%o\n", ft, tr, ordD, tors;

  // --- prime scan: strict certificate ---
  certs := [];      // <p, chi>
  sig := {};        // squarefree cores of the real-subfield discriminant
  sigseq := [];
  for p in PrimesInInterval(11, PMAX) do
    if D0 mod p eq 0 then continue; end if;
    if (Z!LeadingCoefficient(F)) mod p eq 0 then continue; end if;
    Fp := PolynomialRing(GF(p))!F;
    if Degree(Fp) lt 5 or not IsSquarefree(Fp) then continue; end if;
    chi := P!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(C,GF(p))))));
    // RM pre-screen quantity
    dd := Coefficient(chi,3)^2 - 4*(Coefficient(chi,2) - 2*p);
    if dd ne 0 then
      sc := Squarefree(Z!Abs(dd));
      Include(~sig, sc); Append(~sigseq, <p, sc>);
    end if;
    if not IsIrreducible(chi) then continue; end if;
    K := NumberField(chi); pi := K.1; drop := false;
    for nn in [2..12] do
      if Degree(MinimalPolynomial(pi^nn)) lt 4 then drop := true; break; end if;
    end for;
    if not drop then Append(~certs, <p, chi>); end if;
  end for;
  printf "  certifying primes (chi irreducible + no power drop n=2..12): %o\n", [c[1] : c in certs];
  ss := Sort([z : z in sig]);
  printf "  RM pre-screen: #distinct squarefree cores = %o over %o good primes; cores = %o\n",
     #ss, #sigseq, ss[1..Min(10,#ss)];
  if #ss eq 1 then printf "  ** WARNING: constant real-subfield signature %o -> RM suspected **\n", ss; end if;

  // --- two-prime linear disjointness ---
  if #certs lt 2 then
    printf "  NO two-prime certificate (only %o certifying primes up to %o)\n", #certs, PMAX;
  else
    done := false;
    for i in [1..Min(6,#certs)] do
      for j in [i+1..Min(8,#certs)] do
        p0 := certs[i][1]; chi0 := certs[i][2];
        q0 := certs[j][1]; chi1 := certs[j][2];
        L0 := SplittingField(chi0); L1 := SplittingField(chi1);
        d0 := Degree(L0); d1 := Degree(L1);
        Lc := SplittingField(chi0*chi1);
        if Degree(Lc) eq d0*d1 then
          printf "  *** STRICT CERTIFICATE ***\n";
          printf "    p0=%o chi=%o  [splitting field deg %o]\n", p0, chi0, d0;
          printf "    q0=%o chi=%o  [splitting field deg %o]\n", q0, chi1, d1;
          printf "    compositum degree %o = %o * %o  => LINEARLY DISJOINT => End(J_Qbar)=Z\n", Degree(Lc), d0, d1;
          done := true; break;
        end if;
      end for;
      if done then break; end if;
    end for;
    if not done then printf "  no linearly-disjoint pair found among the first certifying primes\n"; end if;
  end if;
 end for;
end for;
printf "\nLANE1_CERTIFY_DONE\n";
quit;
