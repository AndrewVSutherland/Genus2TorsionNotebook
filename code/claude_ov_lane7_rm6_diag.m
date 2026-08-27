// claude_ov_lane7_rm6_diag.m -- Lane 7 (overnight 2026-07-25).
//
// The 195^2 PUZZLE.  claude_ov_lane7_rm6.m established, for the sixth
// contact-7 three-root curve (s,t,u) = (-3,-3/4,-3/5):
//   * torsion exactly [2,2,14], marked order 7;
//   * conductor N = 38025 = 195^2 (PARI genus2red);
//   * the real-subfield-disc core is the CONSTANT {2} at all 164 good primes
//     p < 1000, and Q(sqrt2) embeds in the splitting field of chi_p at 148/148
//     irreducible primes  => End^0(J_Qbar) contains Q(sqrt 2);
//   * BUT S_2^new(Gamma_0(195)) has orbit degrees [1,1,1,1,3] -- NO degree-2
//     Hecke field, so there is NO GL2-type abelian surface over Q of conductor
//     195^2 with real multiplication by Q(sqrt2).
//
// If the RM were defined over Q, J would be GL2-type, hence modular (Ribet +
// Serre's conjecture), hence isogenous to A_f with cond(A_f) = level(f)^2, so
// level(f) = 195 -- contradiction.  If instead Galois acts on End^0(J_Qbar) by
// the nontrivial automorphism of Q(sqrt2), then at every prime p inert in the
// field of definition of the RM, conjugation by sqrt2 sends Frobenius pi to
// -pi, forcing a1(p) = Tr(pi) = 0 -- so a1(p) = 0 must hold at DENSITY 1/2.
//
// This script measures the three quantities that decide between those:
//   (D1) the density of good primes with a1(p) = 0;
//   (D2) whether Q(pi_p) is the SAME quartic field at all ordinary primes
//        (constant => CM by a quartic CM field; varying => RM);
//   (D3) an independent conductor (Magma Conductor + Magma's own warning
//        status), and the full local reduction data.
//   (D4) a cross-check of dim S_2^new(195) and its Hecke fields via modular
//        symbols, independent of the CuspForms/Newforms path.
//
// Markers: D1 / D2 / D3 / D4 / LANE7_RM6DIAG_DONE
SetColumns(0);
SetMemoryLimit(3*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q); Z := Integers();

if not assigned PMAX then PMAX := 1000; elif Type(PMAX) eq MonStgElt then PMAX := StringToInteger(PMAX); end if;
if not assigned FMAX then FMAX := 300; elif Type(FMAX) eq MonStgElt then FMAX := StringToInteger(FMAX); end if;

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

// ---------- D1: a1(p) = 0 density ----------
n0 := 0; ntot := 0; ndd0 := 0; nred := 0; zerop := []; a1list := [];
for p in PrimesInInterval(3, PMAX) do
  if D0 mod p ne 0 then
    chi := P!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(C, GF(p))))));
    if Degree(chi) ne 4 then continue; end if;
    ntot +:= 1;
    a1 := -Coefficient(chi,3);
    Append(~a1list, <p, a1>);
    if a1 eq 0 then n0 +:= 1; if #zerop lt 25 then Append(~zerop, p); end if; end if;
    dd := Coefficient(chi,3)^2 - 4*(Coefficient(chi,2) - 2*p);
    if dd eq 0 then ndd0 +:= 1; end if;
    if not IsIrreducible(chi) then nred +:= 1; end if;
  end if;
end for;
printf "D1 good primes < %o : %o ; a1(p)=0 at %o (%o%%) ; real-subfield disc dd=0 at %o ; chi reducible at %o\n",
  PMAX, ntot, n0, (100*n0) div ntot, ndd0, nred;
printf "D1 primes with a1=0 : %o\n", zerop;
printf "D1 first 25 (p, a1) : %o\n", a1list[1..Min(25,#a1list)];

// ---------- D2: is Q(pi_p) constant? ----------
discs := {}; pairs := [];
for p in PrimesInInterval(3, FMAX) do
  if D0 mod p eq 0 then continue; end if;
  chi := P!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(C, GF(p))))));
  if Degree(chi) ne 4 or not IsIrreducible(chi) then continue; end if;
  K := NumberField(chi);
  dK := Discriminant(MaximalOrder(K));
  Include(~discs, dK);
  if #pairs lt 20 then Append(~pairs, <p, dK>); end if;
end for;
printf "D2 distinct field discriminants of Q(pi_p) over p < %o : %o\n", FMAX, #discs;
printf "D2 discs = %o\n", Sort([d : d in discs]);
printf "D2 sample (p, disc Q(pi_p)) = %o\n", pairs;
printf "D2 VERDICT: %o\n", (#discs eq 1 select "CONSTANT => quartic CM field (CM), not RM"
                                        else "VARYING => RM (End^0 = real quadratic), not CM");

// ---------- D3: independent conductor ----------
Cm := ReducedMinimalWeierstrassModel(C);
fm, hm := HyperellipticPolynomials(Cm);
Dm := Z!Discriminant(Cm);
printf "D3 minmodel y^2 + (%o)*y = %o\n", hm, fm;
printf "D3 |disc(minmodel)| = %o = %o\n", AbsoluteValue(Dm), Factorization(AbsoluteValue(Dm));
try
  NN := Conductor(Cm);
  printf "D3 Magma Conductor = %o = %o\n", NN, Factorization(NN);
catch ee
  printf "D3 Magma Conductor failed: %o\n", ee`Object;
end try;

// ---------- D4: dim S_2^new(195) and Hecke fields, via modular symbols ----------
M := ModularSymbols(195, 2, +1);
S := CuspidalSubspace(M);
NS := NewSubspace(S);
printf "D4 dim S_2(Gamma_0(195)) (plus part) = %o, new part = %o\n", Dimension(S), Dimension(NS);
Dc := NewformDecomposition(NS);
printf "D4 new decomposition dims = %o\n", [Dimension(d) : d in Dc];
for i in [1..#Dc] do
  d := Dc[i];
  fq := qEigenform(d, 20);
  KK := BaseRing(Parent(fq));
  printf "D4 orbit %o dim %o Hecke field deg %o : %o\n", i, Dimension(d), Degree(KK),
     Degree(KK) eq 1 select "Q" else DefiningPolynomial(KK);
end for;

print "LANE7_RM6DIAG_DONE";
quit;
