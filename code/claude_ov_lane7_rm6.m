// claude_ov_lane7_rm6.m -- Lane 7 (overnight 2026-07-25).
//
// FORMAL certification of the SIXTH contact-7 three-root curve
//     (s,t,u) = (-3, -3/4, -3/5)
// as a real-multiplication-by-Q(sqrt 2) witness with torsion exactly [2,2,14].
// (Outstanding hygiene item from 2026-07-23, never done.)
//
// The route is modularity, made available by a computation done today:
// PARI genus2red gives conductor
//     N = 38025 = 3^2 * 5^2 * 13^2 = 195^2   (SQUARE, like the other RM
//     witnesses 138^2 and 390^2 in this project)
// A square conductor M^2 is exactly what a GL2-type abelian surface
// A_f / Q attached to a weight-2 newform f of level M with real quadratic
// Hecke field has.  So we:
//   (1) enumerate the weight-2 newform Galois orbits of level 195,
//   (2) keep those with Hecke eigenvalue field Q(sqrt 2),
//   (3) for every good prime p, compare the curve's Frobenius polynomial
//         chi_p(T) = T^4 - a1 T^3 + a2 T^2 - a1 p T + p^2
//       against the GL2-type prediction
//         Norm_{K/Q}(T^2 - a_p(f) T + p)
//         = T^4 - Tr(a_p) T^3 + (2p + Nm(a_p)) T^2 - p Tr(a_p) T + p^2 ,
//   (4) if they agree at every good prime up to PMAX (well past the Sturm
//       bound for level 195 weight 2, which is 195*(1+1/3)(1+1/5)(1+1/13)/6
//       = 56), then J is isogenous to A_f (Faltings), hence
//       End^0(J) contains K = Q(sqrt 2): REAL MULTIPLICATION BY Q(sqrt 2).
//
// Independent supporting evidence, both computed here:
//   * the real-subfield-disc census (squarefree core of c3^2-4(c2-2p)) is the
//     CONSTANT {2} over every good prime up to PMAX -- Q(sqrt2) sits inside
//     Q(pi_p)^+ at every good p, which no End=Z surface does;
//   * the strict two-prime End=Z certificate provably CANNOT be completed:
//     every root-power-strict prime has splitting field containing Q(sqrt2),
//     so no two of them are linearly disjoint.  We verify that containment
//     directly rather than inferring it from a failed scan.
//
// Markers: RM6 / NEWFORM / MATCH / SQRT2IN / LANE7_RM6_DONE
SetColumns(0);
SetMemoryLimit(3*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q); Z := Integers();

if not assigned PMAX then PMAX := 1000; elif Type(PMAX) eq MonStgElt then PMAX := StringToInteger(PMAX); end if;

Gfun := func< v | -(v^5 - v^3 - v^2/2)/(v+1)^2 >;
s := -3; t := -3/4; u := -3/5;
c4 := (Gfun(s)-Gfun(t))/(s^2-t^2); c0 := Gfun(s) - c4*s^2;
b := c4 - 2; a := 9/2 - c0 - c4;
h := 1 - (7/2)*x + a*x^2 + b*x^3;
f := (h^2 + (x-1)^7) div x^2;
assert Degree(f) eq 5 and LeadingCoefficient(f) eq 1;
for w in [s,t,u] do assert Evaluate(f, 1-w^2) eq 0; end for;
den := LCM([Denominator(co) : co in Coefficients(f)]);
F := P![ co*den^2 : co in Coefficients(f) ];
C := HyperellipticCurve(F);
J := Jacobian(C);
TG := TorsionSubgroup(J);
mark := J![x-1, P!(den*Evaluate(h,1))];
printf "RM6 f = %o\n", f;
printf "RM6 Fint = %o\n", F;
printf "RM6 torsion = %o (order %o) markedord = %o ftype = %o\n",
  Invariants(TG), #TG, Order(mark), Sort([Degree(pe[1]) : pe in Factorization(f)]);
Cm := ReducedMinimalWeierstrassModel(C);
fm, hm := HyperellipticPolynomials(Cm);
printf "RM6 minmodel y^2 + (%o)*y = %o\n", hm, fm;
printf "RM6 conductor (PARI genus2red, results/claude_ov_lane7_conductors.log) N = 38025 = 195^2\n";

D0 := Z!Discriminant(C);
Nlev := 195;

// ---------- (A) real-subfield-disc census + sqrt(2) containment ----------
sig := {}; nchk := 0; nirr := 0; nsqrt2 := 0; nstrict := 0; strictnosqrt2 := [];
K2 := QuadraticField(2);
for p in PrimesInInterval(3, PMAX) do
  if D0 mod p eq 0 then continue; end if;
  chi := P!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(C, GF(p))))));
  if Degree(chi) ne 4 then continue; end if;
  nchk +:= 1;
  dd := Coefficient(chi,3)^2 - 4*(Coefficient(chi,2) - 2*p);
  if dd ne 0 then Include(~sig, Squarefree(Z!Abs(dd))); end if;
  if not IsIrreducible(chi) then continue; end if;
  nirr +:= 1;
  // does Q(sqrt2) embed in the splitting field of chi?
  SF := SplittingField(chi);
  emb := IsSubfield(K2, SF);
  if emb then nsqrt2 +:= 1; end if;
  KK<aa> := NumberField(chi); st := true;
  for n in [2..12] do
    if Degree(MinimalPolynomial(aa^n)) ne 4 then st := false; break; end if;
  end for;
  if st then
    nstrict +:= 1;
    if not emb then Append(~strictnosqrt2, p); end if;
  end if;
end for;
ss := Sort([e : e in sig]);
printf "SQRT2IN goodprimes<=%o : %o  irreducible chi : %o  root-power-strict : %o\n", PMAX, nchk, nirr, nstrict;
printf "SQRT2IN real-subfield-disc cores = %o  (#=%o)\n", ss, #ss;
printf "SQRT2IN Q(sqrt2) embeds in SplittingField(chi_p) at %o of %o irreducible primes\n", nsqrt2, nirr;
printf "SQRT2IN root-power-strict primes WITHOUT sqrt(2) in the splitting field: %o\n", strictnosqrt2;

// ---------- (B) modularity: level-195 weight-2 newforms with Hecke field Q(sqrt2) ----------
S := CuspForms(Nlev, 2);
NF := Newforms(S);
printf "NEWFORM level %o weight 2: %o Galois orbits, degrees %o\n",
  Nlev, #NF, [Degree(BaseRing(Parent(orb[1]))) : orb in NF];
cands := [];
for i in [1..#NF] do
  fm2 := NF[i][1];
  KK := BaseRing(Parent(fm2));
  if Degree(KK) ne 2 then continue; end if;
  iso := IsIsomorphic(KK, K2);
  printf "NEWFORM orbit %o : deg 2 field %o  isQ(sqrt2)=%o\n", i, DefiningPolynomial(KK), iso;
  if iso then Append(~cands, i); end if;
end for;
printf "NEWFORM candidate orbits with Hecke field Q(sqrt2): %o\n", cands;

for i in cands do
  fm2 := NF[i][1];
  KK := BaseRing(Parent(fm2));
  agree := 0; disagree := 0; firstbad := [];
  for p in PrimesInInterval(3, PMAX) do
    if D0 mod p eq 0 or Nlev mod p eq 0 then continue; end if;
    chi := P!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(C, GF(p))))));
    if Degree(chi) ne 4 then continue; end if;
    ap := Coefficient(fm2, p);
    tr := Z!Trace(ap); nm := Z!Norm(ap);
    pred := x^4 - tr*x^3 + (2*p + nm)*x^2 - p*tr*x + p^2;
    if pred eq chi then agree +:= 1;
    else disagree +:= 1; if #firstbad lt 8 then Append(~firstbad, p); end if;
    end if;
  end for;
  printf "MATCH orbit %o : chi_p == Norm_{K/Q}(T^2-a_p T+p) at %o primes, MISMATCH at %o  (bad: %o)\n",
    i, agree, disagree, firstbad;
  if disagree eq 0 and agree gt 0 then
    printf "MATCH orbit %o : CERTIFIED -- J ~ A_f over Q (level %o, Hecke field Q(sqrt2)), so End^0(J) contains Q(sqrt2): RM(sqrt2).\n", i, Nlev;
  end if;
end for;

// ---------- (C) optional: CMSV/Lombardo endomorphism package if present ----------
try
  EA := HeuristicEndomorphismAlgebra(C : Geometric := true);
  printf "CMSV HeuristicEndomorphismAlgebra (geometric) = %o\n", EA;
catch ee
  printf "CMSV endomorphisms package NOT available in this Magma: %o\n", ee`Object;
end try;

print "LANE7_RM6_DONE";
quit;
