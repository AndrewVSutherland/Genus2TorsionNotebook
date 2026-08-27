# Z/37 at level 2190 is not a Jacobian: the narrow-class obstruction (2026-07-30)

Companion to `notes/claude_z31_censuses_and_eisenstein_2026_07_30.md` (which flagged
2190.2.a.v, M_f = 37, as the next record candidate) and correction context for
`notes/claude_z31_rm_witness_2026_07_30.md`. Computations: session agent on
aws-spot-11, scripts/logs in `~/z31/z37/` there (`z37_step1{,b,c}`, `z_pp_driver.m`,
`z_pp_pos.m`, controls `z31_pp_control.log`, `z19_*`), using CHIMP +
ModularAbelianSurfaces clones.

## Result

**No genus-2 curve over Q has Jacobian in the isogeny class of A_f(2190.2.a.v).**
The canonical Q-members (A_f^sub and its maximal-endomorphism variant; the Hecke
order Z[sqrt3] is already maximal) admit **no rational principal polarization**: the
Pfaffian binary form of the rational self-dual hom lattice is x^2 - 3y^2 (disc 12);
unimodular classes exist but all fail positivity. Rational (1,d)-polarizations exist
exactly for d with x^2 - 3y^2 = -d solvable: verified positive for d = 2, 3, 74, 111
and empty for d = 1, 4, 5, 6, **37** (74 = 2*37 works). The root cause is
h^+(Q(sqrt3)) = 2 — the fundamental unit 2+sqrt3 has norm +1, so no totally
positive unimodular class exists.

Class-wide: quotient by a full prime kernel acts by a square (trivial on the
2-torsion narrow group); rank-1 O-kernel moves need reducible rho_bar_q, and a
resultant certificate (gcd = 2^10*3^7*37^2 over good p <= 97) confines candidates to
primes above {2,3,37}; the mixed-signature primes p_2 = (1+sqrt3) and p_3 = (sqrt3)
are proven irreducible from exact eigenvalues; the Eisenstein line at
p_37 = (7+2sqrt3) is real (a_p ≡ 1+p mod p_37 at all good p <= 97) but 7+2sqrt3 is
totally positive, so quotienting does not change the narrow class. Caveats recorded:
ramified primes 2,3 not fully closed for exotic non-O-stable lattices; the rational
37-POINT on a class member is supplied by Alessandri–Coppola Thm 4.2/4.3 (sup over
the isogeny class attains the bound), not recomputed here.

Control validating the machinery: **1830.2.a.q** (RM Q(sqrt2), h^+ = 1, unit norm
-1): exactly ONE rational pp — the known Z/31 curve. Same analysis settles
**1554.2.a.s** (M_f = 19, RM Q(sqrt6), h^+ = 2): no rational pp, no Jacobian /Q;
its two cubic-field olddata entries are two distinct cubic-moduli pp's.

```text
newform      RM field  unit norm  h^+  rational pp   curve /Q
1830.2.a.q   Q(sqrt2)  -1         1    1 (verified)  EXISTS (the Z/31 record)
2190.2.a.v   Q(sqrt3)  +1         2    0 (proven)    OBSTRUCTED
1554.2.a.s   Q(sqrt6)  +1         2    0 (proven)    OBSTRUCTED
```

## Cross-check against Alessandri–Coppola (arXiv 2602.21047v3 + their repo)

Their data repo (github.com/NirvanaC93/Torsion-of-GL2-abelian-varieties,
LMFDBandMAGMA/torsion_g2.txt, 26,440 g=2 forms) records exactly our two headline
forms — "1830.2.a.q: predicted torsion order equal to bound = 31" and "2190.2.a.v:
... = 37" — and its sorted tail equals v3's Conjecture 4.5/4.6 g=2 lists
{1..24,28,31,37,44,56} / {2,3,5,7,11,13,17,19,23,31,37}. Independent convergence
with our M_f sweep. Division of labor now established:

- A-C Thm 4.2/4.3 (proven): some GL2-type abelian variety in the class attains the
  bound — so 31 AND 37 occur for g=2 GL2-type abelian varieties.
- This lab (July 2026): 31 is attained by an explicit geometrically simple JACOBIAN
  (strictly stronger); 37 at level 2190 is NOT attainable by any Jacobian — it lives
  only on (1,2)-/(1,3)-polarized members. Their g=2 lists therefore split into
  Jacobian-realizable vs abelian-variety-only entries, and the discriminating
  invariant is the narrow class group of the Hecke field.

## The new shortlist filter, and where Z/37-on-a-Jacobian could still live

**Eisenstein Jacobian records require h^+(Hecke field) = 1** (equivalently a
norm -1 fundamental unit), else no member of the class is principally polarizable
over Q. Applying to the shortlist: 43920.2.a.cr (possible second 31; Hecke field
Q(sqrt2), h^+ = 1) stays viable; any future 37-candidate needs a dim-2 form with
37 | M_f AND h^+ = 1 Hecke field — none exists at level <= 10^4, so the hunt goes
to Yoo-admissible levels beyond LMFDB coverage (ModularSymbols lane). A (1,2)- or
(1,3)-polarized explicit model of A_f(2190.2.a.v) exhibiting the rational 37-point
would still be a nice object (not a Jacobian; outside this lab's table scope).
