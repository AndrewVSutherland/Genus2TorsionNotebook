# top10 item 5: (2,6,6), order 72 — contact-6 extra-root chart scan extended to height 150

Date: 2026-07-17.  Test-run agent report.  All computation in the scratchpad
(`/tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/a266v/`);
no repo file modified.  Total compute ~7 CPU-min, single-threaded, niced.

## Dossier summary (inputs)

Target [2,6,6] = C2 x C6 x C6, order 72, ranked #5.  Split realizations
abound (HLP positive-rank elliptic surface).  The contact-6 extra-root chart

    h6 = 1 + a x + b x^2 + x^3,   f = h6^2 - (x-1)^6,
    h6(r) = eps (r-1)^3  =>  a = (eps (r-1)^3 - 1 - r^3 - b r^2)/r,

already produces exact [2,6,6] curves, but tier2 (notes/claude_tier2_266_312_contact6.md)
certified that at parameter height <= 60 they are all ONE nonsimple curve.
The simple [6,6] win lives on the same contact-6 chart (core slice a=133/39,
b=-7/13).  p=5 and p=7 kill-tables admit zero good-reduction residues.  The
144 | #J(F_p) filter is empirically validated but unproven.  Spec'd next
step: extend the scan to height ~150 (~10 CPU-min estimate).

## Structural facts established/used (small but load-bearing)

1. f = q1 * (x*q2) with q1 = h6 - (x-1)^3 = (b+3)x^2 + (a-3)x + 2 (a
   QUADRATIC) and x*q2 = h6 + (x-1)^3, q2 = 2x^2 + (b-3)x + (a+3).
   The extra-root condition h6(r) = eps(r-1)^3 is exactly q1(r)=0 (eps=+1)
   or q2(r)=0 (eps=-1).  Hence on this chart "extra rational 2-torsion via a
   split quadratic factor" and "extra rational root" are the SAME condition;
   factor type [1,1,1,2] and 2-rank 3 are automatic, generic chart torsion
   is already [2,2,6], and the only arithmetic accident needed for [2,6,6]
   is one extra rational 3-torsion class.  (The "escape via quadratic factor
   instead of a root" idea in the assignment is therefore vacuous in-chart.)
2. The chart involution (a,b) <-> (b,a) (x -> 1/x, y -> y/x^3) preserves the
   curve and swaps the eps=+1 / eps=-1 branches (verified on the data below).
3. a = -3 <=> x^2 | f: a globally singular rational curve inside the chart
   (for eps=-1 it is the line b = 3-2r).  It passes every residue kill-table
   (bad reduction everywhere) and must be excluded exactly; it accounted for
   ~2000 false sieve candidates before the exact kill was added.

## Validation (calibration, done first)

Script `validate.m` / `validate2.m` (Magma):
- Simple [6,6] curve (a=133/39, b=-7/13): chart reproduces the stated model
  y^2 = 11389248x^5 - ... (equal G2-invariants); TorsionSubgroup = [6,6].
  CORRECTION to the dossier: at p=23, chi = x^4-26x^2+529 is irreducible but
  its 12th-power transform is NOT (alpha^12 has degree 2), so p=23 certifies
  only under the old Lp-irreducibility criterion.  Under the repo's STRICT
  certificate (chi irreducible AND 12th-power transform irreducible of
  degree 4) the first certifying prime is p = 37,
  L_37 = 1369x^4 - 370x^3 + 54x^2 - 10x + 1.  The curve stays certified.
- Nonsimple [2,6,6] curve at (eps,r,b) = (+1, 21, -23/7): a = 187/21
  reproduced; G2-invariants match the stated model; f = x*q1*q2 identity
  holds; factor type [1,1,1,2]; torsion exactly [2,6,6];
  L_11 = (11T^2+1)^2, L_13 = (13T^2-2T+1)^2 (splits, nonsimple); no strict
  certificate for p <= 97; 144 | #J(F_p) at every good p <= 100 (clean model
  -63x^5+1500x^4-4026x^3+6364x^2+2625x visibly has bad reduction at 5, 7).

## Pipeline (rebuilt from scratch; tier2 scripts were not in the repo)

1. `tables.gp` (PARI, 25 s): sound necessary-condition kill-tables p=5..67
   from hyperellcharpoly (2-rank>=3 via disc(q1)/disc(q2) QR; (T-1)^2 | chi
   mod 3; 9 | chi(1)).  Reproduced "p=5 zero allowed good residues" but kept
   12 at p=7 — too weak by a factor ~300 in survivors.
2. `tables_exact.m` (Magma, ~100 s): EXACT-structure tables p=5..67, allowed
   iff J(F_p) has >=3 invariant factors even and >=2 divisible by 3.
   Reproduces the dossier: p=5 AND p=7 admit zero allowed good residues;
   p=31 gives exactly 420 structure-allowed (tier2's number).
3. `sieve266.c` (C, single-threaded): residues of all rationals of height
   <= H per prime (infinity residues allowed), 17 bit-table lookups per pair
   sorted by kill power, then exact int128 Q-level kills: b=-3, r in {0,1},
   h6(1)=0, a=-3, disc(q1)=0, disc(q2)=0.  (The disc-square test is retained
   but provably vacuous — see fact 1.)
   CALIBRATION at H=60: exactly the 12 tier2 candidates (the 8 known-curve
   parametrizations + the 4 killed-only-by-144 rows), nothing else, 0.16 s.
4. `deep.m` (Magma, 47 s): per candidate — exact structure containment at
   good p in (67,149]; 144 | #J(F_p) at all good p <= 149 (log-only);
   strict simplicity certificate hunt p in [11,97]; exact TorsionSubgroup;
   G2-invariant census.
5. `verdict.m` / `classes.m`: GeometricAutomorphismGroup, L_p factorization
   patterns, extended certificate hunt p <= 250, class identification.

Exact commands (from the scratch dir):

    nice -n 15 magma -b validate.m            # then validate2.m
    nice -n 15 gp -q tables.gp
    nice -n 15 magma -b tables_exact.m
    gcc -O2 -march=native -o sieve266 sieve266.c -lm
    nice -n 15 ./sieve266 60  tables_exact_p67.txt sieve_h60x.txt   # calibration
    nice -n 15 ./sieve266 150 tables_exact_p67.txt sieve_h150.txt   # 6.9 s
    nice -n 15 magma -b candfile:=sieve_h150.txt deep.m
    nice -n 15 magma -b verdict.m             # then classes.m

## Main result: census at height <= 150 (both eps, full chart)

1,504,754,940 parameter pairs; 52,126 table survivors; 452 sieve candidates;
17 stage-A survivors — ALL with exact torsion [2,6,6], ALL passing
144 | #J(F_p) at every good p <= 149, NONE with a strict simplicity
certificate (p <= 250 for class representatives).  They form exactly
**4 distinct Q-isomorphism classes**, every one **bielliptic
(GeometricAutomorphismGroup of order 4), hence geometrically split**:

| class | (a,b) representatives (up to (a,b)<->(b,a)) | params | integral model y^2 = |
|---|---|---|---|
| K (known, tier2) | (187/21,-23/7), (-1/42,-13/7), (19/2,-67/25) | 10 | -63x^5+1500x^4-4026x^3+6364x^2+2625x |
| N1 (new) | (173/9,-141/50), (-457/819,-129/91) | 4 | 72900x^5+6357861x^4-17498700x^3+70642900x^2+9000000x |
| N2 (new) | (-1639/2184,-113/91) at r=16/15 | 1 | 530712x^5-2365515x^4+3954474x^3-2938571x^2+819000x |
| N3 (new) | (-435/169,115/18), eps=-1, r=-81/52, -16/117 | 2 | 173765124x^5+191275045x^4-100769292x^3+40744980x^2+7884864x |

(N1's two (a,b) orbits are Q-isomorphic through a nontrivial Moebius map,
not the standard involution; the third model listed under fI in verdict.log,
16773120x^5-..., is the r=9/8 parametrization of N1.)

**Main negative (certified): no geometrically simple [2,6,6] exists on the
contact-6 extra-root chart with parameters (eps, r, b) of height <= 150.**

Side facts:
- The dossier's "one curve up to height 60" is now "four curves up to
  height 150": the [2,6,6] locus is NOT a single curve class, it keeps
  acquiring rational points — but every one lands on the bielliptic locus.
- The unproven 144-filter gains support: 4 distinct curves x ~30 good primes
  each, zero violations.  Notably the filter made ZERO kills in this
  pipeline: exact structure at p <= 149 alone cut 452 -> 17 (the tier2
  [2,2,6]-type near-misses die at structure primes in (67,149]).  A deep
  enough sound filter subsumes the unproven one here.

## Structural reason for the failure, sharpened

The [2,6,6] locus on this chart = (extra-root, automatic 2-part) AND (extra
3-torsion, a degree-40-type cubic-contact cover condition).  Its rational
points at height <= 150 lie without exception on the bielliptic
(Aut = C2xC2) locus, where J ~ E1 x E2 and [2,6,6] points are cheap because
the HLP split machinery produces them.  This is the same phenomenon tier2
proved for [3,12] on M(2,12) (rational points of the irreducible carrier
S12 confined to the split boundary class): an arithmetic thin-set
confinement on an irreducible cover, not a component phenomenon.  Extending
the height is therefore expected to keep producing split curves; the chart
is structurally dry for SIMPLE [2,6,6].

## Strategy (ranked routes)

1. **RM-first construction (recommended escape).**  All realized groups in
   this project came from constructions where the whole group is forced by
   algebraic identities.  For C2xC6xC6 = (Z/2)^3 x (Z/3)^2 the forcing
   decomposes primewise: take a family with rational (3,3)-structure
   (3-rank 2) — e.g. the Bruin–Flynn–Shnidman sqrt3-RM full level-3
   rational moduli (arXiv:2102.04319, already flagged for (3,12)) — and
   impose 2-rank 3 on it.  On a quintic model f = x*g1*g2 the two splitting
   conditions disc(gi) = square are CONICS in the remaining parameters
   (rationally parametrizable fiberwise, cf. fact 1), so the 2-conditions
   cost no arithmetic accidents; gluing 2- and 3-parts forces
   C2xC6xC6 <= J(Q)_tors identically.  RM families are generically
   geometrically simple, so the bielliptic attractor of the contact-6 chart
   is avoided by construction.  First step: write the BFS level-3 chart in
   quintic form and check which factor types its universal curve admits.
2. **Bielliptic-locus quotient of the present chart.**  Compute the
   bielliptic curve B in the (r,b)-chart (Shioda/Igusa condition for
   Aut >= C2xC2 restricted to the family) through the 17 points.  Payoffs:
   (a) parametrize B to mass-produce split [2,6,6] (and settle whether the
   [2,6,6] locus equals a thin set ON B — upgrading the height-150 negative
   to a precise conjecture, as tier2 did for [3,12] via S12); (b) the
   complement chart-minus-B with the 3-torsion cover gives the exact object
   whose rational points would be simple [2,6,6] — test local solvability
   there before any further sweeping.
3. **Height extension H=300-400 of the present sieve** (cheap: the C sieve
   scales at ~5 ns/pair; H=300 is ~2 min, deep phase dominated by candidate
   count).  Low expected yield given the attractor pattern — run only as a
   byproduct/regression while route 1 is being set up, or to harvest more
   split classes for route 2's curve-fitting.

## Next steps

- Set up the BFS sqrt3-RM level-3 chart in odd-degree form; determine the
  conic conditions for 2-rank 3; check whether the combined family is a
  rational surface (route 1, first session).
- Derive the bielliptic-locus equation of the contact-6 chart and verify it
  interpolates all 9 (a,b) points above (route 2; also yields the precise
  "confinement" conjecture statement).
- Prove the 144-lemma for this family (2-rank-3 + marked order-6 structure
  should force v_2 #J(F_p) >= 4, v_3 >= 2); now 4 curves of evidence.
- Bank the three NEW split [2,6,6] curves (N1, N2, N3) in the split census
  data if the project keeps one.

## Files (scratchpad)

Scripts: validate.m, validate2.m, tables.gp, tables_exact.m, sieve266.c,
deep.m, verdict.m, classes.m.
Data: tables_allowed.txt (necessary-condition tables, kept for reference),
tables_exact_p67.txt (exact tables used), sieve_h60x.txt (12-candidate
calibration), sieve_h150.txt (452 candidates), deep_h150.log (census),
verdict.log, validate.log, validate2.log, tables_exact.log.
