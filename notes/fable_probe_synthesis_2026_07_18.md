# Fable session: complementary probes and ranked top-10 (2026-07-18)

Context: run in parallel with (and after a sync of)
the collaborator session that produced
`geometrically_simple_torsion_frontier_2026_07_18.md`.  Everything here is
ADDITIVE: no shared file was edited, and no git commit was made (the `.git`
directory is on a synced working copy and the collaborator session was live-writing;
committing from two machines into a synced `.git` risks corruption — commit
from one side only).

## Independent record check

`scratchpad` re-verification of the collaborator's curve (also see
`code/verify_record_22212_order96.m`): exact torsion `[2,2,2,12]`, order 96,
in 0.31s; `chi_37 = 1369x^4 - 148x^3 + 30x^2 - 4x + 1` irreducible.  Agrees
with the frontier note's two-prime root-power certificates at 37 and 73.

## Four new probes (all small; ≤ 2 CPU-min each; results in results/fable_*)

### P1. Richelot probe of the two simple Z/24 seeds
`code/fable_z24_richelot_probe.m`, `results/fable_z24_richelot_probe.log`

Curve A (r=5, p=-5/2, t=-9/2): **zero** rational Richelot kernels.
Curve B (r=1/3, p=-1/9, t=-1): exactly one; its codomain has exact torsion
`[2,12]` (order 24 — the 2-primary part redistributes DOWN, as in the record
component's `[2,12]` leaves).  Conclusion: the fixed Z/24 seeds are
Richelot-poor; `[2,24]`/`[4,12]` cannot be reached by graph moves from the
known seeds — deform at family level (matches the frontier doctrine).

### P2a. A(12)+3 chart-level locality for [3,12]
`code/fable_a12_a8_localprobe.m`, `results/fable_a12_a8_localprobe.log`

Chart-level containment of `[3,12]` in `J(F_q)` over the full good A(12)
chart: **0 at q=5 (4 good pts), 0 at q=7 (24 good pts), 32/316 at q=11**
(canary `[12]`-containment = 100% of good points, so the transcription is
right).  The M(2,12)+3 wall at 5 is therefore NOT chart-specific: on both
order-12 charts, a global `[3,12]` member must be forced into bad
reduction/boundary at 5 (and, on A(12), also at 7).  Strategy consequence:
boundary-first (CRT lifts with prescribed bad reduction at {5,7}), not
affine scans; and try to PROVE the wall (local-obstructions methodology) —
either outcome is progress.

### P2b. A(8) chart-level locality for [8,8]
Same files.  Containment of `[8,8]`: 6/98 (q=7), 39/589 (q=11), 106/1097
(q=13), with local groups as big as `[16,16]` and `[2,2,8,8]`.  The
`[8,8]`-compatible locus is robustly alive at every probed prime —
consistent with the frontier's good-open points on the reduced
second-halving cover.  `[8,8]` has the healthiest local picture of all open
targets.

### P3. M_1(8,2,2) chart-level locality for [2,2,24] and [2,2,16]
`code/fable_m1822_localprobe.m`, `results/fable_m1822_localprobe.log`

The chart is entirely degenerate mod 5 (good=0).  `[2,2,24]`-containment is
**0 at q=7, 11, 13** and alive only at 17 — the m3222_plus3 (+3) shortcut to
`[2,2,24]` is walled at three consecutive primes; prefer the HLP-anchor
deformation route.  `[2,2,16]`-containment is 0 at q=7 but alive at 11, 13,
17 — so p=7 is a chart-fundamental wall while the known 11-adic wall must be
specific to the halving cover; target bad-reduction-at-7 strata.

### P4. Raw contact-6 chart-level locality for [2,6,6]
`code/fable_contact6_266_localprobe.m`,
`results/fable_contact6_266_localprobe.log`

`[2,6,6]`-containment over the raw (a,b) chart: **0 at q=5 and q=7**, alive
at 11 (9 hits) and 13 (21 hits); canaries `[2,6]`, `[6,6]` healthy.  The
frontier note's cover-level emptiness mod 5/7 is thus already true at the
crude chart level: any global `[2,6,6]` here needs bad reduction at both 5
and 7.  Same "small-prime confinement" shape as `[3,12]`.

## The emergent doctrine

2-primary-heavy targets (`[8,8]`, `[2,2,16]` away from 7) are locally roomy;
targets mixing a 3-part with extra 2-structure (`[3,12]`, `[2,6,6]`,
`[2,2,24]`-on-M_1(8,2,2)) are confined to small-prime boundary strata in
every probed chart.  For those, the productive move is CRT-first: prescribe
the bad-reduction behavior at the walled primes, lift, and only then search —
or prove the wall and close the chart.

## Ranked top-10 open targets (this session's synthesis)

1. `[2,24]` (48) — no detected local obstruction anywhere; three lanes
   (A(8)+3 family-level, A(12) descent funnel, Richelot-vary over
   A(2,2,2,12)); P1 shows fixed-seed graphs are exhausted, so deform.
2. `[8,8]` (64) — healthiest local picture (P2b); do the elimination /
   genus-one pair-fiber ranks on the reduced second-halving cover.
3. `[3,12]` (36) — smallest open order; boundary-confined at {5,7} on both
   order-12 charts (P2a); boundary-first CRT or prove the wall.
4. `[2,2,24]` (96) — ten split HLP anchors; deformation with transverse
   2-parameter planes; +3 shortcut walled (P3).
5. `[6,12]` (72) — Prym rank 1 < 2: complete the certification chain, then
   precision-safe Coleman + Abel–Prym MW sieve; realize or close.
6. `[2,6,6]` (72) — boundary-confined at {5,7} (P4); CRT/boundary program on
   the contact-6 [1,2,2] locus with the simple [6,6] core as anchor.
7. `[2,2,16]` (64) — wall is 7-fundamental, 11-cover-specific (P3);
   bad-reduction-at-7 strata on the halving cover.
8. `[4,12]` (48) — mod-7 boundary of the s=m^2 equation on full M(2,12);
   plus [4,12]-conditions on Richelot codomains over A(2,2,2,12) (P1 shows
   the 2-redistribution mechanism is real).
9. `[2,2,4,8]` (128) — arithmetically feasible (HLP), simple [2,4,8] seed
   exists; pair-fiber analysis and Richelot deformation of the seed.
10. `[2,2,2,24]` (192) — the next-record halving cover over the proven
    A(2,2,2,12) locus; global cover open (actively worked by the
    collaborator session today).

First alternates: `[2,2,4,12]` (192, kernel-first Richelot), `Z/48`, `[60]`,
`Z/35`, `Z/5 x Z/5`, `[4,16]`.  Deliberately excluded: `[2,2,2,14]` — no
machinery meshes a 7-class with full rational 2-torsion, and both sessions
independently ended up without a route.
