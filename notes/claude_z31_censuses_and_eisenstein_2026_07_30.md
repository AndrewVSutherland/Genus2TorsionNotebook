# Censuses, the Eisenstein engine, and the RM torsion frontier (2026-07-30)

Companion to `notes/claude_z31_rm_witness_2026_07_30.md`. All exact computations on
aws-spot-11 (192 cores), Magma; artifacts in `data/claude_z31_*`.

## ERRATUM (supersedes commit 89cc8d5's message)

Commit 89cc8d5 called 1290.2.a.t a "third [2,22] RM witness". **Wrong.** Its
G2-invariants and bad primes {2,3,5,43} coincide exactly with the lab's known second
witness "BLP C4corr" ((p1..p4) = (-34/15,32/45,31/15,-2/9) in the b2p chart): same
moduli point, different model. The correct statement: there are still exactly TWO known
[2,22] moduli points, and both now have modular identities and Eisenstein explanations:

```text
19044.h.2   <-> 138.2.a.d   level 138  = 2*3*23,  RM sqrt5, 11 = num((23-1)/12)  (Mazur at 23)
BLP C4corr  <-> 1290.2.a.t  level 1290 = 2*3*5*43, RM sqrt5, 11 | 43+1           (Yoo at 43)
```

Both were confirmed [2,22], D4-simple, RM-core {5} in this session
(spot-11 `validate_tail.log`, `blp_g2.log`).

## Census 1: Costa olddata, full invariant factors (15,216 curves)

`data/claude_z31_olddata_torsion_invariants.txt` — exact `TorsionSubgroup` for every
curve in `olddata/label_to_curve_QQ.txt` (150-worker sweep, ~4 min, zero errors).
Independently confirms the collaborator's scan: **[31] (1830.2.a.q) is the only
torsion group in the file not already realized**; [2,22] appears exactly for
138.2.a.d and 1290.2.a.t (the two known moduli points, above).

Large-order tail, now ALL D4-certified geometrically simple with RM cores
(`validate_tail.log`):

```text
138.2.a.d [2,22] core 5 | 1290.2.a.t [2,22] core 5 | 1830.2.a.q [31] core 2
39.2.a.b [2,14] core 2  | 390.2.a.h [2,14] core 2  | 4110.2.a.w [23] core 2
322.2.a.f [22] core 3   | 74.2.a.b [19] core 5     | 453.2.a.a [19] core 5
555.2.a.c [19] core 5   | 618.2.a.j [17] core 2    | 670.2.a.f [17] core 2
462.2.a.h [2,2,6] core 3
```

Every one fits the Eisenstein pattern ell | q -+ 1 at a prime q of the level
(37+1=2*19 twice, 151+1=8*19, 103-1=6*17, 67+1=4*17, 137+1=6*23, 23-1=2*11,
43+1=4*11, 13+1=2*7, 61+1=2*31). 13-for-13.

## Census 2: genus2isogenies new models (618,490 curves) — NEGATIVE

van Bommel–Chidambaram–Costa–Kieffer isogeny-completion curves (typical classes,
End = Z): exact torsion for all 618,490 new models
(`data/claude_z31_isog_census_summary.txt`). 42 distinct groups, max order 40
([2,2,10]); **every group already realized generically in g2c_curves_new**. No new
groups from isogeny-twisting the cond < 2^20 world.

## Census 3: Bruin–Stoll small-coefficient curves (196,171 classes) — NEGATIVE

All squarefree deg-5/6 f with coefficients in [-3,3] (no conductor bound — the one
public dataset covering huge-conductor escapees at tiny height): exact torsion for all
(`data/claude_z31_bruinstoll_census_summary.txt`). 21 groups, max order 24. Nothing
new: tiny-height sextics do not carry exotic torsion.

## The Eisenstein prediction engine (Q2's real answer)

`data/claude_z31_eisenstein_shortlist.csv` + `code/compute_mf.py`,
`code/pass_b_enhanced.py`. For every dim-2 weight-2 trivial-character newform in the
LMFDB (all 80,387: 19,129 exact via a_p vectors to p<=199; 61,258 via a validated
trace-chain method to p<=47), compute M_f = gcd_p Norm(1+p-a_p) — a rigorous multiple
of the torsion order of any member of the A_f isogeny class. Results:

- Retro-finds 1830.2.a.q (M_f = 31) deterministically — the discovery required no luck.
- Predictive hit rate 7/8 on shortlisted forms with known Q-curves (the miss,
  834.2.a.j, has M_f = 23 but trivial model torsion: M_f is an isogeny invariant,
  torsion is not — its class member carrying the 23 is an open follow-up).
- Quadratic-character forms (5,269): empty above 13. Large Eisenstein torsion is a
  trivial-character (J_0 cuspidal) phenomenon.
- **Headline: 2190.2.a.v** — level 2190 = 2*3*5*73, Hecke field Q(sqrt3), M_f = 37
  (37 | 73+1; the exact analog of 61+1 = 2*31). Was the Z/37 record candidate;
  **RESOLVED NEGATIVELY 2026-07-30**: no member of its isogeny class is a Jacobian
  /Q — narrow-class obstruction h^+(Q(sqrt3)) = 2, no rational principal
  polarization at all (same for 1554.2.a.s, M_f = 19, Q(sqrt6)). See
  `notes/claude_z37_obstruction_2026_07_30.md`. New shortlist filter: Eisenstein
  Jacobian records need h^+(Hecke field) = 1.
- Trace-chain candidates needing eigenvalue recomputation: 43920.2.a.cr
  (M_odd = 31 — the prime 61 again; a possible second Z/31 level),
  13680.2.a.ch (29), 14246.2.a.f (19).
- Theory: Yoo, Trans. AMS 371 (2019) Thm 1.3(3) — level 1830 with AL signature
  (+,+,+,-) at (2,3,5,61), 61 ≡ -1 (mod 31), is a proven-admissible mod-31 Eisenstein
  configuration (sign conventions still to be re-checked against the published text).
  Mazur-route levels for records were scanned: no dim-2 newform exists at 149 (37),
  83 (41), 311/373 (31) — within LMFDB coverage the shortlist above is complete.
  Beyond level 10^4 is unswept (ModularSymbols at Yoo-admissible levels = future lane).

## Where this leaves Q2

Within all swept data, torsion groups realized by geometrically simple non-generic
Jacobians and by NO generic one: **[2,22], [2,2,14], [31]** (alpha DB census + this
session; [2,2,14]'s witnesses 152100.eb2/eb3 are the 390.2.a.h class). All three are
RM with Eisenstein origin. The engine says where more will appear (and 2190.2.a.v
says the next one may be Z/37).

## Addendum: the 834.2.a.j miss, probed (2026-07-30)

The predicted-23 form's Q-curve has trivial torsion AND its unique Richelot
neighbour (integral rescale of the (2,2)-isogenous Jacobian) also has trivial
torsion (`probe834b.m` on spot-11). So the 23-carrier, if any Jacobian carries
it at all, sits behind an odd-degree isogeny or on a non-Jacobian member of the
class. Calibration lesson for the shortlist: M_f is an isogeny-class invariant;
attainment by the PARTICULAR curve in Costa's file is the 7/8-frequent but not
guaranteed case, and Richelot neighbours need not rescue a miss.
