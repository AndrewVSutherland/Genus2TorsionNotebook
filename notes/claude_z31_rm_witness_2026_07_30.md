# Z/31: first realization — RM witness from 1830.2.a.q (2026-07-30)

## Headline

The genus-2 curve

```text
C : y^2 + (-x^2 - x) y = -839 x^6 + 2841 x^5 - 4587 x^4 + 4300 x^3 - 2466 x^2 + 816 x - 126
```

(Edgar Costa's `ModularAbelianSurfaces` repo, `olddata/label_to_curve_QQ.txt` line 2385,
commit `687c10d`; attached to the modular form **1830.2.a.q**) has

```text
J(Q)_tors = Z/31Z        (exact, Magma TorsionSubgroup)
```

**This is the first known genus-2 Jacobian over Q with a rational point of order 31** —
the smallest previously-unrealized cyclic order (dossier: Elkies 2001-02 had only an
order-31 *subgroup without rational generator*, `notes/claude_torsion_refs_dossier.md`
trap 4). It is also a new record *prime* order for geometrically simple abelian surfaces
over Q (previous: 29, Lep95; the alpha DB's largest prime orders are 23 and 29). J is
geometrically simple but NOT generic: RM by Q(sqrt 2). Found by GPT 5.6 Sol computing
torsion of the curves in Costa's `olddata` directory; independently verified here.

It confirms the entry 31 of the *conjectural* GL2-type torsion-order set
{1..24, 28, 31, 37, 44, 56} recorded from gl2tors (Alessandri–Coppola) in the dossier
(trap 7 there — the list is conjectural, for GL2-type surfaces, not Jacobians per se).

## Verification (2026-07-30, aws-spot-11, this session)

Scripts: session scratchpad `spot_verify31.m` (to be committed as
`code/claude_z31_verify.m`); log `~/z31/verify31.log` on aws-spot-11 (copy below).

- Simplified integral model: `y^2 = F`,
  `F = -3356 x^6 + 11364 x^5 - 18347 x^4 + 17202 x^3 - 9863 x^2 + 3264 x - 504`
  (= 4f + h^2; content-free). `disc(F) = -2^23 * 3^15 * 5^11 * 61^2`. Bad primes
  {2,3,5,61} = support of the level 1830 = 2*3*5*61. F is irreducible over Q
  (2-rank 0, consistent with odd torsion).
- `TorsionSubgroup(J) = [31]` (0.3 s); generator (Mumford)
  `u = x^2 - 21/20 x + 9/20,  v = 101/200 x + 51/200`, order 31.
  **The class is invisible**: disc(u) = -279/400 < 0, so the supporting points are
  complex-conjugate; the curve has no small rational points (none with x-height <= 20;
  leading coeff -3356 < 0 not a square, so no rational points at infinity either).
- **Geometric simplicity certificate** (D4/power-transform, the lab standard):
  witness prime p = 11, `chi = x^4 - 8x^3 + 36x^2 - 88x + 121` irreducible, no degree
  drop for pi^n, n in [2..12]. CERTIFIED geometrically simple.
- **RM by Q(sqrt2)** via the *repaired* RM screen (sound direction "constant core => RM",
  `notes/claude_top10_2026_07_26.md` §1): n_p = a1^2 - 4(a2 - 2p) has squarefree core 2
  at ALL 164 good primes p <= 1000 tested (PARI, local), zero exceptions. Independently:
  LMFDB `mf_newforms` 1830.2.a.q has Hecke field 2.2.8.1 = Q(sqrt2), dim 2, no CM, and
  its Hecke traces match the curve's Frobenius traces at every checked prime
  (a_7=0, a_11=8, a_13=-2, a_17=6, ...). So End(J_Qbar) ⊗ Q = Q(sqrt2) (GL2-type over Q,
  is_gl2_type would be true).
- 31 | #J(F_p) at all 164 good primes p <= 1000 (necessary condition, PARI).

## Why nothing had caught it

- Not in `g2c_curves` (66k) nor `g2c_curves_new` (6.2M), and no torsion order divisible
  by 31 anywhere in the alpha DB. **Structural reason** (critic correction 2026-07-30):
  the curve is GL2-type, so its Jacobian's conductor is 1830^2 = 3,348,900 > 2^20 =
  1,048,576 — the alpha DB is the cond < 2^20 database, so this curve could never have
  been in it. (Same for the in-DB RM witnesses at 138^2 = 19044 and 390^2 = 152100,
  which fit only because their levels are small.) The 2026-07-25 histogram
  (`claude_ov_db_histogram_2026_07_25.md`) was complete *for the DB*; its
  `is_simple_geom` filter would additionally have hidden a *split* [31] had one existed
  (cf. 169.a1, torsion [19], geom_end_alg = M_2(Q)).
- **Priority nuance (record honestly):** the torsion ORDER 31 for the isogeny class
  A_f of 1830.2.a.q has been sitting in Costa's public `olddata/modav-torsion.out`
  (`1830.2.a.q:31`, repo data unchanged since 2024-08-14). Nobody surfaced it. What is
  new in 2026: the exact Jacobian group J(Q)_tors = Z/31Z on the explicit curve, the
  invisible-generator observation, the geometric-simplicity certificate, and the
  first-realization framing against the literature (Platonov 2014 survey: primes known
  only through 29; caution: modav-torsion.out records A_f-class orders, which both
  under- and over-state model torsion — 390.2.a.h:28 vs its DB Jacobians' [2,2,14],
  834.2.a.j:23 vs its Q-curve's trivial torsion).
- **Alessandri–Coppola status (CORRECTED 2026-07-30, after Ari's catch):** arXiv
  2602.21047 has THREE versions. v1/v2 (Feb–Mar 2026) capped g=2 GL2-type torsion
  primes at {2,3,5,7,11,13,19} and orders at {1..16,18..22,24,28,44,56} — this curve
  refutes THOSE (historical footnote only). **v3 (28 Apr 2026), the current version,
  renumbers and expands: Conjecture 4.5 (orders, g=2) = {1..24, 28, 31, 37, 44, 56};
  Conjecture 4.6 (primes, g=2) = {2,3,5,7,11,13,17,19,23,31,37}** — v3 cites Costa's
  dataset [CEH+], so the additions {17,23,31,37} trace to the same isogeny-class
  orders our sweep found. This curve therefore CONFIRMS v3's entry 31, and is
  strictly stronger than the conjecture's setting (Jacobian, not just GL2-type
  abelian variety). Dossier trap 7's transcription {1..24,28,31,37,44,56} matches
  v3's g=2 row exactly (an earlier note here calling it "their g=3 data" was wrong;
  retracted). On v3's entry 37: see notes/claude_z37_obstruction_2026_07_30.md —
  at level 2190 the 37 CANNOT live on a Jacobian /Q (narrow-class obstruction), so
  the entry, if correct, is realized only by non-principally-polarized members.
- The 31-class is invisible AND the curve is pointless: every lab construction
  (contact, D_inf/CF, named charts) forces *visible* classes, so no lab chart could
  have produced it.
- Eisenstein numerology: 61 || 1830 and 61 + 1 = 62 = 2*31 — consistent with a
  cuspidal/Eisenstein origin of the 31-torsion at the level prime 61 (being chased via
  the literature lane; would give a *predictive* criterion for other levels/primes).

## Status of the generic question

Open: find a curve with J(Q)_tors ⊇ Z/31 and End(J_Qbar) = Z. Analogous to the [2,22]
situation (RM-only witnesses; generic open, paper row ^RM). Session program: three fiber
sieves (chart U' universal-Mumford/Kummer-ladder; chart I CF-31; Elkies-Kumar Y_-(8)
RM-sqrt2 exhaustion) + broad 31-divisibility box scan; see
`notes/claude_z31_generic_program_2026_07_30.md` (forthcoming) and code
`code/claude_z31_*.{c,m}`.

## Verbatim verification log (aws-spot-11)

```text
curve y^2 + (-x^2 - x)y = -839*x^6 + 2841*x^5 - 4587*x^4 + 4300*x^3 - 2466*x^2 + 816*x - 126
simplified integral model: y^2 = -3356*x^6 + 11364*x^5 - 18347*x^4 + 17202*x^3 - 9863*x^2 + 3264*x - 504
disc = [ <2, 23>, <3, 15>, <5, 11>, <61, 2> ]
TORSION_INVARIANTS: [ 31 ]   (0.3 s)
TORSION_GEN 1: (x^2 - 21/20*x + 9/20, 101/200*x + 51/200, 2)  order 31
SIMPLICITY: true  witness prime 11  chi = x^4 - 8*x^3 + 36*x^2 - 88*x + 121
BadPrimes(C): [ 2, 3, 5, 61 ]
```
