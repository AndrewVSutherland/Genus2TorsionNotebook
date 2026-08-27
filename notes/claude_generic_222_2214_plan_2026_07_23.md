# Generic `[2,22]` and `[2,2,14]`: strategy portfolio and day-1 results (2026-07-23)

Goal: geometrically simple witnesses with `End(Jac_Qbar) = Z` for the two
torsion groups whose only known geometrically simple witnesses are RM:

- `[2,22]` (order 44): RM witness 19044.h.2 (conductor `138^2`),
  `y^2+(x^2+x)y = x^6-3x^5+9x^4-5x^3+12x^2-6x`;
- `[2,2,14]` (order 56): RM witness 152100.eb.2 (conductor `390^2`),
  `y^2+(x^2+x)y = 9x^6+69x^5+123x^4-95x^3-183x^2+165x-35`.

Square conductors = GL(2)-type: the cheap torsion factory (Eisenstein/modular)
*forces* RM, which is why the DB shows RM-only for these groups.  Generic
witnesses require geometric constructions.  Production LMFDB has NO `[2,22]`
or `[2,2,14]` at all; generic `[22]` = {1192.a, 1312.c}, generic `[2,14]` =
{1416.b} only.

Prior art in repo: commit 34cfb59 (2026-07-22): Richelot BFS from DB seeds
(blocked — 2-rank-1 seeds have almost no rational kernels) and blind
2-rank-gated height scans (hopeless: odd torsion is codim 2).  Losing
pattern per project memory; winning pattern = build the odd torsion into a
family algebraically, then parametrize the remaining 2-condition.  That is
what this plan does.

## Day-1 headline

**`[2,2,14]`: SIX candidate curves found** by a new "three-root" route on the
contact-7 chart (details below).  All six pass `56 | #J(F_p)` at every good
prime tested (`p <= 59`), none forces 112/168 (so no systematic extra
torsion), and each carries a root-power geometric-simplicity certificate at
`p = 17` (computed via PARI `hyperellcharpoly`).  Exact `TorsionSubgroup` +
`End = Z` verification is queued in `code/claude_2214_threeroot_verify.m` —
**blocked only by a transient Magma licence issue** on the compute box
(restored by re-running the local licence setup script).

`[2,22]`: no candidate yet; several sub-routes *measured dead* today; ranked
live portfolio below.

## New structural fact: every witness has its odd class at infinity

CF order of `D_inf = inf+ - inf-` (pell-cf-order machinery, gp port,
self-tests f14→14, f18→18 passed):

| curve | torsion | factor type | CFOrd(D_inf) |
|---|---|---|---|
| 1192.a (generic) | [22] | [1,1,4] | **11** |
| 1312.c (generic) | [22] | [2,4] | **11** |
| 1416.b (generic) | [2,14] | [1,1,1,3] | **7** |
| 19044.h.2 (RM) | [2,22] | [1,1,2,2] (2-rank 2) | **11** |
| 152100.eb.2 (RM) | [2,2,14] | [1,1,1,1,2] (2-rank 3) | **7** |

So both targets live on "(2-rank shape) ∧ (CF-order n at infinity)" surfaces
that provably have rational points — the RM witnesses are points on them
(sitting on the Humbert divisor).  Mod-p census (p = 101, 211; N = 2·10^5
per shape; `results/claude_cf_order_census.log`): on the `[2,2,2]`-shape,
orders 11 and 22 occur with density `c/p^2`, `c ≈ 6–7.5`; on the
`[1,1,1,1,2]`-shape, orders 7 and 14 with `c ≈ 11–17`.  Nonempty, no local
obstruction, healthy constants; mod-p seeds recorded for interpolation
routes.  (Flynn and Daowsud–Schmidt are CFOrd = 11 at every sampled t —
they are curves ON the order-11 surface.)

## `[2,2,14]`: the three-root route (primary; candidates in hand)

Chart (contact-7, `notes/contact7_family.md`): `h = 1-(7/2)x+ax^2+bx^3`,
`f = (h^2+(x-1)^7)/x^2` (monic quintic), marked class `[P - inf]` of order 7,
`P = (1, h(1))`.

Root mechanism: every rational root `r` of `f` satisfies `1-r = v^2` and
`h(r) = v^7`, i.e.

```text
a + b(1-v^2) = A(v),   A(v) = (2v^5+4v^4+6v^3+8v^2+10v+5)/(2(v+1)^2)
```

— LINEAR in `(a,b)` (this is the two-root chart of
`code/contact7_two_root_surface.m`, a rational surface).  Hence:

- two roots `(s,t)`: solve linearly — the known `[2,14]`-surface;
- **three roots `(s,t,u)`: one compatibility equation** — collinearity of the
  three points `(1-v_i^2, A(v_i))` on a rational plane quintic.  Explicitly
  `R(s,t,u) = 0` where `R` = residual cubic after removing `(u-s)(u-t)` from
  the numerator of `A(u) - a(s,t) - b(s,t)(1-u^2)`.

Type `[1,1,1,2]` quintic ⇒ 2-rank 3 ⇒ torsion ⊇ `(Z/2)^3 × Z/7 = [2,2,14]`.
Unlike the "+3" covering condition studied on this chart in 2026-07-11
(degree-40 cover, genus ≥ 15, correctly stopped), the third-root condition is
a **degree-3 cover** of a rational surface.

Structure of `R` (`code/claude_2214_threeroot_surface.gp`): S3-symmetric,
tridegree (3,3,3); coefficient of `u^3` is `2(s+1)^2(t+1)^2(s+t)`; in
elementary symmetric coordinates `e1,e2,e3` it is **quadratic in `e3`**:

```text
R = 2(e2-2)e3^2 + (4e2^2+4e1e2-2e1-4e2-1)e3 + e2(2(e1+e2)^2+e1)
  = 2e2(e1+e2+e3)^2 + (s+t)(t+u)(u+s) - 2e3(e1+2e2+2e3).
```

Sweep (gp, `code/claude_2214_threeroot_sweep.gp`; H = 16 then 32): 16 hits =
**6 distinct curves**, all `[1,1,1,2]`, all gates passed, all with root-power
certificates.  Orbits `(s,t,u)`:

```text
(-3, -3/4, -3/5)        f = x^5 + 4769/400 x^4 + 9009/400 x^3 - 104671/1600 x^2 + 1699/40 x - 42/5
(-10, -1/2, -10/7)
(-5, -15/8, -15/22)
(-1/2, -15/8, -15/19)
(-4/9, 4/17, -4/25)
(4/17, -5/18, -10/49)
```

Every orbit has the form `v_i = c/k_i` (common numerator; e.g. the last is
`20/(85, -72, -98)`) — low-height rational curves on the surface are likely;
find them (Magma geometry pass) and we get an infinite family, not sporadic
hits.

Next steps (in order):
1. **Magma exact verification** (`code/claude_2214_threeroot_verify.m`,
   minutes): exact `TorsionSubgroup`, marked-class order, `End = Z`
   certificate (root-power + disjoint splitting fields + LSSV via order
   56 > 18).  Then the full validate-and-record protocol (fresh-session
   rebuild, integral model, commit + certificates).
2. Exactness fallbacks are cheap: with 6 candidates (and more from a deeper
   sweep) an occasional `[2,4,14]`/`[2,2,28]`/extra-odd member does not hurt.
3. Fit the rational curve(s) through the hit orbits (the `M(24)`-component
   playbook) → family theorem rather than single curves.
4. Bonus watches on the same surface: `SPLITALL` members (5 rational roots)
   ⇒ `[2,2,2,14]` (order 112); halving of a Weierstrass class ⇒ `[2,2,28]`;
   independent 3-torsion ⇒ order-168 territory.  None needed for the table,
   all record-relevant.
5. Independent 7-source for diversity: LPS2004 has *two 2-parameter
   7-families* (one containing X_0(29)) — literature agent is retrieving the
   explicit models; same three-root game applies to them.

## `[2,22]`: ranked portfolio (no candidate yet)

Known 11-sources: Flynn `F_t` and Daowsud–Schmidt `G_u` (both quadratic in
the parameter, both order-11-at-infinity families; DS2018 = J. Number Theory
189 (2018) 115-130, corrected family in arXiv:1708.05511v2), the branch-point
`[22]`-subfamilies (`notes/order22_from_order11.md`), BLP2009's 17 sporadic
curves (agent retrieving), Leprévost 1995 `C_22`, and the two DB curves.

Measured DEAD today (`code/claude_222_qf_incidence.gp`,
`code/claude_222_flynn62_enum.gp`, `code/claude_222_comp43_sweep.gp`):

- The quadratic-factor incidence `{x^2+ux+v | F_t}` splits (t-resultant) into
  a (4,3) plane quartic component and a (6,2) component for BOTH families.
  The (6,2) double covers are elliptic: Flynn `y^2 = x^3+3x^2+2x+1`, DS
  `y^2 = x^3-3x^2+2x+1`, both conductor 92, both analytic rank 1, so both
  families have **infinitely many `[2,4]`-type members** (a second
  `[22]`-component, new).  *Caution recorded: PARI `factor()` drops the
  constant — the discriminants are `-16(u∓1)^2(cubic)`, and the naive
  (untwisted) curve gave rank 0; always re-attach the constant.*
- BUT the `[2,22]`-decisive next step — residual quartic splits into two
  rational quadratics (⇒ `[2,2,2]` or `[1,1,2,2]` type, 2-rank 2) — is
  **empty empirically**: 0/22 Flynn-(6,2) members from elliptic points (heights
  grow fast), 0/46 + 0/46 on the (4,3) components of Flynn/DS (u-height ≤ 48).
  Same fate as the branch-point upgrades (rank-0 elliptic + pointless
  genus-3, note above).  The easy in-family covers are all cold.

Ranked live routes:

- **B1 (primary): CF-backward construction on the shape, RM-anchored.**
  Target variety `V11 ∩ [1,1,2,2]-shape` (contains 19044.h.2!).  Order 11 at
  infinity = degree pattern `[3,1^8]` — only 8 CF steps (the recorded
  dead end, `pell-cf-order` skill, was order *24* = 21 steps, no anchor
  point, and point-reconstruction rather than ideal-reconstruction).  Two
  attacks: (a) quotient-reconstruction (the Daowsud–Schmidt method — build f
  from symbolic partial quotients, then impose the factorization shape);
  (b) p-adic **ideal interpolation**: Newton-lift the recorded mod-p seeds
  to `Z_p` at high precision on the shape-slice and LLL-recover low-degree
  implicit equations of the surface (recover the *ideal*, not a point —
  evades the height wall that killed the order-24 attempt).  Rehearse the
  whole pipeline on order 7 / `V7` first (4 CF steps, and the `[2,2,14]`
  side benefits immediately).
- **B2: BLP2009.** 17 sporadic order-11 Jacobians (Experiment. Math. 18
  (2009) 65-70) — new seeds (factor types! any 2-rank ≥ 1 member is
  gold) and possibly a rerunnable method with a 2-rank gate.  Agent
  retrieving equations.
- **B3: close the branch-point family.** The two extra-2 conditions on the
  `[22]`-branch are a rank-0 elliptic (dead, proven) and a **genus-3 curve
  searched only to height 500**.  Faltings-finite: run rank/quotient
  analysis (Chabauty if rank ≤ 2) on `c_±(s,x)` from
  `notes/order222_from_order11.md`; ANY nonboundary point is a
  `[2,22]`-candidate; else we get a clean impossibility theorem for this
  component.  Same for the DS-side analogue.
- **B4: second-generation Richelot.** We now have infinite streams of
  `[2,4]`-type members. `polgalois` census of their quartic cofactors: any
  member whose quartic has Galois group ⊆ D4 has a rational (2,2)-kernel;
  Richelot maps it (over Q) to a new 11-curve with a NEW 2-type ((2,2)-isogeny
  preserves End^0 and the odd part).  Iterate the graph from infinitely many
  seeds — yesterday's blocker (finitely many DB seeds, 2-rank 1) is gone.
- **B5: Platonov-school explicit order-11/22 sextics** and other CF degree
  patterns (agent).  Note: order **22** at infinity automatically forces a
  rational factorization `f = A·B` with `[A] = 11·D_inf` — that stratum has
  2-rank ≥ 1 built in, one factor-split away from `[2,22]`.
- **B6: exact geometry of the (4,3) components** (genus via Magma; if
  elliptic with rank, enumerate via MW lattice far beyond the naive sweep
  before re-testing quartic splits).

## Dead/blocked ledger (day 1)

- Branch-point `[22]`-upgrade covers: rank-0 elliptic + genus-3 (≤ h500) —
  `notes/order222_from_order11.md` (pre-existing; genus-3 closure pending B3).
- Quartic-split covers over the QF loci: empty at accessible height (above).
- DB-seed Richelot BFS, blind 2-rank+odd-gate scans: commit 34cfb59.
- Modular/GL(2) constructions: cannot give `End = Z` by design.
- PARI `ellrank` may hang (use `ellanalyticrank` at conductor this size);
  PARI `factor()` drops constants (twist trap above).
- Magma licence can drop after a reboot of the compute box; the local
  licence setup script restores it.

## Compute-policy compliance

Everything above ran locally in minutes (gp; ≤ 3 jobs).  NOT launched, needs
sign-off, all packageable for the big box:

1. Three-root sweep at H ≈ 128 (C or parallel gp; ~10^9 pairs; ~1-2 h on 192
   cores) — probably unnecessary if verification succeeds at H ≤ 32.
2. CF-backward-7/11 symbolic (Magma Groebner; unknown runtime; start small
   locally first).
3. V11 ideal-interpolation harvest (embarrassingly parallel mod-p lifts).

## Addendum (same day): literature intake (agent report) + DS correction

**The repo's "Daowsud–Schmidt family" is Flynn's family**: identically,
`G_{t/4}(x+1) = F_t(x)` (verified symbolically; `scratchpad/dscheck.gp`).
The repo's `G_u` matches the WITHDRAWN arXiv:1708.05511v1 (Lorenzini showed
v1 = Flynn; DS corrected in v2, 2022).  So every "DS" measurement above
duplicates Flynn (explains the mirrored conductor-92 curves), and
`notes/order22_from_order11.md`'s two `[22]` parametrizations are one family.
The genuine second order-11 family is **DS v2** (coefficients with
`(u^5+8)`-denominators; transcribed + CF-verified in
`code/claude_222_dsv2_probe.gp`): generic factor type `[1,5]` for ALL 20,736
members with u-height ≤ 32 — a built-in rational Weierstrass point, 2-rank 0,
generic torsion `[11]` (the `[1,5]` type forces `11E = 0`, so no free 22
here).  Its quadratic-factor / extra-root incidence curves are FRESH targets
(Magma genus job queued) — the v2 analogues of the covers that died on Flynn
are different curves with new chances.

Key items from the 7-torsion/Pell literature agent (full report with URLs and
verification logs in the session transcript; all CF claims re-verified):

- **LPS2004** (Leprévost–Pohst–Schöpp, Abh. Hamburg 74 (2004)): the ONLY
  published order-7 genus-2 families — two 2-parameter families, explicit
  (Thm 5.1: product of two cubics in (p0,q0); Thm 5.2: contains X_0(29) at
  (u,d) = (-57/25, 116/25)).  Free preprint:
  page.math.tu-berlin.de/~kant/publications/papers/LPS-rational-torsion.pdf.
  Backup charts for the three-root game ([3,3]-type needs one cubic fully
  split + one root of the other: heavier than contact-7's linear conditions).
- **Order 14: no published family** (Nicholls' thesis tables confirm); only
  Platonov–Petrunin sextics f14,1 = (x²+1)(x⁴+5x²+4x+4), f14,2 =
  (x²+1)(4x⁴+4x³+5x²+1) (D_inf order 14) and f28,1/f28,2 (28-point over
  D_inf order 14).  All CF-verified.
- **Order-22-at-infinity: does NOT exist in the literature** (explicit
  negative).  Every published order-22 object — Leprévost's
  f22 = (2x²-2x+1)(2x⁴-2x³+x²-4x+4), Petrunin's quintic
  (16x⁴-16x³+x+1)(x+1) (ChebSb 16:4 (2015) Thm 7, S-unit with a degree-1
  finite valuation), both LMFDB [22]s — carries the 22 on a finite-support
  class over an 11-at-infinity class.  Confirms the B1 target should be
  V11 ∩ [1,1,2,2] (11 at infinity + independent 2s), NOT CF-22.  Note the
  deg-5-and-6 periodicity classification is explicitly OPEN (Platonov RMS
  79:6 (2024)); deg ≤ 4 fully classified.
- **BLP2009**: 19 rows = 18 distinct curves (rows 3 = 6 duplicate), shape
  y² = R(x)² - 4c²S(x)², R = x³-x²+ax+b, S = x²+d — f = product of two
  cubics, 11-class = (P)+(P̄)-(∞⁺)-(∞⁻) with u-poly x²+d: a THIRD stratum
  of the 11-locus (Mumford class, not a difference class).  Row 5 =
  X_0(23); **row 4 (a,b,c,d) = (1159/81, 261607/2187, 40/9, 13/27) is an
  uncorrectable typo — fails 11-divisibility at 23/28 good primes**
  (independent erratum).  The norm-equation method behind it (cf. Kronberg's
  2018 thesis, who also shows the analogous order-7 norm-equation family
  collapses to a 4-dim variety with "rather large" equations) could be
  re-run with a 2-rank gate.
- **Fedorov's master equation** (ChebSb 21:1 (2020), free): every genus-2
  order-m torsion class with canonical representative D - (∞⁺+∞⁻) arises
  from ω₁² - ω₂²f = γh^m; solvability criteria + complete parametrization
  for m ≤ 5.  The order-28 INFINITE family (Platonov–Fedorov, Dokl. Math. 98
  (2018), paywalled) came from the deg-h = 2 generalization (machinery free
  in Fedorov ChebSb 19:3 (2018)).  This is the systematic scaffolding for
  route B1 (m = 11, and rehearsal m = 7).
- **RM torsion beyond J_0(23)/J_0(29)** (useful context, LMFDB-verified):
  X_0(39)/⟨w13⟩ = 1521.a.41067.1 has torsion [14]; X_0(67)⁺ = 4489.c.4489.1
  has [11]; 15129.b.408483.1 (newform 123.2.a.c) has [7].
- **[2,2,14] novelty confirmed**: no torsion group containing (Z/2)²×Z/7
  exists anywhere in the literature or production LMFDB, for ANY endomorphism
  type; QM is excluded outright (Laga–Schembri–Shnidman–Voight: QM torsion is
  12-torsion, order ≤ 18); split is excluded in practice by Mazur (an
  elliptic 7-torsion factor has torsion exactly Z/7, so gluings cannot add
  2s on the 7-side: HLP 2000 Table 1 has no 14).  So even the extended-DB RM
  witness is beyond the published record, and a generic witness would be
  doubly new.
- Bonus explicit sextics now in hand (CF-verified): Platonov's f7, f11, f13,
  f14s, f18s, f23, f29, f33 (order 33 at infinity!), f36/f48s, Howe's f70
  with D_inf order 35; Petrunin's ChebSb catalog f20–f39; Elkies' order-32
  family and order-34/39/40 curves.  Seeds for shape-intersection probes.

## Addendum 2 (same day): 11-torsion literature agent + BLP census

Full dossier in the session transcript; verification scripts in scratchpad
(`blprm.gp` promoted to `code/claude_222_blp_rm_census.gp`).  Highlights:

- **BLP2009 fully retrieved** (all 19 table rows verbatim; 18 distinct —
  C̃3 = C̃6).  Their method: NOT a family — a search on the 2-dimensional
  resultant variety Σ = {ord(D_inf) = 11} inside the 4-parameter ansatz
  `y² = R(x)² − 4c²S(x)²`, `R = x³−x²+ax+b`, `S = x²+d` (so f = product of
  two cubics `(R−2cS)(R+2cS)`; the 11-class satisfies D_δ = 5·D_inf).  BLP
  themselves tried and failed to extract a family through X_0(23).
- **Row C4 as printed is wrong** (fails 11 | #J(F_7) = 80; internally
  consistent with the printed Igusa invariants, so a bad stored curve, not a
  transcription slip; no published erratum exists).  The unique correction on
  its (a,c,d)-fiber is **b = −277/243**, giving (integral model)

  ```text
  y² = x⁶−18x⁵−4001x⁴−22524x³+859039x²−1926258x−9043839
     = (x−9)(x+21)(x²−80x+439)(x²+50x+109)
  ```

  with ord(D_inf) = 11 exact, factor type [1,1,2,2] ⇒ 2-rank 2, torsion
  **exactly [2,22]** (44 | #J(F_p) always, #J(F_7) = 44), absolutely simple
  (D4 at p = 7,17,37,53), conductor odd part 645² — a SECOND explicit
  [2,22] curve (beyond staging-DB 19044.h2, outside the DB's range)… but
  Frobenius real-quadratic field constantly Q(√5): **again RM(√5)**.
  *(2026-07-25 erratum follow-up: the "odd part 645²" came from PARI
  `genus2red`, which omits the 2-part.  Magma on the reduced minimal
  Weierstrass model gives the FULL conductor `N = 1664100 = 1290² =
  2²·3²·5²·43²` — a perfect square, exactly as GL2-type/RM(√5) predicts, and a
  consistency check on the 2-part in its own right, since 2^k·645² is a square
  only for even k.  Torsion re-verified `[2,22]`.
  `code/claude_ov_erratum_blp22.m` → `results/claude_ov_erratum_blp22.log`;
  see `notes/claude_ov_erratum_2026_07_25.md`.)*
- **RM census over all 18 BLP curves** (real-Weil-subfield disc kernel
  `core(c3²−4(c2−2p))` across good p ≤ 97; control 1192.a varies as
  expected): only C4-corrected, C5 = X_0(23), and **C6** are constant-{5}
  (RM(√5) signature); the other 15 vary — End=Z signature.  So the BLP
  stratum is NOT inside Humbert-5; its known RM points are sporadic.
  **Route B2′ (new, concrete):** the [1,1,2,2]-locus on Σ (both cubics
  acquire a rational root — a finite cover of Σ with known rational point
  (C4corr, r₁ = 9, r₂ = −21)); rebuild Σ's two defining equations by the
  `F·U² − V² = c·S` elimination (13 coefficient equations, net codim 2),
  then study/parametrize the cover and stream [2,22]-candidates, filtered by
  the disc census (cheap) before full End=Z certificates.  This is now the
  most concrete [2,22] route alongside B1.
- The **real-subfield-disc census is the standard cheap RM-vs-generic
  prefilter** from now on (microseconds per curve via hyperellcharpoly; RM
  witnesses read as a constant set like {5}, generic as a scatter).
- Nicholls (thesis + `phd-code/torsion/found_torsion_curves.m`): explicit
  order-11, order-22 (f = (x²+1)(4x⁴−12x³+13x²−4x+1), type [2,4], D_inf
  order 11), and two order-33 curves (ord(D_inf) = 33!) — more seeds; his
  genus-2 search sets include 22 but NOT 44.
- Extended-DB counts (6,216,959 curves): [11]×171, [22]×37, [33]×1,
  [2,22]×1 (19044.h2, RM(√5)); no [44]/[2,2,11]/[66].
- Context: Alessandrì–Coppola (arXiv:2602.21047) Conjecture 4.5 puts the
  GL2-type abelian-surface torsion orders at {1..24, 28, 31, 37, 44, 56} —
  44 AND 56 are conjecturally possible in the RM world too, consistent with
  our RM witnesses; no order-44 example was known to them (the corrected C4
  and 19044.h2 are two).
- Flynn's 11-family provenance pinned: J. Number Theory 36 (1990) 257–265,
  Application 3.1 (also Invent. Math. 106 (1991) Result 3 with the general
  even-genus family; for g = 2, r = 0: `Y² + t(X−1)Y = (X³+X²+X)² −
  t(X⁴+X)`, exact order 11).

## Addendum 3 (same day): disc census on the six [2,2,14] candidates

Applying the real-subfield-disc census (`code/claude_222_blp_rm_census.gp`
method; ~30 good primes ≤ 149) to the six three-root candidates:

```text
cand1 (-3,-3/4,-3/5):      {2}                       <- RM(sqrt2) SIGNATURE
cand2 (-10,-1/2,-10/7):    {2,7,11,15,29,39,...}     <- End=Z signature
cand3 (-5,-15/8,-15/22):   {2,11,21,46,51,57,65}     <- End=Z signature
cand4 (-1/2,-15/8,-15/19): {2,7,15,22,30,57,65,113}  <- End=Z signature
cand5 (-4/9,4/17,-4/25):   {2,11,15,21,23,37,...}    <- End=Z signature
cand6 (4/17,-5/18,-10/49): {2,14,23,29,65,113,...}   <- End=Z signature
```

So the smallest candidate is (almost certainly) RM by Q(sqrt2) — parallel to
J_0(29) — and would be a *new RM* [2,2,14] witness at tiny height (itself
worth recording; possibly modular of square conductor — check).  The passing
root-power certificates are consistent (root-power proves geometric
simplicity, which RM satisfies; only the disjoint-splitting-fields step
separates End=Z).  **Five End=Z-signatured candidates remain** for the Magma
exact verification; prioritize cand2-cand6 in
`code/claude_2214_threeroot_verify.m` output reading.  Pattern note: the
three-root surface carries sporadic RM points (cand1 here; 152100.eb.2 is
also on the surface after normalization) alongside generic ones — same
structure as the BLP surface on the [2,22] side.

## Addendum 4 (same day): Flynn43 family parametrized; quadratic-split cover killed (conductor-19 rank 0)

The genus-0 (4,3) quadratic-factor component is parametrized
(`code/claude_222_flynn43_param.m`, log committed):

```text
u(tau) = (-11/100 tau^3 + 919/50 tau^2 - 8017/25 tau + 40154/25)/((tau-6)^2 (tau-11))
v(tau) = (2401/1600 tau^4 - ... )/((tau-6)^2 (tau-11)^2)      [see log]
t(tau) = degree-6/degree-6 rational function                   [see log]
```

Structure discoveries:
- The quartic cofactor `Q4` is **reducible over Q(tau)** (Galois group S3 of
  its cubic factor): the family's generic factor type is `[1,2,3]`, i.e. a
  rational Weierstrass root appears for free, and the generic member has
  torsion `[22]` (four Magma-verified specializations, all exactly `[22]`).
- Since disc-census signatures scatter (End=Z-like) at the sampled members,
  this gives a (probable) **infinite generic-[22] family** — the production
  DB had only TWO generic [22] curves.  Worth a verification pass of its own.
- Extra-root cover = tautological genus-0 section + a **genus-3** curve
  (deg 9, absolutely irreducible) — the same shape as the branch-point
  family's blocker; Faltings-finite, unresolved.
- **Quadratic-split cover (the [1,1,1,3] ⇒ [2,22] upgrade): DEAD.**  First
  attempt hit the factor()-drops-content trap AGAIN (the untwisted quartic's
  points produced only [1,2,3]-members — four Magma-verified [22]s); with the
  exact numerator (leading coefficient −3744) the correct model is

  ```text
  Y^2 = -234 tau^4 + 7346 tau^3 - 81284 tau^2 + 336216 tau - 214544,
  ```

  whose Jacobian is the conductor-19 elliptic curve (19.a), **analytic rank 0,
  torsion Z/3**; its only points to height 2·10^5 are tau = 6 and tau = 11 —
  BOTH poles of the parametrization.  So the splitting locus is boundary-only
  (rigorous closure = routine 2-cover bookkeeping against E(Q) = Z/3; queued).

Updated [2,22] state: every low-genus 2-upgrade cover over Flynn's family is
now killed by a rank-0 elliptic curve (branch-point quadratic-split: rank-0
curves, per notes/order222_from_order11.md; QF quadratic-split: rank-0
cond 19; QF (6,2)/(4,3) quartic-splits: empirically empty).  The surviving
in-Flynn hopes are exactly TWO genus-3 Faltings-finite curves (branch-point
cubic-root cover; QF extra-root cover) — both now first-class Chabauty
targets.  Outside Flynn: B1 (CF-backward, RM-anchored), B2' (BLP
[1,1,2,2]-locus), DS-v2 incidence geometry (Magma job still running at
commit time).

## Addendum 5 (same day, post-merge): all parametrized 11-families are 2-upgrade-locked

- **Sieve to parameter height 6000** (three parallel gp jobs,
  `code/claude_222_cubicroot_sieve.gp` machinery, 43.77M parameters per
  family, 30-prime bitmap): the cubic-root condition has NO non-boundary
  rational points on famA± (branch-point families; 12x the old h=500 search)
  or famB (Flynn43-parametrized).  The two genus-3 covers remain pointless.
- **DS-v2 is locked for the opposite reason** (Magma 3-var Groebner hit the
  3GB cap after ~2.5 CPU-h; redone by resultant elimination,
  `code/claude_222_dsv2_qf_resultant.gp`): after boundary factors
  `u^15 (u^5+8)^6`, the quadratic-factor incidence `Res_B(e1,e0)` is a SINGLE
  irreducible curve of bidegree (66,10) in (u,A) — `deg_A = 10 = C(5,2)`, so
  the Galois action on root-pairs of the quintic cofactor is transitive and
  no low-genus stratum splits off.  The extra-root locus is likewise a single
  irreducible (29,5) curve (full quintic transitivity).  Unlike Flynn, whose
  incidence split into a genus-0 and an elliptic component, v2 has no
  tractable strata at all: generic member stays [1,5], torsion [11].
- **Consolidated [2,22] verdict:** upgrading any known parametrized 11-family
  to 2-rank 2 is structurally blocked everywhere we have looked — Flynn's
  low-genus covers are all rank-0-dead (three independent rank-0 kills:
  cond-54-class quadratic-split, cond-19 QF quadratic-split, rank-0 twist on
  (6,2)-quartic-split... plus empirical emptiness), its two genus-3 covers
  are pointless to height 6000, and v2 has no low-genus covers to try.  The
  live strategy is therefore to BUILD the 2-structure in from the start:
  B1 (CF-backward on the [1,1,2,2] ∩ D_inf-order-11 surface, anchored at RM
  witness 19044.h.2; rehearse with order 7 = 4 CF steps), B2' (the BLP
  resultant surface's [1,1,2,2]-locus, anchored at corrected C4), plus cheap
  insurance (C-port sieve to 1e5-1e6) and Chabauty closure of the two
  genus-3 curves (would make the Flynn lock a theorem).

## Addendum 6 (evening): Flynn lock finished at all depths; the [2,2,14]-family program

**[2,22] / Flynn: the lock is now complete at every accessible depth.**

- Genus-3 blocker structure (code/claude_222_genus3_*.m): both X1 (branch-point
  cubic-root cover) and X2 (Flynn43 extra-root cover) have exactly one
  involution; both quotients are genus 1 with Jacobian THE SAME curve:
  y² = x³ − x + 1, conductor 92 = 4·23, rank 1, trivial torsion.  This one
  elliptic curve also underlies both (6,2)-components — the whole Flynn
  2-upgrade tower runs on a single level-92 engine (23 = the X_0(23)/Flynn
  level; the thinness looks modular in origin).  X1's Prym is
  generic-signatured (scattered real-subfield discs), so no RM/rank shortcut.
- **Mordell–Weil fiber test** (code/claude_222_mw_pullback.m, final verified
  version): enumerating |n| ≤ 60 multiples of the rank-1 generator and testing
  each fiber of the double covers exactly: X1 contains ONLY its two boundary
  points (s,x) = (1,1), (0,0) (at n = −5, −2); X2 ONLY its pole point
  (τ,r) = (11,1) (at n = 3).  Parameter heights ~e^{3600·ĥ} — hundreds of
  digits, far beyond any sieve.  METHODOLOGY WARNING recorded: Magma's
  Pullback through CurveQuotient/EllipticCurve maps picks up patch-
  indeterminacy junk (base points appear in EVERY fiber; and point-transport
  via naive pullbacks silently returns the base point) — verify images by
  evaluation at EVERY transport stage.
- **C sieve at height 10⁶** (60 primes; code/claude_222_cubicroot_sieve.c,
  validated against the gp pilot 14/14/10): all three families CLEAN — the
  40/40/38 reduced survivors are exactly the boundary orbits (famA: par 0, 1;
  famB: par 6, 11, 58/3) plus cubic-irrational stragglers; zero rational
  cubic roots (results/claude_222_cubicroot_sieve_h1e6.log).
- Verdict: [2,22] will not come from upgrading Flynn.  Remaining routes
  stand: B1 (CF-backward with 2-structure built in), B2' (BLP locus),
  proof-grade closure of X1/X2 (quadratic-Chabauty-class).

**[2,2,14] infinite-family program (B1 rehearsal) — first structural results:**

- The RM anchor's CF pattern in the [1,1,1,1,2]-chart is **[3,1,2,1]** (all 12
  root-normalizations; pattern is x-affine-invariant).  Mod-p census at
  p = 101 (code/claude_v7_pattern_census.gp): among 413 order-7 shape points,
  62% have generic pattern [3,1,1,1,1], 38% the anchor pattern [3,1,2,1] —
  two strata, both substantial.
- Stage A (code/claude_v7_stratum2.m): the deg-drop hypersurface E0 :=
  lc(Q_2) has 581 terms, degrees (10,10,10,5) in (a,b,c,d), IRREDUCIBLE,
  and contains all 12 anchor points (data/claude_v7_stratum2_eqs.txt).
  The 12 anchors affinely span A^4 (no common hyperplane).
- Stage B (code/claude_v7_stratum2b.m, running): CF over the function field
  Q(a,b,c)[d]/(E0) with the degree drop built in; the closure coefficients of
  Q_4 cut out the anchor stratum S.  Rational curves on S through the anchors
  = the infinite generic [2,2,14] family.  The generic-stratum derivation
  (code/claude_v7_cfbackward.m) is still computing in parallel (its closure
  pair would give the second, denser stratum).
- Negative worth recording: anchor-centered coordinate-line slicing of V(E0)
  (12 anchors x 3 slices x 4407 sweep values, all rational d-fibers CF-tested)
  found ZERO new S-points — as it must: such slices meet the codim-1 surface
  S in finitely many, generically irrational points.  S(Q)-harvesting needs
  S's equations (the running symbolic jobs) or mod-p interpolation of the
  closure polynomial (basis: deg_d <= 4 representatives; S(F_p) ~ p^2 points
  per prime from pattern-filtered hypersurface scans — feasible if
  deg(E_c) <= ~20).  Next session: whichever of (symbolic completion |
  interpolation) is ready first, then curve-fitting through the 12 anchors.

## Files

- `code/claude_2214_threeroot_sweep.gp` (+ `results/claude_2214_threeroot_h16.log`, `_h32.log`)
- `code/claude_2214_threeroot_surface.gp` (R, e-basis form, hit checks)
- `code/claude_2214_threeroot_verify.m` (Magma exact verification, READY)
- `code/claude_222_qf_incidence.gp` (+ `results/claude_222_qf_incidence.log`) — incidence components, witness CF table
- `code/claude_222_flynn62_enum.gp` (+ log) — elliptic-point member enumeration + split test
- `code/claude_222_comp43_sweep.gp` (+ `results/claude_222_qf_covers.log`) — (4,3) sweeps + DS62 disc
- `code/claude_cf_order_census.gp` (+ `results/claude_cf_order_census.log`) — mod-p densities + seeds
- this note.
