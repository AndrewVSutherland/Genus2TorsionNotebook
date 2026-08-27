# B4: the Elkies-Kumar disc-8 (RM sqrt2) 31-torsion sieve (2026-07-30)

Lane 3 of `notes/claude_z31_generic_program_2026_07_30.md`: exhaust the
RM-by-Q(sqrt2) surface for rational 31-torsion and decide whether 1830.2.a.q
is (essentially) the unique disc-8 member in a searched box.

## The model (Elkies-Kumar, disc 8)

Source: the arXiv e-print tarball of arXiv:1209.3527 IS the ancillary data
(one directory per discriminant); files for D=8 recorded in
`data/claude_z31_ek8_source/` (8.txt, igusa8.txt + SOURCE.md with URL/sha).
The rational (r,s)-plane maps to Igusa-Clebsch invariants

```text
I2  = -4(3s+8r-2)                        deg (r,s) = (1,1)
I4  =  4(9rs+4r^2+4r+1)                  deg (2,1)
I6  = -(16/3)(9rs+4r^2+4r+1)(3s+8r-2) + (4/3)(54r^2 s+81rs-16r^3-24r^2-12r-2)
                                         deg (3,2)
I10 = -8 r^3 s^2                         deg (3,2)
```

(simplified from igusa8.txt's `[-24 B1/A1, -12 A, 96 (A/A1) B1 - 36 B,
-4 A1 B2]`; equality of the two forms asserted in Magma —
`IG_SIMPLIFICATION_OK` in `results/claude_z31_ek8_witness.log`.)
Y_-(8) itself is the double cover z^2 = 2(16rs^2+32r^2 s-40rs-s+16r^3+24r^2
+12r+2), but the moduli point depends only on (r,s), so the sieve works on the
plane. Degenerate loci: rs = 0 (I10 = 0) and the special-automorphism locus
(detected as Mestre-conic det = 0).

## The witness lives at height 96721 in this chart

`code/claude_z31_ek8_witness.m` computes the IC invariants of the RM witness
(integral model y^2 = -3356x^6 + ... - 504) and solves the weighted-projective
match on the chart (Groebner, saturated by rs != 0):

```text
(r, s) = (6000/96721, 3557520/96721),   96721 = 311^2   — UNIQUE preimage
r = 2^4*3*5^3 / 311^2,  s = 2^4*3^6*5*61 / 311^2   (61 = the Eisenstein prime
of level 1830; 311 = 10*31 + 1 = 1 mod 31)
```

Verified: exact weighted IC match (lambda^2 = 1/1547536), conic det != 0,
rebuild over Q via HyperellipticCurveFromIgusaClebsch reproduces the moduli.
**Consequence:** the natural slice box |height(r)| <= 60 does NOT contain the
witness; RM points of modular origin can hide at slice-height ~10^5 and
solved-height ~4*10^6. A box sweep answers emptiness of the box only; the
end-to-end validation instead plants the exact witness slice.

## Sieve architecture (`code/claude_z31_ek8_sieve.c`)

Fiber sieve in the b2p style, slice = r (rational), solved coordinate = s.
Per candidate (r,s) in F_P (P = 3 mod 4, P > 5, P != 31):

1. IC invariants (formulas above) -> Clebsch (A,B,C,D) (Magma igusa.m
   conversion) -> Mestre conic A_ij + cubic a_ijk (Magma mestre.m formulas;
   the van Wamelen U/DP rescalings are twist normalizations over Q and are
   omitted — mod P we accept both twists).
2. Deterministic conic point (x1 = 0,1,2,... , sqrt via a^((P+1)/4)),
   parametrization, cubic evaluation -> sextic f(m).
3. Exact chi via N1 = #C(F_P), N2 = #C(F_{P^2}) (F_{P^2} = F_P[i], i^2 = -1;
   norm-character lookup; degree-6 finite differences in the real part).
   Accept iff 31 | chi(1) or 31 | chi(-1)  — twist-class invariant.
   Twist mask recorded per prime (bit0/bit1); DEGEN (rs=0, det=0, no conic
   point, deg f < 5, non-squarefree f) forwarded as wildcard mask 7.
4. Sweep: per-prime fiber rows cached by r mod P (needed rows only),
   CRT product across npr in {4,5,6} primes with hoisted partial CRT,
   rational reconstruction of s (bound = recbound), then membership check at
   nver extra verify primes.  Work flattened (slice, j0[, j1]) for full-core
   utilization even on a single slice.

Twist-sign consistency across primes is NOT enforced in C: the local Mestre
model carries no canonical twist, and any sign pattern on <= 6 primes is
realizable by some quadratic character, so the constraint has no sieve power;
the global twist d is recovered in the exact stage instead.

## Exact stage (`code/claude_z31_ek8_exact.m`)

CAND (r,s) -> curve /Q from IC -> reduced minimal model -> chi_p(+-1) mod 31
at ~30 good primes: each prime with exactly one working twist class
contributes a constraint chi_d(p) = eps_p; a prime with NO working class
kills the candidate (DEAD_AT_p). Squarefree d supported on bad primes + p <=
20 enumerated against the constraints; surviving d (ordered by |d|) get
TorsionSubgroup(Jacobian(y^2 = d*F0 minimal-reduced)). Validated on the
witness: the constraints pinned the UNIQUE twist d (52-digit, matching the
Mestre model's junk support) out of ~2^16 candidates; TorsionSubgroup on the
reduced twist confirms Z/31 (HIT31 marker) — see
`results/claude_z31_ek8_exact_witness.log`.

## Validation summary (all on aws-spot-11 unless noted)

- `IG_SIMPLIFICATION_OK` — chart formulas equal EK's raw ones as rational
  functions.
- Magma reference vectors (`code/claude_z31_ek8_ref.m` ->
  `data/refvec_ek8_P{103,211,1019,1031}.tsv`): C `ref` mode PASS 734/734
  (250+245+120+119), 0 FAIL, 0 skip; DEGEN classification agrees exactly
  (Magma-side criterion identical: I10 = 0 or conic det = 0; 4 DEGEN rows).
- `selftest`: witness fiber anchors at P in {103,211,1019,1031,1039,1051,
  1063,1087}: 8/8 PASS — C chi multisets {chi(1),chi(-1)} equal Magma
  LPolynomial values at the witness reduction, masks nonzero at every prime.
- End-to-end: `sweepr 6000 96721` with CRT primes {1019,1031,1039,1051,1063,
  1087}, verify {1103,1123,1151,1163}, recbound 4*10^6 reconstructs
  `CAND r=6000/96721 s=3557520/96721` exactly
  (`results/claude_z31_ek8_witness_e2e.log`), among spurious near-recbound
  reconstructions (expected ratrecon leakage at B^2/M ~ 1.5e-5/leaf).
- Exact stage on that CAND: unique d, TORSION [31], HIT31.

## Throughput (measured)

- Per-candidate: 3.6 ms at P = 1019 on aws-spot-11 (281 candidates/sec/core;
  local Ryzen: 220/s/core); cost ~ P^2 (N2 count dominates). Accept rate
  0.071-0.081 =~ 2/31 with a mild GL2-type enrichment; mean fiber 89-91 + ~3
  DEGEN per prime at P in [1019, 1223].
- Fiber rows (P candidates each): ~3.7 core-s/row (spot-11), 4.7 (local).
  Full 12-prime row set for the H=60 box (~13.3k rows): ~4.6 min wall on 192.
- Match stage (optimized per-level-inverse CRT, clean local measurements,
  `results/claude_z31_ek8_meas_local.log`):
  deep config (5 CRT + 7 verify, recbound 2e6): 1.757e10 combos / 61.0 s wall
  on 16 phys cores = 56 ns/leaf (~43 ns projected Zen4);  pass-1 config
  (4 CRT + 8 verify, recbound 7e5): 6.64e8 combos / 3.8 s = 92 ns/leaf.
  Earlier un-optimized shard on the (contended) spot box:
  `results/claude_z31_ek8_sweep_meas.log` (~550 ns/leaf).
- Spurious CAND volume: deep config ~0.9-2.5/slice; pass-1 ~0.5/slice
  (incl. the special-locus wildcard points, e.g. the det=0 rational points
  found at r=-60).  Witness e2e run: 5.14e11 combos -> 430 CAND lines, the
  witness among them (`results/claude_z31_ek8_witness_e2e.log`).
- Projections for aws-spot-11 (192 cores, exclusive):
  pass 1: ~16 core-s/slice -> ~6 min per prime set for the full H=60 box
  (~1.0M slices/day, fiber-bound);
  deep pass: ~200 core-s/slice -> ~1.3 h per prime set (~83k slices/day).
  Full two-pass x two-prime-set campaign: ~2.8 h wall.

## Caveats

- A solved s whose denominator is divisible by a CRT prime is invisible at
  that prime (fiber lists affine residues only): production should re-run
  the box with a second disjoint prime set (same cost).
- Slices with r = 0 mod P or den(r) = 0 mod P skip that prime; the sweep
  counts them (`skipped`) — rerun those slices with alternate primes.
- The witness's own slice (r-height 96721) is outside any feasible slice box;
  box sweeps answer the box question only. Smarter slicing (e.g. square-
  denominator families suggested by 311^2, or a Cremona change of chart) is
  an open improvement.
- Exact stage factors disc of the raw Mestre model; a hard 100+-digit
  semiprime there would stall it (did not occur on the witness; add
  TrialDivision fallback if it bites).
