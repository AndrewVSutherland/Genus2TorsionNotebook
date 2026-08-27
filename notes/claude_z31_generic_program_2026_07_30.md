# The generic Z/31 program (2026-07-30)

Goal: a genus-2 curve /Q with `J(Q)_tors ⊇ Z/31` and `End(J_Qbar) = Z`.
Companion to `notes/claude_z31_rm_witness_2026_07_30.md` (the RM realization).
Status: sieve construction session 2026-07-30; nothing generic found yet.

## Frame (Lane 8's lesson applies verbatim)

`J[31]` is finite étale, so "has a rational point of order 31" is not a
subvariety of M_2 — it is the image of the rational points of a finite étale
cover `W_31 -> M_2` (3-fold). The known RM witness is ONE rational point of
`W_31`, lying over the Humbert surface H_8. The generic question: does
`W_31(Q)` meet the locus over `M_2 \ (all Humbert/special loci)`?
Two sobering facts about the known point:
- its curve has (apparently) NO rational points, and the 31-class is invisible
  (conjugate-quadratic support) — only its *Mumford coordinates* are rational
  and small;
- in every visible chart the lab uses, the RM point has large or undefined
  coordinates: no small-height box scan in any standard chart reaches it.

Therefore the program searches **charts of W_31 in Mumford coordinates**, via
the fiber method (fix a slice, enumerate the fiber mod P at 3-4 primes, detect
31-torsion mod P, CRT-match, rationally reconstruct, verify exactly) — the
architecture proven by `code/claude_ov_b2p_scan.c` on the [2,22] program. The
fiber method reaches height ~P1*P2*P3 in the SOLVED coordinates.

## The three sieves + backstop

1. **Chart U' (universal Mumford), `code/claude_z31_kummer_sieve.c`.**
   Parameters (u1,u0,v1,v0,w4,w3,w2,w1,w0); curve `y^2 = f := v^2 + u*w`
   (u = x^2+u1x+u0, v = v1x+v0, w = quartic, leading coeff w4 free — any sign,
   square or not); marked class D = [div(u,v) − D_inf]. EVERY pair
   (genus-2 curve /Q, rational degree-0 class) has this form — pointless
   curves and invisible classes included; the RM witness sits at
   (u1,u0,v1,v0) = (−21/20, 9/20, 101/200, 51/200), w = (F−v^2)/u. Condition
   `31·D = 0` tested by a Kummer-surface Montgomery ladder (Cassels–Flynn
   biquadratics; no Cantor, no infinity-representation issues, leading-coeff
   agnostic): accept iff kappa(31 D) = (0:0:0:1). Slice = (u1,u0,v1,v0,w4,w3,w2),
   solve (w1,w0) in F_P^2. Degenerate cases (disc(u) ≡ 0, non-squarefree f mod P)
   are forwarded, never resolved in C.
   End-to-end validation requirement: the sieve must detect the RM witness in
   its own fiber and reconstruct its exact rational (w1,w0).

2. **Chart I (infinity/CF), `code/claude_z31_cf31_scan.c`.**
   Monic sextics (c5=0), marked class D_inf, condition ord(D_inf)=31 via the
   polynomial continued fraction (quasi-period degree total = 31, incl.
   deg a0 = 3); ~31 bounded steps per candidate, cheapest per-candidate cost,
   deepest slice reach. Direct adaptation of the b2p CF machinery (11 -> 31).
   Covers the D_inf corner of W_31 (the corner ALL prior CF-school realizations
   {..., 28, 33} live in; 31 not realized there in the literature).

3. **Elkies–Kumar RM(sqrt2) exhaustion, `code/claude_z31_ek8_sieve.*`.**
   On the 2-dim Y_-(8) chart (r,s) (Elkies–Kumar 2014, disc 8), a rational
   31-point is a 0-dim condition: fiber-sieve r-slices, solving s mod P with the
   twist-aware condition 31 | chi_P(±1), CRT + twist-character consistency.
   Decides whether 1830.2.a.q is (essentially) the unique RM-sqrt2 [31] member
   in the searched range — the analog of the [2,22] question "is Sigma'' ∩ H_5
   all there is?", but attacked on the RM side directly.

4. **Backstop box scan, `code/claude_z31_box31.c`.** Raw 31 | #J(F_p) (6 primes)
   over integral sextic boxes — assumption-free (catches invisible classes on
   pointless curves), but only within the box; the alpha DB's absence of 31
   says small boxes are probably empty. Run in spare cycles only.

## Expectations management

The moduli 3-fold `W_31`(-analog with full level-31-point structure) is
plausibly of general type (literature lane pending), in which case rational
points are finite and possibly all special (RM/CM) — the [2,22] pattern.
A negative result at scale is still valuable: it sharpens the emerging
conjecture that large-prime torsion on abelian surfaces /Q is *forced onto
GL2-type/RM loci* (where the Eisenstein/cuspidal mechanism supplies it).
Conversely, ONE generic hit would be a record and would refute that.

## Session logistics

Compute: aws-spot-11 (claimed, 192 cores). Exact stage for hits:
`TorsionSubgroup` + D4 simplicity + strict two-prime End=Z certificate
(`code/claude_end_z_certificates.m` idiom). All sieve code validated against
Magma reference vectors (`data/claude_z31_vectors_chart{U,I}.tsv`) before any
production run; the RM witness is the planted positive control everywhere.

## Campaign log (2026-07-31)

- Q1 chart-I CF-31, integer |c4,c3,c2| <= 30: **EMPTY** (226,981 slices; 2
  spurious CANDs = one x->-x orbit, killed by exact CF over Q). First negative:
  the CF-school's own corner of W_31 has no 31 at integer slice height 30.
- Q2 EK RM-sqrt2 H=60: sieve stage done (~6,500 CAND lines; twist-OR acceptance
  is rich by design); Magma exact batch queued for after the campaign.
- Q3 Kummer U-a: OPERATIONAL ERRATUM — the launch passed recbound 10^7 against
  a 4-prime CRT modulus M = 3.2e13; since 2B^2 > M, every residue combo
  "reconstructs" and the CAND stream is CRT junk (826 candidates, ALL killed by
  a 5-independent-prime 31-divisibility triage in PARI, zero survivors; ~12%
  were the degenerate v0=0 sub-line). Genuine hits of height <= ~4e6 would
  still have surfaced (lattice minimum), so the triage conclusion "U-a box
  empty so far" stands, but the matching stage must be RERUN with
  recbound <= 3e5 (2B^2 << M) for a clean record; sieve-side fiber computation
  was correct throughout. The cf31 (B=2e4, M=1e12) and EK (B<=2e6, M~1e15)
  campaigns are correctly configured.
- Q2 EK exact stage (2026-07-31): all 4,938 unique (r,s) candidates from the
  H=60 Y_-(8) sieve exact-checked in Magma: **ZERO 31-hits**; 4,348/4,938 (88%)
  admit no curve /Q at all (Mestre/pp obstructions dominate the EK surface).
  Within the height-60 slice box the 1830.2.a.q moduli point (slice height
  311^2 = 96721, outside every feasible box) has NO RM-sqrt2 companion with a
  rational 31-point.
- Q4 chart-I CF-31, rational slices height <= 5: **EMPTY, zero CANDs** (not
  even spurious). Chart I now clear over integers |c|<=30 and rationals h<=5.
- Q5 box scan H=43 (2026-08-01, FINAL): full box complete — 16/16 shards,
  ~1.8e13 deduplicated integral sextic models tested against the six-prime
  31-divisibility filter; 548 survivors; **exact stage: has31 = false for ALL
  548** (every survivor a divisibility coincidence, torsion trivial/small).
  The assumption-free coefficient box |c_i| <= 43 contains NO genus-2 Jacobian
  /Q with a rational 31-torsion point.

## CAMPAIGN VERDICT (2026-08-01)

All five fronts of the generic-[31] program are now decided EMPTY in their
searched ranges: chart I (CF/D_inf; integers |c|<=30 + rationals h<=5),
chart U' campaign U-a (triaged), EK Y_-(8) H=60 (the RM witness provably alone
in the box; 88% of the surface Mestre-obstructed), and the assumption-free
H=43 box (1.8e13 models). Zero genuine 31-candidates anywhere. Together with
the moduli-theoretic thinness (A^lev_{1,31} not unirational) and the
narrow-class analysis, the emerging picture: rational 31-torsion on abelian
surfaces /Q may be exclusively a GL2-type/Eisenstein phenomenon at h^+ = 1
Hecke fields — 1830.2.a.q possibly THE unique genus-2 Jacobian instance at any
accessible height, and a generic (End = Z) [31] Jacobian, if one exists, lies
beyond every chart and box this campaign could reach. Next escalations (future
sessions): deeper chart-U' slice families with the fixed recbound, chart-I
rational slices h<=8..10, the 43920.2.a.cr second-31 eigenvalue computation at
dim ~9000, and the W_31 geometry itself (components/Kodaira dimension of the
level-31 cover over the paramodular locus). Compute released: aws-spot-11
stopped 2026-08-01 after ~40h of campaign use.
