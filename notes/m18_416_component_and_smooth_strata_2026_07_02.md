# M_1(8,4) [4,16]: projected components and smooth p=7 strata, 2026-07-02

Goal: follow up the p=7 blowup/smoothness work for the `[4,16]` target in
the `M_1(8,4)` family along two fronts:

1. isolate the repeated low-degree component in the projection
   `Res_b(E1core,E0core)` of the reduced second-halving surface
   `Sigma' : E1core = E0core = 0` in `A^4_{R,w,a,b}`;
2. run a rational search restricted to the six p=7 strata where the
   implicit-function smoothness scan found genuine `Q_7` points:

```text
<3,3>, <3,4>, <4,0>, <5,0>, <5,2>, <5,5>.
```

This note is meant to be a working log.  It deliberately records negative
checks too, so future agents can avoid repeating the same blind computations.

## Starting state

The relevant recent files before this run were:

```text
code/agent_m18_416_e1e0_reduction.m        (untracked)
code/agent_m18_416_branch_discriminant.m   (untracked)
data/agent_m18_416_branch_discriminant.log (untracked)
notes/agent_m18_416_p7_blowup_notes.md
code/agent_m18_416_live_stratum_search.m
```

The smoothness note supersedes the earlier "no Hensel-smooth point" language:
the fixed-`(R,w)` aux-Jacobian test was too strong.  Letting `(R,w)` move,
the `[4,16]` cover has smooth `Q_7` points on six boundary strata, so the
current absence of rational hits through height 70 is a height/global-geometry
question, not a p=7 local obstruction.

## Planned computations

### A. Projected branch/resultant component

The existing `branch_discriminant` log found, modulo both `101` and `103`,
the same factor profile for

```text
Res_b(E1core,E0core) in Q[R,w,a]:
  (w-1)^3 (w+1)^3 (R+1)^8 * (degree 5)^8 * (degree 52).
```

The degree-5 factor with multiplicity 8 is the concrete object to isolate.
The next script should compute the rational resultant or, if that is too
large, compute and print the degree-5 factor modulo several primes and compare
the lifted shape.

## Symbolic result: the repeated degree-5 factor is `d4=0`

New scripts/logs:

```text
code/agent_m18_416_resb_component.m
data/agent_m18_416_resb_component.log
data/agent_m18_416_resb_component_d4saturated.log
code/agent_m18_416_degree5_component_gcd.m
data/agent_m18_416_degree5_component_gcd.log
```

The rational resultant was easy enough to factor over `Q`:

```text
Q: Res_b degree 106 terms 36154
Q Res_b factor profile:
  (w-1)^3
  (w+1)^3
  (R+1)^8
  F5^8, where F5 has degree 5 and 15 terms
  F52
```

The degree-5 factor is

```text
F5 =
R^5 + 1/4*R*w^4 - R^4 + 3/2*R^2*w^2 + 3/4*w^4
- R*w^2*a - 2*R^3 + R*w^2 - w^2*a
+ 1/2*R^2 - 3/2*w^2 + R*a - 1/4*R + a - 1/4.
```

But this is not a new geometric branch of the valid square-quartic cover.
The leading coefficient condition in the reduced square test requires
`d4 != 0`, and the cleared leading coefficient factors as

```text
d4 = (R-1) * 4*F5.
```

After the previously known `R-1` factor is stripped, the entire repeated
degree-5 resultant factor is exactly `d4core^8`:

```text
Q: d4core^8 divides Res_b? true
Q: Res_b / d4core^8 degree 66 terms 7501
Q d4-saturated Res_b factor profile:
  (w-1)^3
  (w+1)^3
  (R+1)^8
  F52
```

After also stripping the obvious boundary factors `w=+-1` and `R=-1`, the
remaining projected component is a single factor:

```text
degree 52, terms 4253, multiplicity 1
degrees in (R,w,a) = (42,28,16)
```

I also tested the graph `F5=0`, i.e. solved for `a` and substituted into
`E1core,E0core`.  Their gcd in `b` has degree `3`, but this lives on the
excluded `d4=0` locus: `E1core` is entirely absorbed by that gcd and `E0core`
only leaves a linear quotient.  This confirms the degree-5 factor is a
leading-coefficient artifact, not the intended `<0,0>` square-root branch.

Conclusion: the projected resultant component to study after saturation is
the remaining degree-52 factor, or else one should return to the localized
`E416`/residual-`G` machinery at `<0,0>`.  The degree-5 factor should be
ignored in future rational searches unless one deliberately studies the
degenerate `d4=0` boundary.

### B. Smooth-strata rational search

The existing `agent_m18_416_live_stratum_search.m` uses the full 449-residue
mod-49 survivor gate.  The narrower search should preserve its exact chain:

```text
FirstCoverPossible -> TangentCandidates -> exact Tx half -> exact P_R half -> torsion.
```

The difference is the CRT gate: include only mod-49 residues above the six
smooth mod-7 strata, optionally with a second mode that keeps only the known
mod-343 survivors for `<5,0>` and related refined strata when available.

## Smooth-strata search setup

I patched

```text
code/agent_m18_416_live_stratum_search.m
```

with an optional flag:

```text
smooth_strata_only:=true
```

The default behavior is unchanged.  With the flag enabled, the script still
recomputes the 449-residue mod-49 blowup gate, then keeps only residues above
the six smooth mod-7 strata.  A height-8 validation run gave:

```text
mod-49 survivor gate: 449 residues
smooth-strata-only gate: 171 residues over 6 mod-7 strata
  <3,3>: 23
  <3,4>: 23
  <4,0>: 28
  <5,0>: 49
  <5,2>: 24
  <5,5>: 24
height=8, pair budget=433
DONE cand=188 auxkill=228 smooth=168 firstposs=1 tangent=0 exact=0 firstver=0 prhalf=0 hits=0
```

I then launched a height-100 search in six parts:

```text
magma -b height:=100 aux_primes:="11,13" smooth_strata_only:=true \
  NParts:=6 Part:=i progress_interval:=50000 \
  code/agent_m18_416_live_stratum_search.m \
  > data/agent_m18_416_smooth_strata_h100_part{i}.log
```

Completed aggregate:

```text
height                         100
parameters                     12175
mod49 smooth-strata gate       171 residues
mod49-gated pair budget        8081977
candidates after aux filters   4951090
aux-killed                     3128508
smooth curves                  4947498
first-cover possible           177
tangent bases                  5
exact tests                    5
first halves verified          5
P_R halves                     0
hits                           0
```

No `PR_HALF` or `HIT_416` lines occurred in any of the six logs.

Per-part summaries:

```text
part0: cand 812012, auxkill 512588, smooth 811727, firstposs 21, tangent 0, exact 0, firstver 0, prhalf 0, hits 0
part1: cand 841323, auxkill 531002, smooth 840816, firstposs 35, tangent 0, exact 0, firstver 0, prhalf 0, hits 0
part2: cand 839327, auxkill 532270, smooth 838102, firstposs 27, tangent 2, exact 2, firstver 2, prhalf 0, hits 0
part3: cand 837766, auxkill 532589, smooth 837194, firstposs 32, tangent 3, exact 3, firstver 3, prhalf 0, hits 0
part4: cand 809227, auxkill 506124, smooth 808941, firstposs 25, tangent 0, exact 0, firstver 0, prhalf 0, hits 0
part5: cand 811435, auxkill 513935, smooth 810718, firstposs 37, tangent 0, exact 0, firstver 0, prhalf 0, hits 0
```

Conclusion: the smooth p=7 strata are locally real but rationally very sparse.
Even on these positive local strata, the same obstruction appears at the exact
second-halving step: whenever a rational first half was verified, `P_R` still
failed to halve.  A blind height increase is possible, but this motivated the
mod-343 gate recorded below, especially because the `<5,0>` stratum already
had level-2 data available.

## C. Smooth-strata mod-343 gate and height-100 rerun

New script/log/data:

```text
code/agent_m18_416_smooth_strata_mod343_gate.m
data/agent_m18_416_smooth_strata_mod343_gate.log
data/agent_m18_416_smooth_strata_mod343_gate.txt
data/agent_m18_416_smooth_strata_mod343_h100_part0.log
data/agent_m18_416_smooth_strata_mod343_h100_part1.log
data/agent_m18_416_smooth_strata_mod343_h100_part2.log
data/agent_m18_416_smooth_strata_mod343_h100_part3.log
data/agent_m18_416_smooth_strata_mod343_h100_part4.log
data/agent_m18_416_smooth_strata_mod343_h100_part5.log
```

I added a standalone gate emitter:

```text
magma -b output_file:=data/agent_m18_416_smooth_strata_mod343_gate.txt \
  code/agent_m18_416_smooth_strata_mod343_gate.m \
  > data/agent_m18_416_smooth_strata_mod343_gate.log
```

The script reuses the exact `E416` lift test from
`agent_m18_416_p7_blowup_level2.m`, but restricts to the six smooth mod-7
strata.  It writes one surviving `(R,w) mod 343` residue per line, followed
by its underlying `(R,w) mod 7` stratum.  Output format:

```text
# comments
R W R0 W0
```

As a validation, I first ran only the `<5,0>` stratum:

```text
magma -b only:="5,0" \
  output_file:=data/agent_m18_416_smooth_strata_mod343_gate_50_test.txt \
  code/agent_m18_416_smooth_strata_mod343_gate.m \
  > data/agent_m18_416_smooth_strata_mod343_gate_50_test.log
```

It reproduced the known level-2 count:

```text
stratum <5,0>: aux7=6 survive=1715 capped=0 killed=686 total=2401
DONE survive=1715 capped=0 killed=686
```

I also checked that the live search could parse an external mod-343 gate file
with a tiny height-8 run:

```text
height=8, aux_primes=11,13, gate_mod=343
gate_residues=1715, mod343-gated pair budget=98
DONE cand=40 auxkill=50 smooth=40 firstposs=0 tangent=0 exact=0 firstver=0 prhalf=0 hits=0
```

The full six-strata mod-343 gate has no capped/unknown cases:

```text
<3,3>: survive 766, capped 0, killed 1635, total 2401
<3,4>: survive 766, capped 0, killed 1635, total 2401
<4,0>: survive 833, capped 0, killed 1568, total 2401
<5,0>: survive 1715, capped 0, killed 686, total 2401
<5,2>: survive 806, capped 0, killed 1595, total 2401
<5,5>: survive 806, capped 0, killed 1595, total 2401

total survive 5692, capped 0, killed 8714
```

The gate file has 5695 lines: three header/comment lines plus 5692 live
residues.

I then patched `agent_m18_416_live_stratum_search.m` so the existing exact
search chain can read an external residue gate:

```text
gate_mod:=343
gate_file:=data/agent_m18_416_smooth_strata_mod343_gate.txt
```

If `gate_file` is supplied and `gate_mod` is omitted, the script defaults to
`gate_mod:=343`; otherwise the original mod-49 behavior remains the default.
The script buckets rational `(R,w)` pairs modulo the active gate modulus, so
the same path also supports any future residue file at another modulus.  After
the height-100 run, I also made the external-gate path skip the unused inline
mod-49 blowup enumeration; the same height-8 parser check still returns
`gate_residues=1715`, `mod343-gated pair budget=98`, and `DONE cand=40`.

Height-100 run:

```text
magma -b height:=100 aux_primes:="11,13" gate_mod:=343 \
  gate_file:=data/agent_m18_416_smooth_strata_mod343_gate.txt \
  NParts:=6 Part:=i progress_interval:=50000 \
  code/agent_m18_416_live_stratum_search.m \
  > data/agent_m18_416_smooth_strata_mod343_h100_part{i}.log
```

Completed aggregate:

```text
height                         100
parameters                     12175
mod343 smooth-strata gate      5692 residues
mod343-gated pair budget       5495485
candidates after aux filters   3365620
aux-killed                     2128248
smooth curves                  3363112
first-cover possible           121
tangent bases                  4
exact tests                    4
first halves verified          4
P_R halves                     0
hits                           0
```

No `PR_HALF` or `HIT_416` lines occurred in any of the six mod-343 logs.

Per-part summaries:

```text
part0: cand 564454, auxkill 356474, smooth 564001, firstposs 17, tangent 1, exact 1, firstver 1, prhalf 0, hits 0
part1: cand 563972, auxkill 356011, smooth 563682, firstposs 16, tangent 0, exact 0, firstver 0, prhalf 0, hits 0
part2: cand 561052, auxkill 354538, smooth 560625, firstposs 22, tangent 1, exact 1, firstver 1, prhalf 0, hits 0
part3: cand 557184, auxkill 353764, smooth 556654, firstposs 19, tangent 1, exact 1, firstver 1, prhalf 0, hits 0
part4: cand 557304, auxkill 353112, smooth 556930, firstposs 19, tangent 0, exact 0, firstver 0, prhalf 0, hits 0
part5: cand 561654, auxkill 354349, smooth 561220, firstposs 28, tangent 1, exact 1, firstver 1, prhalf 0, hits 0
```

Compared with the mod-49 smooth-strata run, the mod-343 gate reduced:

```text
pair budget        8081977 -> 5495485
aux-filtered cand  4951090 -> 3365620
firstposs          177     -> 121
exact tests        5       -> 4
hits               0       -> 0
```

Conclusion: the sharper local gate is working and removes about one third of
the height-100 work, but it does not change the qualitative outcome.  On the
six smooth `Q_7` strata, the exact first-half checks remain rare and every
verified rational first half still fails the final `P_R` halving step.

## D. Aux-prime scanner, per-stratum counters, and `37,41` rerun

New script/logs:

```text
code/agent_m18_416_aux_prime_sieve_scan.m
data/agent_m18_416_aux_prime_scan_h30.log
data/agent_m18_416_mod343_aux37_41_h100_part0.log
data/agent_m18_416_mod343_aux37_41_h100_part1.log
data/agent_m18_416_mod343_aux37_41_h100_part2.log
data/agent_m18_416_mod343_aux37_41_h100_part3.log
data/agent_m18_416_mod343_aux37_41_h100_part4.log
data/agent_m18_416_mod343_aux37_41_h100_part5.log
```

I added a standalone aux-prime sieve scanner.  It reads the mod-343 gate,
builds the height-bounded rational `(R,w)` buckets, applies the current base
aux primes, then estimates the extra cut from each candidate prime and each
candidate-prime pair.  It deliberately stops before the expensive exact
halving chain.

Command used:

```text
magma -b height:=30 base_primes:="11,13" \
  candidate_primes:="17,19,23,29,31,37,41" \
  gate_mod:=343 \
  gate_file:=data/agent_m18_416_smooth_strata_mod343_gate.txt \
  code/agent_m18_416_aux_prime_sieve_scan.m \
  > data/agent_m18_416_aux_prime_scan_h30.log
```

Finite-field profiles:

```text
p=11: good_open 28, first_Tx_half 20, PR_half 6, target416 4
p=13: good_open 56, first_Tx_half 28, PR_half 6, target416 4
p=17: good_open 136, first_Tx_half 64, PR_half 16, target416 12
p=19: good_open 164, first_Tx_half 104, PR_half 22, target416 16
p=23: good_open 288, first_Tx_half 184, PR_half 42, target416 32
p=29: good_open 528, first_Tx_half 256, PR_half 60, target416 32
p=31: good_open 640, first_Tx_half 376, PR_half 86, target416 56
p=37: good_open 968, first_Tx_half 516, PR_half 104, target416 64
p=41: good_open 1240, first_Tx_half 644, PR_half 172, target416 112
```

Height-30 empirical scan, after the mod-343 gate and base primes `11,13`:

```text
gate budget 45789
raw checked 45634
base candidates 29112
base aux-killed 16522
base first-cover possible 244
```

Per-stratum base counts at height 30:

```text
<3,3>: base_cand 3865, base_firstposs 55
<3,4>: base_cand 3865, base_firstposs 62
<4,0>: base_cand 4240, base_firstposs 1
<5,0>: base_cand 8550, base_firstposs 6
<5,2>: base_cand 4296, base_firstposs 64
<5,5>: base_cand 4296, base_firstposs 56
```

Best single-prime cuts at height 30:

```text
p=41: cand 8480, firstposs 225
p=37: cand 8636, firstposs 229
p=31: cand 9772, firstposs 224
```

Best two-prime cuts at height 30:

```text
p=37,41: cand 2600, firstposs 223
p=31,41: cand 2766, firstposs 224
p=31,37: cand 2958, firstposs 222
```

Interpretation: the larger aux primes cut the raw candidate count very well,
but most first-cover possibilities pass these extra primes.  This means the
extra primes are still useful for speed, but the remaining obstruction is
closer to the rational first-cover/tangent stage than to the coarse finite
good-reduction target sieve.

I also patched `agent_m18_416_live_stratum_search.m` to print a
`P7_STRATUM_COUNTS` block.  For each stratum, the block records the full gate
residue/budget counts, the residue/budget counts processed in the current
part, and the stage counters:

```text
auxkill, cand, smooth, firstposs, tangent, exact, firstver,
prfail, prhalf, hits
```

Smoke tests:

```text
height=8, gate_file=<5,0> test gate, aux_primes=11,13:
  DONE cand=40 auxkill=50 smooth=40 firstposs=0 ...
  P7_STRATUM_COUNTS prints one <5,0> row.

height=8, inline mod49 smooth_strata_only=true, aux_primes=11,13:
  DONE cand=188 auxkill=228 smooth=168 firstposs=1 ...
  P7_STRATUM_COUNTS prints all six smooth rows.

height=8, full mod343 gate, aux_primes=11,13,37,41:
  DONE cand=18 auxkill=246 smooth=6 firstposs=0 ...
```

I then ran height 100 with the best scanner pair:

```text
magma -b height:=100 aux_primes:="11,13,37,41" gate_mod:=343 \
  gate_file:=data/agent_m18_416_smooth_strata_mod343_gate.txt \
  NParts:=6 Part:=i progress_interval:=50000 \
  code/agent_m18_416_live_stratum_search.m \
  > data/agent_m18_416_mod343_aux37_41_h100_part{i}.log
```

Completed aggregate:

```text
height                         100
parameters                     12175
mod343 smooth-strata gate      5692 residues
mod343-gated pair budget       5495485
candidates after aux filters   442068
aux-killed                     5051800
smooth curves                  439560
first-cover possible           13
tangent bases                  0
exact tests                    0
first halves verified          0
P_R halves                     0
hits                           0
```

No `PR_HALF` or `HIT_416` lines occurred in any of the six `37,41` logs.

Per-part summaries:

```text
part0: cand 74192, auxkill 846736, smooth 73739, firstposs 3, tangent 0, exact 0, firstver 0, prhalf 0, hits 0
part1: cand 74013, auxkill 845970, smooth 73723, firstposs 2, tangent 0, exact 0, firstver 0, prhalf 0, hits 0
part2: cand 73489, auxkill 842101, smooth 73062, firstposs 1, tangent 0, exact 0, firstver 0, prhalf 0, hits 0
part3: cand 73518, auxkill 837430, smooth 72988, firstposs 3, tangent 0, exact 0, firstver 0, prhalf 0, hits 0
part4: cand 73353, auxkill 837063, smooth 72979, firstposs 2, tangent 0, exact 0, firstver 0, prhalf 0, hits 0
part5: cand 73503, auxkill 842500, smooth 73069, firstposs 2, tangent 0, exact 0, firstver 0, prhalf 0, hits 0
```

Aggregated per-stratum counts for the `11,13,37,41` height-100 run:

```text
<3,3>: cand 61514, smooth 60873, firstposs 0, tangent 0, exact 0, hits 0
<3,4>: cand 61514, smooth 60873, firstposs 4, tangent 0, exact 0, hits 0
<4,0>: cand 63748, smooth 63748, firstposs 0, tangent 0, exact 0, hits 0
<5,0>: cand 126740, smooth 126740, firstposs 4, tangent 0, exact 0, hits 0
<5,2>: cand 64276, smooth 63663, firstposs 5, tangent 0, exact 0, hits 0
<5,5>: cand 64276, smooth 63663, firstposs 0, tangent 0, exact 0, hits 0
```

Compared with the previous mod-343 run using only `11,13`:

```text
candidates        3365620 -> 442068
aux-killed        2128248 -> 5051800
firstposs         121     -> 13
tangent/exact     4       -> 0
hits              0       -> 0
```

Conclusion: `37,41` is a good default add-on for further height pushes.  At
height 100 it made the run about 7.6 times smaller after aux filtering and
removed all exact-halving work.  The remaining rational first-cover points
are concentrated in `<3,4>`, `<5,0>`, and `<5,2>`, and they all fail before
the tangent-candidate stage.

## E. Tracing first-cover survivors and mod-l^2 tangent diagnostics

New code/logs/data:

```text
code/agent_m18_416_firstposs_l2_diag.m
data/agent_m18_416_trace_firstposs_h8_test.txt
data/agent_m18_416_trace_firstposs_aux37_41_h100_part0.log
data/agent_m18_416_trace_firstposs_aux37_41_h100_part1.log
data/agent_m18_416_trace_firstposs_aux37_41_h100_part2.log
data/agent_m18_416_trace_firstposs_aux37_41_h100_part3.log
data/agent_m18_416_trace_firstposs_aux37_41_h100_part4.log
data/agent_m18_416_trace_firstposs_aux37_41_h100_part5.log
data/agent_m18_416_trace_firstposs_aux37_41_h100_part0.txt
data/agent_m18_416_trace_firstposs_aux37_41_h100_part1.txt
data/agent_m18_416_trace_firstposs_aux37_41_h100_part2.txt
data/agent_m18_416_trace_firstposs_aux37_41_h100_part3.txt
data/agent_m18_416_trace_firstposs_aux37_41_h100_part4.txt
data/agent_m18_416_trace_firstposs_aux37_41_h100_part5.txt
data/agent_m18_416_firstposs_l2_diag_37_41_h100.log
data/agent_m18_416_firstposs_l2_diag_17_41_h100.log
data/agent_m18_416_firstposs_l2_diag_43_71_h100.log
```

I patched `agent_m18_416_live_stratum_search.m` with:

```text
trace_firstposs:=true
trace_file:=...
```

When enabled, every rational pair that passes `FirstCoverPossible` is printed
as a human-readable `FIRSTPOSS` line and, if `trace_file` is supplied, written
to a compact table:

```text
R W st0 st1 gate_R gate_W plus minus tangent
```

The height-8 smoke test with `11,13,37,41` had no first-cover survivors, so
the trace file contained only the header.  I then reran the height-100
`11,13,37,41` search with tracing:

```text
magma -b height:=100 aux_primes:="11,13,37,41" gate_mod:=343 \
  gate_file:=data/agent_m18_416_smooth_strata_mod343_gate.txt \
  trace_firstposs:=true \
  trace_file:=data/agent_m18_416_trace_firstposs_aux37_41_h100_part{i}.txt \
  NParts:=6 Part:=i progress_interval:=50000 \
  code/agent_m18_416_live_stratum_search.m \
  > data/agent_m18_416_trace_firstposs_aux37_41_h100_part{i}.log
```

It reproduces the previous aggregate:

```text
candidates after aux filters   442068
first-cover possible           13
tangent bases                  0
exact tests                    0
hits                           0
```

The 13 traced first-cover points are:

```text
1.  R=11/26,    w=-35/76, t=-4699/55432,  stratum <5,0>, minus cover
2.  R=-71/11,   w=-41/11, t=1539/572,     stratum <5,2>, plus cover
3.  R=-2/43,    w=-17/58, t=531/9245,     stratum <5,2>, plus cover
4.  R=-27/17,   w=-63/43, t=34371/61268,  stratum <5,0>, plus cover
5.  R=-47/10,   w=-66/29, t=1599/500,     stratum <3,4>, plus cover
6.  R=-48/89,   w=2/39,   t=-2834/293077, stratum <3,4>, plus cover
7.  R=-43/2,    w=-58/17, t=531/20,       stratum <3,4>, plus cover
8.  R=-16/25,   w=-2/3,   t=598/3125,     stratum <3,4>, plus cover
9.  R=82/43,    w=41/24,  t=-9102/31433,  stratum <5,2>, plus cover
10. R=-13/31,   w=5/13,   t=85/961,       stratum <5,2>, plus cover
11. R=-25/16,   w=-35/26, t=6075/7808,    stratum <5,0>, plus cover
12. R=-80/19,   w=50/11,  t=4590/4693,    stratum <5,2>, plus cover
13. R=-18/37,   w=-63/82, t=21960/39701,  stratum <5,0>, plus cover
```

All 13 have `tangent=0` over `Q`.

The focused mod-l^2 diagnostic reads the compact trace tables and checks the
tangent-candidate congruence used by `TangentCandidates`.  For each traced
pair and prime `ell`, it counts square-compatible tangent solutions modulo
`ell` and `ell^2`.  A `DIES_MOD_L2` verdict would mean the tangent data
exists mod `ell` but not mod `ell^2`.

First diagnostic:

```text
magma -b primes:="37,41" trace_files:="..." \
  code/agent_m18_416_firstposs_l2_diag.m \
  > data/agent_m18_416_firstposs_l2_diag_37_41_h100.log
```

Summary:

```text
ell=37: SURVIVES_MOD_L2 2, boundary 10, denom 1
ell=41: SURVIVES_MOD_L2 2, boundary 10, denom 1
```

So the stronger `37,41` aux primes cut the raw search very well, but they are
not good primes for explaining these 13 first-cover survivors: most of the
survivors reduce to boundary or denominator-bad points at `37` and `41`, and
the few affine-open cases survive the mod-l^2 tangent diagnostic.

I then ran the same diagnostic over the scanner primes:

```text
primes = 17,19,23,29,31,37,41
```

Summary:

```text
ell=17: NO_TANGENT_MOD_L 2, SURVIVES_MOD_L2 2, boundary 7, denom 2
ell=19: NO_TANGENT_MOD_L 3, boundary 8, denom 2
ell=23: NO_TANGENT_MOD_L 6, SURVIVES_MOD_L2 5, boundary 2
ell=29: NO_TANGENT_MOD_L 7, SURVIVES_MOD_L2 2, boundary 2, denom 2
ell=31: NO_TANGENT_MOD_L 6, SURVIVES_MOD_L2 4, boundary 2, denom 1
ell=37: SURVIVES_MOD_L2 2, boundary 10, denom 1
ell=41: SURVIVES_MOD_L2 2, boundary 10, denom 1
```

No `DIES_MOD_L2` verdict occurs.  Among these primes, 12 of the 13 traced
points have a tangent obstruction already modulo `ell`; the lone exception is
the minus-cover point `(R,w)=(11/26,-35/76)`.

I probed a few more primes:

```text
primes = 43,47,53,59,61,67,71
```

The outlier is killed already modulo `47`:

```text
idx=1, ell=47: modl_roots=4, modl_tangent=0,
               modl2_roots=4, modl2_tangent=0,
               verdict=NO_TANGENT_MOD_L.
```

Again, no `DIES_MOD_L2` verdict occurs in this range.

Conclusion: going one level deeper to mod `ell^2` did not reveal a new
obstruction for the 13 height-100 survivors.  The useful local obstruction is
already visible modulo `ell` in the tangent-candidate congruence, provided one
chooses primes where the point is in the affine good-open chart.  The sample
is covered by mod-`ell` tangent obstructions at small primes; for example,
the set `{17,23,29,47}` covers all 13 traced points.

## F. Promoting the tangent-congruence diagnostic into the live search

New logs/data:

```text
data/agent_m18_416_tangent_sieve_h100_part0.log
data/agent_m18_416_tangent_sieve_h100_part1.log
data/agent_m18_416_tangent_sieve_h100_part2.log
data/agent_m18_416_tangent_sieve_h100_part3.log
data/agent_m18_416_tangent_sieve_h100_part4.log
data/agent_m18_416_tangent_sieve_h100_part5.log
data/agent_m18_416_tangent_sieve_h100_part0.txt
data/agent_m18_416_tangent_sieve_h100_part1.txt
data/agent_m18_416_tangent_sieve_h100_part2.txt
data/agent_m18_416_tangent_sieve_h100_part3.txt
data/agent_m18_416_tangent_sieve_h100_part4.txt
data/agent_m18_416_tangent_sieve_h100_part5.txt
```

I patched `agent_m18_416_live_stratum_search.m` with a new optional stage:

```text
tangent_sieve_primes:="17,23,29,47"
```

This runs only after a rational pair passes `FirstCoverPossible`.  For each
listed prime `ell`, the sieve:

1. skips the prime if `(R,w)` has denominator divisible by `ell`;
2. skips the prime if `(R,w)` is on the finite boundary modulo `ell`;
3. otherwise counts square-compatible solutions to the same tangent
   congruence used by `TangentCandidates`, but only modulo `ell`;
4. kills the rational pair as soon as one good-open prime has no compatible
   tangent modulo `ell`.

The default behavior is unchanged because `tangent_sieve_primes` defaults to
the empty list.  I also added counters:

```text
tansievekill
```

globally, in progress/DONE lines, and per `P7_STRATUM_COUNTS` row.

Smoke test without the new sieve:

```text
magma -b height:=8 aux_primes:="11,13,37,41" gate_mod:=343 \
  gate_file:=data/agent_m18_416_smooth_strata_mod343_gate.txt \
  code/agent_m18_416_live_stratum_search.m
```

Result:

```text
DONE cand=18 auxkill=246 smooth=6 firstposs=0 tansievekill=0 tangent=0 exact=0 firstver=0 prhalf=0 hits=0
```

Then I reran height 100 with the tangent sieve and tracing enabled:

```text
magma -b height:=100 aux_primes:="11,13,37,41" \
  tangent_sieve_primes:="17,23,29,47" \
  gate_mod:=343 \
  gate_file:=data/agent_m18_416_smooth_strata_mod343_gate.txt \
  trace_firstposs:=true \
  trace_file:=data/agent_m18_416_tangent_sieve_h100_part{i}.txt \
  NParts:=6 Part:=i progress_interval:=50000 \
  code/agent_m18_416_live_stratum_search.m \
  > data/agent_m18_416_tangent_sieve_h100_part{i}.log
```

Completed aggregate:

```text
height                         100
parameters                     12175
mod343 smooth-strata gate      5692 residues
mod343-gated pair budget       5495485
candidates after aux filters   442068
aux-killed                     5051800
smooth curves                  439560
first-cover possible           13
tangent-sieve killed           13
tangent bases                  0
exact tests                    0
first halves verified          0
P_R halves                     0
hits                           0
```

No `PR_HALF` or `HIT_416` lines occurred.

Per-part summaries:

```text
part0: cand 74192, auxkill 846736, smooth 73739, firstposs 3, tansievekill 3, tangent 0, exact 0, hits 0
part1: cand 74013, auxkill 845970, smooth 73723, firstposs 2, tansievekill 2, tangent 0, exact 0, hits 0
part2: cand 73489, auxkill 842101, smooth 73062, firstposs 1, tansievekill 1, tangent 0, exact 0, hits 0
part3: cand 73518, auxkill 837430, smooth 72988, firstposs 3, tansievekill 3, tangent 0, exact 0, hits 0
part4: cand 73353, auxkill 837063, smooth 72979, firstposs 2, tansievekill 2, tangent 0, exact 0, hits 0
part5: cand 73503, auxkill 842500, smooth 73069, firstposs 2, tansievekill 2, tangent 0, exact 0, hits 0
```

Kill-prime distribution:

```text
ell=17: 2
ell=23: 6
ell=29: 4
ell=47: 1
```

Aggregated per-stratum counts:

```text
<3,3>: cand 61514,  smooth 60873,  firstposs 0, tansievekill 0, tangent 0, exact 0, hits 0
<3,4>: cand 61514,  smooth 60873,  firstposs 4, tansievekill 4, tangent 0, exact 0, hits 0
<4,0>: cand 63748,  smooth 63748,  firstposs 0, tansievekill 0, tangent 0, exact 0, hits 0
<5,0>: cand 126740, smooth 126740, firstposs 4, tansievekill 4, tangent 0, exact 0, hits 0
<5,2>: cand 64276,  smooth 63663,  firstposs 5, tansievekill 5, tangent 0, exact 0, hits 0
<5,5>: cand 64276,  smooth 63663,  firstposs 0, tansievekill 0, tangent 0, exact 0, hits 0
```

Conclusion: the promoted tangent-congruence sieve exactly reproduces the
diagnostic coverage of the 13 height-100 first-cover survivors and prevents
them from reaching the rational tangent/exact-halving stage.  For future
height pushes, the recommended search command is now:

```text
magma -b height:=H aux_primes:="11,13,37,41" \
  tangent_sieve_primes:="17,23,29,47" \
  gate_mod:=343 \
  gate_file:=data/agent_m18_416_smooth_strata_mod343_gate.txt \
  NParts:=N Part:=i \
  code/agent_m18_416_live_stratum_search.m
```

## G. Height-150 run with the promoted tangent sieve

I pushed the same stack to height 150:

```text
magma -b height:=150 aux_primes:="11,13,37,41" \
  tangent_sieve_primes:="17,23,29,47" \
  gate_mod:=343 \
  gate_file:=data/agent_m18_416_smooth_strata_mod343_gate.txt \
  trace_firstposs:=true \
  trace_file:=data/agent_m18_416_tangent_sieve_h150_part{i}.txt \
  NParts:=6 Part:=i progress_interval:=100000 \
  code/agent_m18_416_live_stratum_search.m \
  > data/agent_m18_416_tangent_sieve_h150_part{i}.log
```

The run completed cleanly in six parts.  No `PR_HALF` or `HIT_416` lines
occurred.

Completed aggregate:

```text
height                         150
parameters                     27431
mod343 smooth-strata gate      5692 residues
mod343-gated pair budget       27889578
candidates after aux filters   2278948
aux-killed                     25606998
smooth curves                  2273318
first-cover possible           18
tangent-sieve killed           17
tangent bases                  0
exact tests                    0
first halves verified          0
P_R halves                     0
hits                           0
```

Per-part summaries:

```text
part0: cand 380451, auxkill 4282953, smooth 379438, firstposs 4, tansievekill 4, tangent 0, exact 0, hits 0
part1: cand 379787, auxkill 4275803, smooth 379127, firstposs 4, tansievekill 4, tangent 0, exact 0, hits 0
part2: cand 380418, auxkill 4272423, smooth 379470, firstposs 1, tansievekill 1, tangent 0, exact 0, hits 0
part3: cand 379499, auxkill 4258491, smooth 378275, firstposs 4, tansievekill 4, tangent 0, exact 0, hits 0
part4: cand 378605, auxkill 4248123, smooth 377751, firstposs 3, tansievekill 2, tangent 0, exact 0, hits 0
part5: cand 380188, auxkill 4269205, smooth 379257, firstposs 2, tansievekill 2, tangent 0, exact 0, hits 0
```

Kill-prime distribution among the 18 first-cover positives:

```text
ell=17: 3
ell=23: 6
ell=29: 5
ell=47: 3
not killed by these primes: 1
```

The single first-cover positive not killed by the four-prime tangent sieve was:

```text
R=121/13, w=119/15, t=-4692/2197, stratum=<5,0>, mod343=<194,168>,
plus=true, minus=false, tangent=0, tansievekill=false
```

Thus it survives the selected modulo-`ell` tangent sieve primes, but the
rational tangent-candidate computation is already empty.  It still does not
reach exact halving.

Full first-cover-positive trace:

```text
part0 R=11/26    w=-35/76   stratum=<5,0> mod343=<40,329>  tangent=0 tansievekill=true  kill_ell=47
part0 R=-71/11   w=-41/11   stratum=<5,2> mod343=<243,121> tangent=0 tansievekill=true  kill_ell=29
part0 R=50/117   w=65/102   stratum=<3,4> mod343=<276,4>   tangent=0 tansievekill=true  kill_ell=47
part0 R=-2/43    w=-17/58   stratum=<5,2> mod343=<327,254> tangent=0 tansievekill=true  kill_ell=23
part1 R=-16/57   w=-70/111  stratum=<5,0> mod343=<96,21>   tangent=0 tansievekill=true  kill_ell=29
part1 R=-27/17   w=-63/43   stratum=<5,0> mod343=<180,182> tangent=0 tansievekill=true  kill_ell=29
part1 R=-29/30   w=-116/125 stratum=<3,4> mod343=<262,32>  tangent=0 tansievekill=true  kill_ell=47
part1 R=-47/10   w=-66/29   stratum=<3,4> mod343=<304,116> tangent=0 tansievekill=true  kill_ell=23
part2 R=-48/89   w=2/39     stratum=<3,4> mod343=<38,88>   tangent=0 tansievekill=true  kill_ell=29
part3 R=-74/135  w=21/40    stratum=<5,0> mod343=<96,112>  tangent=0 tansievekill=true  kill_ell=17
part3 R=-43/2    w=-58/17   stratum=<3,4> mod343=<150,158> tangent=0 tansievekill=true  kill_ell=23
part3 R=-16/25   w=-2/3     stratum=<3,4> mod343=<164,228> tangent=0 tansievekill=true  kill_ell=29
part3 R=82/43    w=41/24    stratum=<5,2> mod343=<313,16>  tangent=0 tansievekill=true  kill_ell=23
part4 R=121/13   w=119/15   stratum=<5,0> mod343=<194,168> tangent=0 tansievekill=false kill_ell=0
part4 R=-13/31   w=5/13     stratum=<5,2> mod343=<243,317> tangent=0 tansievekill=true  kill_ell=23
part4 R=-25/16   w=-35/26   stratum=<5,0> mod343=<320,91>  tangent=0 tansievekill=true  kill_ell=17
part5 R=-80/19   w=50/11    stratum=<5,2> mod343=<68,254>  tangent=0 tansievekill=true  kill_ell=23
part5 R=-18/37   w=-63/82   stratum=<5,0> mod343=<222,154> tangent=0 tansievekill=true  kill_ell=17
```

Aggregated per-stratum counts:

```text
<3,3>: cand 319656, smooth 318225, firstposs 0, tansievekill 0, tangent 0, exact 0, hits 0
<3,4>: cand 319656, smooth 318225, firstposs 6, tansievekill 6, tangent 0, exact 0, hits 0
<4,0>: cand 314820, smooth 314820, firstposs 0, tansievekill 0, tangent 0, exact 0, hits 0
<5,0>: cand 650150, smooth 650150, firstposs 7, tansievekill 6, tangent 0, exact 0, hits 0
<5,2>: cand 337333, smooth 335949, firstposs 5, tansievekill 5, tangent 0, exact 0, hits 0
<5,5>: cand 337333, smooth 335949, firstposs 0, tansievekill 0, tangent 0, exact 0, hits 0
```

Conclusion: through height 150, the exact halving stage remains completely
quiet.  The promoted tangent sieve removes 17 of the 18 first-cover positives;
the remaining point is already killed by the rational tangent-candidate test.
So the current computational evidence still points away from a hidden
`[4,16]` example in the searched region.

## Suggested next steps

1. For symbolic work, ignore the repeated degree-5 factor in
   `Res_b(E1core,E0core)`: it is exactly the excluded `d4=0` boundary.  Either
   study the residual degree-52 projection, or go back to the localized
   residual-`G` square-structure analysis at `<0,0>`.
2. For search work, use the mod-343 gate and aux primes `11,13,37,41` for any
   further height increases, together with the tangent-congruence sieve:
   `gate_mod:=343`,
   `gate_file:=data/agent_m18_416_smooth_strata_mod343_gate.txt`,
   `aux_primes:="11,13,37,41"`,
   `tangent_sieve_primes:="17,23,29,47"`.
3. A height-200 push is now the natural next computational test.  The
   tangent-congruence stage only runs on first-cover survivors, so it adds
   negligible overhead while keeping the exact chain quiet.
4. If future height pushes produce several `tansievekill=false` survivors,
   first add a few more small tangent-sieve primes.  Revisit mod-`ell^2` only
   if points survive many good modulo-`ell` tangent tests with nonempty
   rational tangent candidates.
5. Study the residual degree-52 projection in parallel, especially the
   first-cover/tangent locus, since the computational obstruction has moved
   from coarse aux-prime filtering to the tangent-candidate stage.
