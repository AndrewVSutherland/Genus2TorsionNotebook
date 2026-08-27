# Parallel update, 2026-07-01

Goal: find a geometrically simple genus-2 Jacobian over `Q` with rational
torsion containing either `[4,16]` or `[2,2,16]`.

This note records the parallel follow-up after importing Filip's `torsion_jac`
repo into the shared working copy.  The main conclusion so far is that
the obvious A8/A(2,8)/HPL boxes are not hiding an already-known simple target.
The live routes are boundary/local-chart computations.

## Existing target-looking hits

The workspace audit found two genuine target-invariant hits, but both are
split/non-simple.

- Plain A8 gives `[4,16]` at `rv=pv=tv=-8`, but the sextic is even in `x`, so
  the Jacobian has an elliptic quotient.
- W-split A8 gives `[2,2,16]` at
  `rv=-5/4, tv=1, bv=+-3/4, pv=-7/20`, but certification again shows
  `even_in_x=true` and reducible Frobenius polynomials.

There are simple A8 examples with cyclic `[16]`, but their rational 2-rank is
too small for `[4,16]` or `[2,2,16]`.

## M_1(8,4) route to `[4,16]`

Filip's notes show that the good-open `[4,16]` locus is empty over `F_7`.
Thus rational `[4,16]` points in this chart must reduce to the `p=7` boundary.

### Affine p=7 boundary CRT search

New script:

```text
code/agent_m18_416_search_crt.m
```

This recomputes the cleared `p=7` boundary closure for the `[4,16]` target and
uses those residues as a CRT gate in the affine `(R,w)` chart.  It optionally
adds auxiliary good-open filters at primes such as `11` and `13`.

The p=7 closure counts are:

```text
boundary residues       41
first-closure residues  24
target416 residues      19
```

Height-20 affine results:

```text
p=7 only:
  exact first-halving candidates verified 20
  P_R halves 0
  hits 0

p=7 plus auxiliary primes 11,13:
  exact first-halving candidates verified 5
  P_R halves 0
  hits 0
```

Logs:

```text
data/agent_m18_416_search_h20_p7only.log
data/agent_m18_416_search_h20_aux11_13.log
```

This affine pass does not enumerate parameters whose denominator is divisible
by `7`; those are boundary/infinity chart cases.

### Denominator-p=7 boundary charts

New script:

```text
code/agent_m18_416_den7_boundary_search.m
```

This covers the complementary p=7 denominator charts:

- `Rinf_w`: `R = infinity`, finite `w mod 7`.
- `r_Winf`: finite `R mod 7`, `w = infinity`.
- `Rinf_Winf`: simultaneous infinity.

Only one p=7 denominator chart survives the target closure:

```text
r_Winf:1, i.e. R == 1 mod 7 and w = infinity.
```

Height-30 result on this complement:

```text
total_pairs 1234321
den7_pairs 287592
p7_chart_pass 18906
aux_pass 6052
smooth 5914
exact_tests 5914
Tx_halves 0
PR_halves 0
both_halves 0
torsion_tests 0
hits 0
```

Logs:

```text
data/agent_m18_416_den7_h20.log
data/agent_m18_416_den7_h30.log
```

Conclusion: both the affine p=7 boundary CRT search and the p=7 denominator
complement are negative at these first bounds.  The next `[4,16]` step should
not be another blind height increase; it should be a deeper blowup/Hensel
analysis of the surviving p=7 boundary strata, especially the affine target
closure and `r_Winf:1`.

## M_1(8,2,2) route to `[2,2,16]`

Filip's notes already show the direct rational-open route is strongly
constrained: open halving is obstructed at `p=7` and `p=11`, and the
height-50 rational-open scan through `p<=73` found no survivors.

A new boundary-CRT worker started an optimized integer/CRT enumerator, with
candidate files named `agent_m3222_boundary_*`.  On the coordinator side these
files currently appear as unmaterialized placeholders in the synced working copy: they have logical
sizes but zero allocated blocks and read as empty.  Therefore they should not
yet be treated as reliable candidate data.  A separate exact-test worker was
asked to regenerate or materialize them before reporting counts.

Pending follow-up: trusted exact-test summary for the regenerated
`agent_m3222_boundary_*` candidates.

### Exact-test blocker

The exact-test worker could not run trusted checks on the new
`agent_m3222_boundary_*` candidate files because the synced working copy exposed
them as placeholders.  Several files had nonzero logical sizes, for example
`agent_m3222_boundary_candidates_h80_p73.txt`, but normal reads via `wc`,
`head`, `xxd`, and `file` saw empty content.  The sync client had not materialized them locally, and materialization
could not be triggered from the sandbox here.

Therefore the current reliable exact-test count for these new candidate files
is:

```text
readable candidates 0
exact halving checks 0
new [2,2,16] hits none established
```

This is a technical/materialization blocker, not a mathematical negative.
The next practical step is to materialize the files in the synced working copy
or regenerate them locally into non-placeholder files, then exact-test
`agent_m3222_boundary_candidates_h80_p73.txt` first.

## HPL / M(2,2,4,8) pair fibers

The HPL route has been converted from a height-box search into genus-one
pair-fiber arithmetic.  New scripts/logs live under:

```text
code/agent_m2248_pair_*.m
code/agent_m2248_pair_*.py
data/agent_m2248_pair_*.log
```

Partial pair-rank data for the HPL `F1/F2` fiber:

- The pair fiber is a genus-one quartic.
- The reciprocal quartic check passes.
- A minimal elliptic model was computed.
- Torsion invariants on the minimal model are `[2,2]`.
- The root number is `-1`.
- The small slope/multiple search recovered the HPL full-cover point only.
- The recovered full-cover point is non-simple.

The `q^2=1` projective local check found no local obstruction to depth `6` at
`p=3` or `p=5` on the broad projective closure; many surviving lifts are
boundary-like.  This means the earlier affine good-chart obstruction should be
interpreted carefully: the projective boundary is not empty.

Pending follow-up: final pair-rank worker summary for `F1/F4` and `F2/F4`, or
a recommendation to deprioritize HPL if those remain too slow.

The pair-rank worker did not return a final message before shutdown, but it did
write materialized files:

```text
code/agent_m2248_pair_rank_probe.m
data/agent_m2248_pair_rank_search_h220_m80.log
```

The completed portion covers `F1/F2` at HPL.  It recovers only the HPL
full-cover point in the tested slope/multiple search.  No simple full-cover
point was found.  The uncompleted `F1/F4` and `F2/F4` rank probes should be
treated as still open.

## A(2,2,4,4) / M(2,2,4,8) K3 boundary sweep

The improved B=20000 forced-boundary pipeline was verified and extended.

Parameters:

```text
B = 20000
boundary = 11:N,23:N
local depth = 4
early common A(2,2,4,4) partition prefilter enabled
```

Previous documented improved-filter coverage was `ax=1..100`, all zero.
New chunks:

```text
ax=101..110: rows 0, runtime 73.467s
ax=111..120: rows 0, runtime 91.984s
ax=121..130: rows 0, runtime 67.174s
```

Updated coverage:

```text
ax=1..130, all zero rows.
```

Conclusion: this is a good resumable/distributed sieve, but it is not the
fastest near-term route to a target example.  Continue in independent chunks
if spare compute is available; otherwise prioritize the `[4,16]` p=7 boundary
geometry and the `[2,2,16]` simultaneous `7/11` boundary route.

## Current priority order

1. Deeper p=7 boundary blowups/Hensel charts for the `M_1(8,4)` `[4,16]`
   equations.
2. Regenerate/materialize and exact-test the `m3222` simultaneous `7/11`
   boundary candidates for `[2,2,16]`.
3. Finish the HPL pair-fiber arithmetic enough to decide whether it has any
   non-HPL rational full-cover points.
4. Continue the A2244 B=20000 boundary sweep as a background distributed job.
