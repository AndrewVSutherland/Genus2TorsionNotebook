# A(2,12) -> A(2,24) height-5 low-branch scan

Date: 2026-07-02.

This broadens the split-fiber branch-factor scan from height 4 to height 5,
while avoiding the four already closed height-4 fibers except for one small
sanity run.

New code:

```text
code/agent_A2_24_height5_lowbranch_scan.m
```

Logs:

```text
results/A2_24_height5_full.log
results/A2_24_height5_height4_compare.log
results/A2_24_height5_height4_include_closed_sanity.log
results/A2_24_height5_syntax_smoke.log
```

The scanner keeps split A(2,12) fibers, translates the visible order-12 class
by rational 2-torsion, removes the `s4=0` square-quartic boundary, and records
the saturated affine factor degrees in `M`.  It also checks every saturated
affine factor for rational `M`-roots; if any occur, it lifts to `N`, prints the
candidate line, and runs exact integral-model halving/torsion certification.

## Commands

Syntax smoke:

```text
magma -b Height:=4 MaxChecked:=200 Progress:=0 \
  code/agent_A2_24_height5_lowbranch_scan.m \
  > results/A2_24_height5_syntax_smoke.log 2>&1
```

Height-4 comparison with the four closed fibers skipped:

```text
magma -b Height:=4 Progress:=0 \
  code/agent_A2_24_height5_lowbranch_scan.m \
  > results/A2_24_height5_height4_compare.log 2>&1
```

Full height-5 run, skipping the four closed fibers:

```text
magma -b Height:=5 Progress:=5000 \
  code/agent_A2_24_height5_lowbranch_scan.m \
  > results/A2_24_height5_full.log 2>&1
```

Height-4 sanity run with the four closed fibers included:

```text
magma -b Height:=4 SkipClosedBest:=false Progress:=0 \
  code/agent_A2_24_height5_lowbranch_scan.m \
  > results/A2_24_height5_height4_include_closed_sanity.log 2>&1
```

## Height-5 result

The full height-5 run enumerated all

```text
checked=56316
```

triples in the height-5 box.  Excluding the four closed fibers, it found

```text
split_fibers=84
order12_split_fibers=84
translated_order12_rows=352
errors=0
positive_dim_rows=0
```

The saturated affine factor-degree distribution was:

| saturated affine `M` degrees | translated rows |
|---:|---:|
| `[16]` | 280 |
| `[8, 8]` | 72 |

Equivalently, the minimum affine factor degree distribution was:

| minimum degree | translated rows |
|---:|---:|
| 16 | 280 |
| 8 | 72 |

There were no new low branches and no rational projected branches:

```text
low_rows_le_4=0
degree4_or_less_rows=0
rational_M_root_rows=0
rational_M_root_factors=0
rational_point_rows=0
exact_divisible_rows=0
torsion_cert_rows=0
```

Since no saturated affine factor had a rational `M`-root, no new `N`-lift or
exact A(2,24) torsion certification was triggered.

Four height-5 fibers had residual quartic factor degrees
`[ <1,1>, <1,1>, <2,1> ]` and hence eight rational 2-torsion translations
each:

```text
p=-1/5, z=-1/5, r= 5/4
p=-1/5, z= 1/5, r= 5/4
p= 1/5, z=-1/5, r=-5/4
p= 1/5, z= 1/5, r=-5/4
```

They still contributed only to the `[16]`/`[8,8]` saturated distribution; none
had degree `<=4` or a rational `M`-root.

## Height-4 sanity

With the four known closed fibers skipped, the height-4 comparison gave:

```text
checked=11132
split_fibers=24
order12_split_fibers=24
translated_order12_rows=96
```

and distribution

```text
[16]  : 64 rows
[8,8] : 32 rows
```

With the four closed fibers included, the sanity run gave:

```text
checked=11132
split_fibers=28
order12_split_fibers=28
translated_order12_rows=112
low_rows_le_4=8
rational_M_root_rows=0
```

The only degree-4 rows were exactly the two extra 2-torsion translations on
each of the four previously closed fibers:

```text
(-1/3,-1,  4/3)
(-1/3, 1,  4/3)
( 1/3,-1, -4/3)
( 1/3, 1, -4/3)
```

Each has saturated affine degrees

```text
[ <4,1>, <4,1>, <8,1> ]
```

and `rational_M_roots=[]`, agreeing with
`agent_A2_24_branch_closure.md`.

## Box comparison

Subtracting the height-4 skipped-closed comparison from the height-5
skipped-closed run, the genuinely new height-5 shell contributes

```text
60 split/order-12 fibers
256 translated order-12 rows
```

with distribution

```text
[16]  : 216 rows
[8,8] : 40 rows
```

Thus the full height-5 box, after adding back the known closed rows, has:

```text
[16]    : 288 rows
[8,8]   : 72 rows
[4,4,8] : 8 rows
```

The eight `[4,4,8]` rows are precisely the two extra translations on the four
height-4 best fibers, already closed by the rational-root, boundary, and exact
divisibility checks in `agent_A2_24_branch_closure.md`.

## Verdict

No new split A(2,12) fiber in the height-5 box exposes a saturated affine
halving branch of degree `<=4`, and no saturated affine branch of any degree
has a rational `M`-root.  The four previously closed height-4 best fibers are
the only degree-4 examples in the height-5 box.  Therefore this height-5
expansion produces no A(2,24) candidate and no exact divisibility/torsion
certification target beyond the already closed four fibers.
