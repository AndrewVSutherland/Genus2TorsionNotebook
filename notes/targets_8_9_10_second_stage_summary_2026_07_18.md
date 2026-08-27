# Targets #8--#10: second-stage comparison

Date: 2026-07-18

## Executive result

No new geometrically simple realization was found.  The three parallel
attacks nevertheless produced substantially sharper search boundaries:

| rank | target | strongest completed result | exact hit |
|---:|---|---|---:|
| 8 | `[2,2,4,8]` | all 26,653 `tor2244` rows, all 799,590 ordered Weierstrass charts: zero full-cover witnesses | 0 |
| 9 | `[2,2,2,24]` | exhaustive finite-group target masks at 29,31,37,41 and boundary-aware q-square search through height 187 | 0 |
| 10 | `[2,2,4,12]` | exhaustive fixed-partition mod-`11^2` boundary classification and q-square CRT search through height 160 | 0 |

The negative statements have different scopes.  The #8 chart sweep is
exhaustive for a finite bank.  The #9 and #10 height searches are exhaustive
only on their stated two-dimensional q-square slices and boxes, not on the
full three-dimensional order-12 family.

## #8: `[2,2,4,8]`

### Direct cover

The exact full-cover equations for

```text
M(2,2,4,8) -> A(2,2,4,4)
```

were applied to every row of `tor2244.txt`.  Testing only the displayed
`0,infinity` chart gives zero witnesses but is not complete.  The completed
driver tests all 30 ordered choices of two of the six Weierstrass points to
move to `0,infinity`, followed by all finite-root permutations and cover sign
choices:

```text
rows                              26,653
ordered Weierstrass charts       799,590
square-normalized charts         160,018
full-cover rows/charts/witnesses       0
```

The known HPL split examples are positive controls: each has six square
normalizations, four full charts, and 32 full witnesses.

### Richelot source search

The first 100 bank rows and the corrected chunks 101--1300 give:

```text
rows 1--100:
  reverse edges 1500; one-split 600; exact tests 558
  exact [2,4,8] sources 2 (row 51, edges 13/14)

rows 101--1300:
  distinct bases 1075; duplicate bases 125
  reverse edges 16125; Jacobian codomains 16093
  one-split/gate passes 6410; unique exact tests 6348
  exact [2,4,8] sources 4
    row 136, edges 13/14
    row 1140, edges 13/14
```

All six exact `[2,4,8]` sources are proven geometrically nonsimple by a
rational Richelot quotient of type `SetCart`, explicitly a product of two
elliptic curves.  Their full-Weierstrass forward neighbors were all tested
exactly and none has torsion `[2,2,4,8]`.  The independent row-51 audit also
finds two rational degree-2 elliptic subcovers of the base curve.

Artifacts:

- `notes/target_2248_second_stage_2026_07_18.md`
- `notes/target_2248_source51_audit.md`
- `code/target_2248_tor2244_all_charts.m`
- `results/target_2248_tor2244_all_charts.log`
- `results/target_2248_tor2244_rows101_500.log`
- `results/target_2248_tor2244_rows501_900.log`
- `results/target_2248_tor2244_rows901_1300.log`

## #9: `[2,2,2,24]`

On the q-square slice of the direct `[2,2,2,12]` construction, exhaustive
finite calculations at `p=29,31,37,41` checked the exact marked order-12
class and its divisibility in the complete finite Jacobian group.  Prime 31
is decisive on the smooth open:

```text
smooth/open marked presentations       576
exact order 12                         576
divisible by 2 in J(F_31)                0
```

This includes halves outside the affine `(M,N)` chart.  Therefore any
rational target on this slice must reduce at 31 to a zero/pole, sheet, or
branch-collision/bad-reduction boundary.

The conservative boundary-aware height-187 scan gave:

```text
positive rational values             21,407
ordered pairs                    95,427,692
four-prime modular survivors        982,078
exact nontrivial double-square survivors 0
```

The known `[2,2,2,12]` record is an exact-square height-187 positive control
and is rejected specifically because it is smooth and nontarget at 31 and
37.

Artifacts:

- `notes/target_22224_second_stage_2026_07_18.md`
- `code/target_22224_qsquare_crt_sieve.m`
- `code/target_22224_qsquare_crt_sieve.cpp`
- `code/target_22224_qsquare_crt_exact.m`
- `results/target_22224_qsquare_crt_sieve.log`
- `results/target_22224_qsquare_crt_sieve_h187.log`

## #10: `[2,2,4,12]`

For one fixed second-half partition, every projective mod-11 cover base and
every normalized lift to `11^2` was classified by the four cross-difference
squareclasses:

```text
mod-11 boundary cover bases            536
normalized mod-11^2 lifts          713,416
killed                              449,432
resolved compatible                240,620
deep unresolved                     23,364
```

Thus 11-adic boundary structure is useful but is not a local obstruction.
On the genuine two-dimensional q-square order-12 slice, the combined
order-3/second-half cover is branch-boundary-only through `p=41`, with
smooth positive controls at 43, 47, and 53.

The height-160 CRT search gave:

```text
ordered pairs                    50,744,445
finite-mask survivors                    32
exact double-square points               31
smooth double-square points               0
exact second-half points                  0
```

All 31 exact points lie on the elementary collision component
`(A,B,C)=(r,1,r^-1)` and fail the second-half condition.

Artifacts:

- `notes/target_22412_second_stage_2026_07_18.md`
- `code/target_22412_boundary_crt.py`
- `code/target_22412_boundary_crt_verify.m`
- `results/target_22412_boundary_crt_h160.log`
- `results/target_22412_boundary_crt_verify_h160.log`

## Recommended order of continuation among these three

1. **#9:** resolve the six explicit `p=31` collision divisors and toric
   valuation chambers.  The complete finite-group obstruction has already
   reduced every possible q-square target to these bad-reduction strata.
2. **#10:** leave the q-square slice and feed the compatible mod-`11^2`
   cross-difference masks into the full `(u,t,v)` quartic-splitting solver,
   intersected with the mod-13 contact boundary.  This retains one additional
   geometric dimension but is now locally structured.
3. **#8:** do not enlarge the `tor2244` box or the HPL neighborhood blindly.
   Every completed direct-bank chart fails, and all six new exact sources are
   decomposable.  A future attack should use a genuinely different family or
   a symbolic simultaneous-cover fibration, not more of the same bank.
