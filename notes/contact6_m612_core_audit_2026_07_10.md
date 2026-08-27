# Audit of the contact-6 `[6,12]` dual split-core search

Date: 2026-07-10

## Conclusion

The reported `core_lifts > 0` but `verified_core = 0` did **not** reveal a
failure of the cubic-contact equations or the exact Jacobian verifier.  It
also did not mean that dozens of smooth contact points represented the
built-in/dependent order-3 direction.

The actual explanation is:

1. almost all formal lifts came from the degree-drop fiber `b=-3`;
2. the remaining square-cover lifts were singular curves; and
3. every smooth independent core point through height 10 fails all three
   necessary dual discriminant covers.

Thus the negative `[6,12]` result through this bound remains valid, but the
old counts and their interpretation were misleading.  There is no open
candidate that the exact verifier incorrectly discarded.

## Positive controls

The shared verifier in
`code/contact6_m612_tb_core_tools.m` was tested independently using the
known simple `[6,6]` source

```text
(a,b)=(133/39,-7/13).
```

All three known cubic-contact presentations pass exactly:

```text
M=9/256,   U=5/12,   v=5/6,
M=841/256, U=-9/4,   v=5/2,
M=169/16,  U=-17/2,  v=-5/4.
```

In each case the marked class has order 6, the reconstructed class has
order 3, and the latter is independent of the marked order-3 direction.
The standalone driver `code/contact6_m612_tb_core_verifier.m` gives the same
three hits.

The new exact dual-class control
`code/contact6_m612_dual_class_exact.m` is also consistent.  On this seed:

```text
Delta = 5972/507,
DB = -29696/1521,       nonsquare,
DC = -19652/507,        nonsquare,
DB*DC = 583585792/771147, nonsquare.
```

The `+` dual twist has exact torsion `[6,6]` and none of `R1,R2,R3` is
divisible by 2.  The `-` twist has exact torsion `[2,2]`, again with no
dual-kernel half.  In particular, the necessary square covers and the exact
class tests agree on the known seed.

## Classification of the old formal lifts

The diagnostic driver
`code/contact6_m612_dual_split_core_audit.m` repeats the slice recovery and
classifies every square-`M` result before the cover filter.

For the current pre-fix driver the classifications are:

| height | formal lifts | on a square cover | degree drop | singular | verified independent core |
|---:|---:|---:|---:|---:|---:|
| 6  | 49  | 45  | 42  | 5  | 2 |
| 10 | 138 | 130 | 122 | 10 | 6 |

At height 6 the cover count consists of all 42 degree-drop points and 3 of
the 5 singular points.  Neither independent core point is on a cover.  At
height 10 it consists of all 122 degree-drop points and 8 of the 10 singular
points.  Again, none of the 6 independent core points is on a cover.

The archived summary reports `136` rather than `138` formal lifts at height
10, and `15875` rather than `16002` checked slices.  Those numbers are
consistent with the earlier revision that omitted the entire special slice
`v=1`:

```text
127 parameter values * 125 values with v != 0,1 = 15875.
```

The current special-slice solver adds two height-10 formal cover points;
both are singular.  It does not alter the open conclusion.

## The boundary bug

For

```text
f=(1+a*x+b*x^2+x^3)^2-(x-1)^6,
```

the coefficient of `x^5` is `2*(b+3)`.  Consequently `b=-3` is a degree
drop and is not a genus-2 source in this chart.

The generic fixed-`(b,v)` solver specialized `b` before saturation.  On the
fiber `b=-3`, multiplying a saturation boundary by `b+3` cannot remove the
fiber: the scalar has already become zero.  In the current generic routine
the factor was absent altogether.  This produced 42 of 49 formal lifts at
height 6 and 122 of 138 at height 10.

There was a second reporting issue: Magma gives dimension `-1` for an empty
ideal, but the driver reported every dimension other than zero as an
exceptional slice.  The special slice `(b,v)=(-1,1)` is empty, not
exceptional.

The fixes in `code/contact6_m612_dual_split_core_search.m` are deliberately
narrow:

- skip `b=-3` explicitly before solving any slice;
- treat dimension `-1` as empty and reserve the negative sentinels `-2,-3`
  for genuine solver failures;
- distinguish `formal_lifts` and `formal_split_cover` from exact
  `verified_core` and `split_cover`; and
- run the exact core verification before interpreting a cover survivor as a
  smooth independent `[3,3]` source.

The corrected runs are:

```text
height 6:
  checked 2116, exceptional 0,
  formal_lifts 7, formal_split_cover 3,
  verified_core 2, split_cover 0, hits 0.

height 10:
  checked 15876, exceptional 0,
  formal_lifts 16, formal_split_cover 8,
  verified_core 6, split_cover 0, hits 0.
```

## The genuine core points

Through height 10 the six verified core presentations are two presentations
on each of the three already known sources:

```text
(a,b)=(-19/9,3/2):
  (M,U,v)=(784/729,-64/27,-2/3),
          (100/729,-14/27, 2/3),
  DB=-800/81, DC=-175/36.

(a,b)=(-43/25,1/8):
  (M,U,v)=(8281/15625,-491/250,-4/5),
          (6561/15625,-329/250, 4/5),
  DB=-1701/625, DC=-3159/1600.

(a,b)=(-15/8,5/9):
  (M,U,v)=(3025/4096,-199/96,-3/4),
          (1225/4096,-101/96, 3/4),
  DB=-2695/576, DC=-245/81.
```

For each source, `DB`, `DC`, and `DB*DC` fail the required square tests.
Moreover all three sources have `DeltaR=0`, so their distinguished Richelot
codomains are degenerate.  They therefore cannot contribute an open simple
dual `[6,12]` target.

## Recommendation

Do not extend this affine height search merely to improve the old raw
counts.  Its corrected height-10 conclusion is clean: there is no smooth
independent core point on any of the three necessary dual square covers.
The useful next work is the projective/bad-reduction boundary analysis, plus
the exact `R1/R2/R3` tests on any genuine boundary family that survives.
