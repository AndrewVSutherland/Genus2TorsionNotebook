# U2 at `p=13`: the missing infinity component dies on the full-cover blow-up

Date: 2026-07-18

## Outcome

The U2 infinity component that is absent from the sampled 144-row tangent
bank was analyzed directly.

There are genuine smooth branches of the **cubic-contact equations alone**:
two unit-`L` incidences modulo 13 lift exhaustively through `13^5` and cover
every projective infinity-parameter class.  However, neither branch lies on
the full `A(2,2,2,8)` cover.  Two normalized cover radicands are nonsquare
units modulo 13.  Consequently,

```text
actual U2 infinity full-cover branches over Q_13 = 0.
```

This final statement is a complete local obstruction on the infinity
component and does not depend on the sampled tangent bank.

## Infinity chart and full-cover obstruction

Start from

```text
U2 = [A,B,C,D*T^2]
   = [-1071,-1054,1116,1134*T^2].
```

For `v_13(T)<0`, put `z=1/T`, so `z` is a nonzero element of `13 Z_13`.
After projective scaling the branch tuple is

```text
[A*z^2,B*z^2,C*z^2,D].
```

The four full-halving radicands factor as

```text
R0 = z^6 * A*B*C*D,
R1 = z^6 * A*(A+B)*(A+C)*(D+A*z^2),
R2 = z^6 * B*(B+A)*(B+C)*(D+B*z^2),
R3 = z^6 * C*(C+A)*(C+B)*(D+C*z^2).
```

Because `6*v_13(z)` is even, squarehood is controlled by the four unit
factors.  For every `z` in `13 Z_13` their residues are

```text
[9,7,11,4] modulo 13.
```

The residues `7` and `11` are nonsquares.  Hence `R1` and `R2` are
nonsquares in `Q_13` for every nonzero infinity-chart parameter `z`.
At `z=0` the tuple is the degenerate projective endpoint, not a genus-2
curve.  This exhausts the whole infinity disk.

Equivalently, after dividing a hypothetical square root by `z^3`, its first
digit would have to solve `Y^2=7` or `Y^2=11` over `F_13`, which is
impossible.  The all-zero square-root solution visible on the unblown-up
special fibre is therefore an extraneous boundary component.

## Contact-only lift

For completeness, the corrected cleared cubic-contact equations in
`(z,L,U,v)` were analyzed before imposing the cover.  At `z=0` they have 17
incidences:

```text
smooth unit-L:   (0,4,1,3), (0,9,1,3)
singular unit-L: (0,6,0,0), (0,7,0,0)
L=0:             13 singular incidences
```

The two smooth branches have Jacobian rank three in four variables.  Their
digit lifts were exhaustively enumerated:

| exponent | contact nodes | distinct `z` classes |
|---:|---:|---:|
| 1 | 2 | 1 |
| 2 | 26 | 13 |
| 3 | 338 | 169 |
| 4 | 4394 | 2197 |
| 5 | 57122 | 28561 |

Thus contact alone admits every `z` in the infinity chart, with two smooth
incidences over each class.  The full-cover nonsquare obstruction kills all
of them.  The singular contact branches were not promoted or assumed to
lift; this does not affect the conclusion because the cover obstruction is
independent of `(L,U,v)`.

```text
code/target_22224_U2_infinity_contact_lift.py
results/target_22224_U2_infinity_contact_lift.log
results/target_22224_U2_infinity_contact_lift.tsv
```

The TSV records contact-only classes and explicitly marks nonzero classes as
`full_cover_possible=0`.

## Intersection with the 281 Mordell--Weil survivors

The matching-basis file

```text
results/target_22224_U2_rank2_modular_N10000_fast.tsv
```

was reevaluated at precision 100 using the same saturated generators.  All
281 classes have finite, integral `T`; there are no nonintegral or exact-pole
classes.

The earlier generic deep-filter log reported `280 finite + 1 infinity`.
The exceptional coefficient class is

```text
(m,n,torsion_coset) = (-1,3,1).
```

Its p-adic map evaluation lost precision in a removable presentation and
division raised an exception, which the generic code conservatively labeled
as infinity.  Exact rational evaluation is small and gives

```text
P = (531016/25,-387597636/125),
T = 664293/527543,
v_13(T)=0, T=2 modulo 13.
```

Therefore none of the 281 imported classes even enters the infinity disk.

```text
code/target_22224_U2_infinity_MW_check.m
results/target_22224_U2_infinity_MW_check.log
results/target_22224_U2_infinity_MW_check.tsv
```

## Scope: complete component versus sampled bank

Two logically distinct statements should be retained:

1. The earlier global projective-line scan showed that U2 misses the
   **current sampled 144-key tangent bank** already modulo 13.  That was a
   bank-coverage statement, because the bank omitted the infinity key.
2. The present blow-up proves that the **actual full-cover infinity
   component is empty over `Q_13`**.  This is complete and independent of
   how the tangent bank was sampled.

This note does not claim that every possible denominator-nonunit contact
chart over finite integral `T` has been resolved.  The existing corrected
affine/projected mask has no finite U2 class, and the 281 bounded MW
survivors are all finite, but a claim closing all finite projective contact
charts would require the corresponding chart-completeness argument.  What
is now fully closed is precisely the missing infinity component requested
here.
