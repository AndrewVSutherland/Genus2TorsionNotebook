# `[2,2,4,12]`, second stage: boundary normalization and q-square CRT search

Date: 2026-07-18

## Outcome

No rational `[2,2,4,12]` realization was found, but this stage materially
sharpens both the local and global search problem.

- The fixed-partition second-half cover was classified exhaustively modulo
  `11^2`.  Among 713,416 normalized lifts of the 536 mod-11 boundary bases,
  240,620 resolve to a compatible `Q_11` squareclass, 23,364 remain deep,
  and 449,432 are killed.  Thus the mod-11 boundary is highly structured but
  is not a local obstruction.
- The order-3 equations were coupled directly to the second-half cover on the
  genuine two-dimensional `q=(x+t)^2` slice of `A(2,2,2,12)`.  On this slice
  the combined cover is branch-boundary-only at every tested prime through
  41, but has smooth positive controls at 43, 47, and 53.
- An exact CRT-guided rational search of height 160 checked 50,744,445
  symmetry-reduced pairs.  Only 32 passed all finite masks.  Of these, 31
  satisfy the two exact q-square equations, all 31 lie on the elementary
  branch-collision component, and all 31 fail the exact second-half cover.
  There are no smooth candidates.

This computation is independent of the old 26,653-row `tor2244.txt` bank.
It is exhaustive only on the specified height-160 box in the two-dimensional
q-square subfamily, not on the full three-dimensional `A(2,2,2,12)` family.

## 1. The global fixed-partition boundary calculation

On

```text
y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2),
```

fix the partition `[a,b | c,d]`.  Put

```text
Delta_13 = c^2-a^2,   Delta_14 = d^2-a^2,
Delta_23 = c^2-b^2,   Delta_24 = d^2-b^2.
```

The four second-half radicands are pairwise products of these four cross
differences.  Away from zero, all four radicands are squares if and only if
the four `Delta_ij` have the same class in `Q_p^*/Q_p^{*2}`.

This formulation resolves the forced valuations on a boundary lift.  For a
cross difference known modulo `p^2`, the classifier divides by its known
power of `p` and records

```text
(valuation mod 2, Legendre symbol of the unit).
```

If a difference is zero modulo `p^2`, it is marked `deep`; it is not silently
declared square.  A lift is:

- `resolved` if all four determined squareclasses agree;
- `killed` if two determined squareclasses differ;
- `deep` if the determined classes agree but at least one difference still
  vanishes modulo `p^2`.

Projective points are normalized by making their first nonzero mod-`p`
coordinate exactly 1.  Each mod-11 base then has exactly `11^3=1331`
normalized lifts to `11^2`.

### Complete mod-11 and mod-`11^2` counts

```text
projective bases                         1464
ambient open bases                        192
mod-11 second-half cover bases            536
open cover bases                            0
boundary cover bases                      536
normalized mod-11^2 lifts              713416

killed lifts                           449432
deep lifts                              23364
resolved compatible lifts              240620
all compatible or unresolved lifts     263984

mod-11 bases with all lifts killed          72
deep-only bases                            276
resolved-only bases                        180
resolved-and-deep bases                      8
```

The exact signed component marginals are:

| component type | bases per signed component | killed lifts | deep lifts | resolved lifts | base outcomes `(killed, deep, resolved, mixed)` |
|---|---:|---:|---:|---:|---:|
| `E12+/-`, `E34+/-` | 93 | 51,840 | 3,542 | 68,401 | `(0,38,51,4)` |
| `E13+/-`, `E14+/-`, `E23+/-`, `E24+/-` | 77 | 96,489 | 5,478 | 520 | `(19,54,0,4)` |
| each `Zi` | 53 | 29,964 | 5,973 | 34,606 | `(4,23,26,0)` |

These are nonexclusive marginals: a base on an intersection contributes to
each component containing it.  The complete signed-signature table is in the
log.  It shows the expected dichotomy.  Collisions within one side of the
partition (`E12`, `E34`) often resolve immediately because they do not occur
among the cross differences.  A simple cross collision usually forces a
deeper chart.  Seventy-two multi-collision bases are killed at the first
normalized lift stage, while the eight full-collision signatures have both
resolved and deep lifts.

### Scope of this result

This part is global with respect to the standard second-half cover: it makes
no q-square or old-bank assumption.  It is for the fixed partition
`12|34`; the other two partitions are obtained by permuting labels.  It does
not impose the order-3/contact equations.  Since 263,984 mod-`11^2` classes
remain compatible or deep, it is not a `Q_11` nonexistence proof.

## 2. Pullback to the q-square `A(2,2,2,12)` slice

The direct slice from `target_22224_order12_halving.m` is

```text
C = 1/(A*B),
rho^2   = A^2+B^2+C^2-3,
sigma^2 = A^-2+B^-2+C^-2-3,
s=1/(2*rho),
(a,b,c)=s(A,B,C),
d=2s^2 sigma.
```

After dividing the four branch square-root parameters by the common factor
`s`, their projective class is simply

```text
[A, B, 1/(AB), sigma/rho].
```

Thus the order-3/contact condition and the second-half cover can be tested
with only the two rational variables `A,B`.

Let

```text
R = A^2+B^2+(AB)^-2-3 = rho^2,
S = A^-2+B^-2+(AB)^2-3 = sigma^2,
N = A^2 B^2 R = (AB rho)^2.
```

For the representative partition `[A,B | (AB)^-1,sigma/rho]`, the four raw
radicands have squareclasses

```text
[Z*W, X*Z, Y*W, X*Y],
```

where

```text
X=A^2-1,
Y=B^2-1,
Z=A^4 B^2-1,
W=A^2 B^4-1.
```

For example, before discarding known square factors they are

```text
R1 = Z*W/(A^4 B^4),
R2 = X^3*Z/(A^2*N),
R3 = Y^3*W/(B^2*N),
R4 = A^2 B^2 X^3 Y^3/N^2.
```

Their product is a square, so this remains a rank-three, degree-eight cover.
The other two partitions are its images under the `S_3` symmetry permuting
`A,B,C`.

This factorization is a useful next algebraic entry point: off the boundary,
the target condition says that `X,Y,Z,W` all have the same squareclass.

## 3. Exact finite intersection with the 3-primary slice

For each finite field, the driver enumerates all oriented nonzero `(A,B)`,
requires nonzero square values of `R` and `S`, forms the projective branch
tuple above, and tests all three second-half partitions.

| `p` | q-square slice bases | branch-smooth slice bases | bases with a second half | smooth bases with a second half |
|---:|---:|---:|---:|---:|
| 5  | 12  | 0   | 12  | 0 |
| 7  | 24  | 0   | 0   | 0 |
| 11 | 48  | 0   | 24  | 0 |
| 13 | 60  | 0   | 60  | 0 |
| 17 | 84  | 0   | 84  | 0 |
| 19 | 120 | 0   | 48  | 0 |
| 23 | 144 | 0   | 72  | 0 |
| 29 | 276 | 96  | 156 | 0 |
| 31 | 240 | 48  | 96  | 0 |
| 37 | 404 | 192 | 212 | 0 |
| 41 | 468 | 192 | 228 | 0 |
| 43 | 480 | 192 | 192 | 48 |
| 47 | 600 | 288 | 360 | 192 |
| 53 | 804 | 432 | 468 | 144 |

At `p=11`, the 24 cover bases have unsigned collision signatures

```text
E14: 8, E24: 8, E34: 8.
```

At `p=13`, the 60 cover bases have signatures

```text
E14: 16, E24: 16, E34: 16,
E14+E23: 4, E13+E24: 4, E12+E34: 4.
```

This is the promised direct intersection with the 3-primary constraints:
the q-square identity supplies an actual rational order-3 contact class, not
merely the necessary divisibility test `3 | #J(F_p)`.  On this slice the
combined target is forced to the branch boundary at both 11 and 13.

This is still a statement about the affine q-square chart.  Points with
`rho=0`, `sigma=0`, or poles in `A,B` require other p-adic charts and were
retained, not rejected, in the rational residue masks.

## 4. Positive controls

The Magma verifier checks that the computations are not vacuous.

### Global standard-cover controls

```text
p=17, roots [1,4,2,8], partition 12|34:
    second half true, #J(F_17)=576, 3 | #J(F_17).

p=19, roots [1,4,3,7], partition 12|34:
    second half true, #J(F_19)=576, 3 | #J(F_19).
```

These reproduce the earlier good-open global controls.  They do not lie on
the q-square target open.

### Smooth q-square controls

```text
p=43: A=2, B=3, C=36, D^2=15, partition 14|23,
      Order(D3)=3, second half true, #J(F_43)=2304.

p=47: A=2, B=17, C=18, D^2=37, partition 12|34,
      Order(D3)=3, second half true, #J(F_47)=2304.
```

Thus the q-square pullback of the target cover has a genuine smooth open and
the formula for the order-3 class is being tested directly.

## 5. Height-160 CRT-guided rational search

The search uses all positive reduced rationals `n/d` with

```text
1 <= n,d <= 160.
```

There are 15,611 such values.  The `S_3` symmetry is quotiented by requiring

```text
A <= B <= C=1/(AB).
```

This leaves exactly 50,744,445 pairs.

For each prime

```text
11,13,17,19,23,29,31,37,41,43,47,53,
59,61,67,71,73,79,83,89,97,
```

an exact `(A mod p,B mod p)` mask requires both q-square equations and at
least one second-half partition.  Masks are applied in increasing density
order.  If a numerator or denominator is divisible by `p`, or if `R=0`, the
point is retained as a chart-boundary case.  The search therefore does not
discard candidates merely because this affine reduction is undefined.

Result:

```text
ordered pairs                         50,744,445
finite-mask survivors                         32
fail an exact q-square equation                1
exact double-square points                    31
smooth double-square points                    0
branch-collision double-square points         31
exact second-half points                       0
smooth candidates                              0
```

All 31 exact double-square survivors are

```text
(A,B,C)=(r,1,r^-1),   D^2=1,
```

with collision signature `E24`.  This is the elementary rational component
already seen in the earlier height-20 q-square search.  It cannot support
the second half over the reals: for the only collision-compatible partition,
one required radicand is

```text
(1-r^2)(1-r^-2)=-(r-r^-1)^2,
```

which is negative for `r != 1`; `r=1` is more singular.  The exact tests
accordingly reject all 31.

The candidate file contains only its header, so no exact torsion or
root-power call was necessary.  The Magma verifier nevertheless implements
the required policy: any future survivor is reconstructed on the direct
`A(2,2,2,12)` chart, exact `TorsionSubgroup` must be `[2,2,4,12]`, and
geometric simplicity is certified by an irreducible Frobenius quartic whose
Frobenius powers `2,...,12` all retain degree four.

## 6. What this does and does not prove

There are three distinct scopes.

1. **Global fixed-partition standard cover:** the mod-`11^2` component table
   is exhaustive for all projective branch tuples and one second-half
   partition.  It does not impose order 3.
2. **Global `A(2,2,2,12)` family:** not exhausted.  Its direct `(u,t,v)`
   chart has one more dimension than the q-square slice.
3. **Q-square subfamily:** the finite tables are exhaustive over the stated
   finite fields, and the rational search is exhaustive in the stated
   positive height-160 box.  This is a new search, not the 26,653-row
   `A(2,2,4,4)+3` bank test.

In particular, branch-boundary-only behavior through `p=41` is a strong CRT
constraint on this slice, not a global nonexistence theorem.

## 7. Recommended next attack

The q-square factorization suggests a sharper algebraic calculation than a
larger rectangular search.  Put

```text
X=A^2-1, Y=B^2-1, Z=A^4B^2-1, W=A^2B^4-1.
```

Off the boundary, impose

```text
Y/X = r1^2,  Z/X = r2^2,  W/X = r3^2
```

together with the two q-square equations.  Eliminate one of `A,B`, or use
the symmetric variables

```text
x=A+B+C,  y=AB+AC+BC,
rho^2=x^2-2y-3,
sigma^2=y^2-2x-3,
```

plus the splitting condition for
`T^3-xT^2+yT-1`.  The smooth residues at 43 and 47 give Hensel seed fibers
for this calculation.  This should expose elliptic or higher-genus fibers
and separates the elementary real-obstructed component automatically.

In parallel, the full-family route should feed the exact compatible
mod-`11^2` cross-difference masks into the direct quartic-splitting solver
and intersect them with the existing mod-13 contact boundary masks.  That is
the natural way to leave the q-square slice without returning to a blind
three-dimensional height search.

## 8. Reproduction and artifacts

```bash
python3 code/target_22412_boundary_crt.py \
  --height 160 \
  --log results/target_22412_boundary_crt_h160.log \
  --candidates results/target_22412_boundary_crt_candidates_h160.txt

magma -b \
  candidate_file:=results/target_22412_boundary_crt_candidates_h160.txt \
  log_file:=results/target_22412_boundary_crt_verify_h160.log \
  code/target_22412_boundary_crt_verify.m
```

Artifacts:

```text
code/target_22412_boundary_crt.py
code/target_22412_boundary_crt_verify.m
results/target_22412_boundary_crt_h160.log
results/target_22412_boundary_crt_candidates_h160.txt
results/target_22412_boundary_crt_verify_h160.log
notes/target_22412_second_stage_2026_07_18.md
```

