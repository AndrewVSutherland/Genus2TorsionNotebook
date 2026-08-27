# Cyclic order `44` from the order-`22` families: bounded scout

## Outcome

No cyclic order-`44` example was found.

The class-specific local cover is nonempty at several good primes, so the
route is not globally disproved.  However, the good open is empty at both
`p=5` and `p=17`.  Consequently every rational half must enter a bad or
projective boundary disk at each of those primes.  The only genuine smooth
member hidden among the displayed boundary points is also killed at both
primes.  The remaining disks specialize to square, non-genus-`2` fibers and
are only formally live to first order.

A conservative rational search through height `1000`, followed by exact
Magma divisibility tests, found no half of the marked order-`22` class.
This makes cyclic `44` a poor immediate example-producing lane unless the
degenerate boundary geometry is resolved first.

All computations were bounded by `8 GiB` virtual memory and at most five
minutes per process.  No elimination or Groebner-basis computation was run.

## The exact marked class

On an even sextic `y^2=f(x)`, let

```text
E = W - infinity_+,       W=(r,0).
```

The divisor calculation in the original order-`22` note gives

```text
2E = infinity_- - infinity_+,
```

so `E` has order `22`.

To remove any convention about the two infinities, set

```text
X = 1/(x-r),       Y = y*X^3,
g(X) = X^6*f(r+1/X).
```

Then `g` is an odd quintic with `g(0)=1`; `W` becomes infinity and
`infinity_+` becomes `(0,1)`.  Thus `E` maps to the negative of

```text
D = (0,1) - infinity,
```

and `E` is divisible by `2` exactly when `D` is.  Magma represents this
class unambiguously as `J![X,1]`.

As a control, on the Flynn specialization `s=2`, `eps=-1`, `t=-36`, the
even-model divisor `(4,0)-infinity_-` has exact order `22` and
`IsDivisibleBy(E,2)` returns `false`.

## A useful cancellation and the equivalence of the two presentations

For source sign `e=+/-1`, the parameter simplifies from the formula in the
original note to

```text
t_e(s) = -[s*(s^2+e*s+1)/(s+e)]^2,
u_e(s) = t_e(s)/4.
```

Thus `s=e` is a removable point, while `s=-e` is the actual pole.

The two source families are not merely similar.  Direct symbolic expansion
gives

```text
G_u(x) = F_(4u)(x-1).
```

The coordinate map is `x_Flynn=x_DS-1`, with `y` unchanged, and sends the
marked root `s^2` to `1+s^2`.  The identity is certified by
`code/order44_from_order22_family_equivalence.m`.  It justifies using only
the Flynn presentation in the larger search without losing any DS curves.

## Class-specific finite-field scout

The script

```text
code/order44_from_order22_finite_scout.m
```

enumerates the honest affine open `s != 0,+/-1`.  For every good model it:

1. moves the marked root to infinity;
2. verifies that `D=J![X,1]` has order `22`;
3. maps `D` to an explicit finite abelian group;
4. tests whether its coordinates lie in `2J(F_p)`.

The four projective base points and every singular or collapsed-order
reduction are kept, never rejected.

The Flynn and DS counts agree exactly:

| `p` | open sheet pairs | good | singular kept | order `22` | divisible by `2` |
|---:|---:|---:|---:|---:|---:|
| 3  | 0  | 0  | 0  | 0  | 0  |
| 5  | 4  | 4  | 0  | 4  | 0  |
| 7  | 8  | 4  | 4  | 4  | 0  |
| 13 | 20 | 12 | 8  | 12 | 4  |
| 17 | 28 | 28 | 0  | 28 | 0  |
| 19 | 32 | 18 | 14 | 18 | 6  |
| 23 | 40 | 28 | 12 | 28 | 4  |
| 29 | 52 | 46 | 6  | 46 | 16 |
| 31 | 56 | 46 | 10 | 46 | 8  |

For example, the Flynn cover has open halves modulo `13` at

```text
(s,e,t) = (5,-1,4), (5,1,9), (8,-1,9), (8,1,4).
```

Therefore the cover is genuinely locally alive.  But at `p=5` and `p=17`
all good open reductions are killed.

## Boundary audit at `p=5` and `p=17`

The script

```text
code/order44_from_order22_boundary_audit.m
```

uses the four regular charts

| chart | base substitution | odd-model scaling |
|---|---|---|
| zero | `s=q` | `X=Z` |
| cancelled point | `s=e+q` | `X=Z` |
| pole | `s=-e+q` | `X=qZ` |
| infinity | `s=1/q` | `X=q^2 Z` |

At the cancelled point, both presentations give the same genuine smooth
curve:

```text
Flynn:             t=-9/4,  r=1
Daowsud-Schmidt:   u=-9/16, r=2.
```

Its marked class has order `22`, but is not divisible by `2` at either
prime:

```text
p=5:  J(F_p) invariants [2,22]
p=17: J(F_p) invariants [2,198]
```

Hence this genuine boundary member has no mod-`p` or mod-`p^2` seed.

The other three scaled special polynomials are

```text
s=0:       (X^2+X+1)^2
s=-e:      (1-X^2)^2
s=infty:   (2X^2+3X+1)^2.
```

They are square, reducible special fibers rather than genus-`2` curves.  If
`h(q,X)` denotes the regular scaled family and `h_5` its leading
coefficient, their first-order data are

| chart | `ord_q(h_5)` | `L_1` in `h=(H_0+qL_1)^2+O(q^2)` |
|---|---:|---|
| zero | 5 | `0` |
| pole, `e=-1` | 2 | `-4X^2+4X` |
| pole, `e=+1` | 2 | `4X^2+4X` |
| infinity | 3 | `0` |

Thus at both `p=5` and `p=17`, setting `q=p` gives a formal solution of
the exact halving identity modulo `p^2`:

```text
ell = half_sign*(H_0+pL_1),       u arbitrary modulo p.
```

These are singular closure seeds, not Hensel seeds on a smooth Jacobian.
They explain why first-order local testing cannot close the order-`44`
route.  A real continuation would require weighted higher-order blowups and
saturation of these three disks.

## Exact halving-cover equations

The non-eliminated cover is generated by

```text
code/order44_from_order22_halving_cover.m
```

After moving the marked root to infinity, write

```text
u(X)   = X^2+aX+b,
ell(X) = mX^2+nX+k,       k=+/-1.
```

The exact condition is

```text
g(X)-ell(X)^2 = g_5*X*u(X)^2.
```

The `X^5` and `X^4` coefficients solve linearly for `a,b`; the remaining
two equations are

```text
K_2(s,m,n)=0,       K_1(s,m,n)=0.
```

In the displayed chart they have `(m,n)` bidegrees `(6,2)` and `(8,2)`.
The script derives all eight combinations of family, source sign, and half
sign.  It intentionally performs no resultant, saturation, or primary
decomposition.

## Conservative rational search

The script

```text
code/order44_from_order22_rational_sieve.m
```

builds a mask on `P^1(F_p)` for each sign.  A residue is allowed if it is:

- an exact open marked-class half;
- a singular or collapsed-order open reduction;
- one of `s=0,+/-1,infinity`.

Only a good open residue whose exact marked class is not divisible by `2`
is rejected.  Every rational survivor is converted to an integral square
model and tested with Magma's exact `IsDivisibleBy(D,2)`.

### Height `300`

Both presentations and both signs were searched using primes through `31`:

```text
reduced rational s values       109591
family/sign trials              438364
final conservative survivors      2172
degenerate exact boundaries            8
smooth exact divisibility tests      2164
exact errors                           0
exact halves                           0
```

### Height `1000`

By the exact translation identity, it suffices to use Flynn with both
signs.  Adding the conservative masks at `37,41,43,47` kept the exact stage
small:

```text
reduced rational s values      1216767
sign trials                    2433534
survivors after p=31             13270
survivors after p=37              6282
survivors after p=41              1820
survivors after p=43               628
final survivors after p=47         282
degenerate exact boundaries            4
smooth exact divisibility tests       278
exact errors                            0
exact halves                            0
```

This completed in `27.4` seconds under the `8 GiB` hard limit.  A prior run
with only primes through `31` reached the five-minute timeout during exact
checks and is not counted.

Full command lines and compact output are recorded in
`data/order44_from_order22_scout_2026_07_13.txt`.

## Recommendation

Do not spend the next example-producing sprint on a larger undirected
order-`44` height search.  The good open is already forced into boundary
disks simultaneously at `5` and `17`, and height `1000` is empty after
exact checks.

If this target is revisited, the next justified calculation is not a larger
box.  It is a weighted second- and third-order saturation of the three
square boundary charts, followed by normalization of any surviving
one-dimensional component.  Until that geometry is available, the
`[2,60]` and general order-`60` lanes remain better candidates for producing
a new high-torsion example.
