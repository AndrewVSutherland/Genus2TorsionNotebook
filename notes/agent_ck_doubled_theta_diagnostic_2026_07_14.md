# Doubled-theta loci on the Clebsch--Klein family

Date: 2026-07-14.

## Outcome

The proposed `2P-K` restriction gives three genuine, nonempty special loci
inside the universal source-halving covers.  None has a rational point in the
primitive Clebsch--Klein height-`100` box, but all three have exact smooth-open
points over some tested finite field.  Thus the idea is not killed by an
identity or by a universal good-reduction obstruction.

The exact diagnostic is

```text
code/agent_ck_doubled_theta_diagnostic.m
```

and the compact run ledger is

```text
data/agent_ck_doubled_theta_diagnostic_2026_07_14.txt.
```

Every run used at most `4 GiB` of address space and at most five minutes.

## Explicit marked classes

On the Clebsch--Klein surface put

```text
sum_i r_i = 0,       sum_i r_i^3 = 0,
C: y^2 = x*product_i(x-r_i^2).
```

Send the distinguished branch point `x=0` to infinity by

```text
X=1/x,       Z=y*X^3.
```

The odd model is

```text
Z^2 = g(X) = product_i(1-r_i^2*X).
```

Elkies's atypical order-`5` class becomes the completely explicit class

```text
T = (0,1)-infinity = J![X,1].
```

The program asserts that `T` is nonzero and killed by `5` on every smooth
fiber it uses.  If `a_i=r_i^2`, representative rational `2`-classes are

```text
S01 = J![X-1/a1,0],
S12 = J![(X-1/a1)*(X-1/a2),0].
```

An order-`10` generator has the form `a*T+S`, with `a=1` or `2` up to
sign.  The case `T+S01` cannot be doubled theta: geometrically its unique
effective degree-two representative is `P+W_i`, with distinct support.
The three remaining cases up to the `S_5` action are therefore

```text
T+S12,       2*T+S01,       2*T+S12.
```

For `D=[u,v]` on the odd model,

```text
D=2*P-K=2*(P-infinity)
```

is equivalent to

```text
Degree(u)=2,       Discriminant(u)=0.
```

For every reported finite hit, the code recovers the double root `alpha`,
puts `beta=v(alpha)`, constructs `Q=(alpha,beta)-infinity`, and asserts

```text
beta^2=g(alpha),       2*Q=D,       Order(D)=10.
```

Thus these are actual doubled-theta points, not merely discriminant zeros in
an unreduced presentation.

## Relation with the existing halving covers

If `D=a*T+S=2*Q`, then

```text
S=5*D=2*(5*Q).
```

Consequently every doubled-theta point lies on the exact Stoll--Zarhin
halving cover for `S`.  The old finite restrictions apply verbatim:

```text
2*T+S01:              forced to CK boundary at 11,19,29,31;
T+S12 and 2*T+S12:    forced to CK boundary at 11,19,23.
```

This implication is independently checked at every finite doubled-theta hit
using both the four-radicand criterion and the identity `2*(5Q)=S`.

These are forced bad-reduction disks, not established local impossibilities.
For orbit `12`, the prior work certifies smooth-open points on the full
halving cover over `Q_11`, `Q_19`, and `Q_23`; however, those displayed
Hensel branches were not tested against the additional doubled-theta
condition.  For orbit `01`, the earlier projective calculation proves local
solubility at `3` and `7`, but the four forced primes `11,19,29,31` have not
been resolved p-adically for this stricter locus.  Hence no existing result
proves or disproves local solubility of the doubled-theta curves at all of
their forced primes.

The radius-`10^6` CRT/lattice run concerned one selected product of three
orbit-`12` Hensel disks.  Its failure to find even one exact square radicand
also excludes both orbit-`12` doubled-theta cases in that selected rational
coset, but says nothing about the other local states.

## Finite-field census

The first two count columns below reproduce the exact marked source-halving
counts.  The last three count the stricter doubled-theta cases.

| `p` | open CK | general `01` | general `12` | `T+S12` | `2T+S01` | `2T+S12` |
|---:|---:|---:|---:|---:|---:|---:|
| 11 | 24   | 0  | 0   | 0  | 0  | 0  |
| 19 | 120  | 0  | 0   | 0  | 0  | 0  |
| 23 | 240  | 24 | 0   | 0  | 24 | 0  |
| 29 | 480  | 0  | 24  | 12 | 0  | 24 |
| 31 | 504  | 0  | 24  | 12 | 0  | 0  |
| 37 | 720  | 48 | 60  | 24 | 48 | 48 |
| 41 | 1104 | 72 | 120 | 48 | 24 | 96 |
| 43 | 960  | 24 | 72  | 36 | 0  | 0  |

In particular, all three doubled-theta cases occur on smooth fibers over a
tested finite field.  Examples are

```text
p=23:  2T+S01,  r=(1,4,12,14,15),  P=(10,14);
p=29:  T+S12,   r=(1,5,2,10,11),   P=(2,9);
p=29:  2T+S12,  r=(1,5,2,10,11),   P=(8,12).
```

## Exact bounded rational census

Rational mode uses the repository's validated quadratic enumeration of
primitive integral CK tuples, removes permutation/global-sign duplicates by
the unordered set of five squares, and tests all `5+10+10=25` relevant
classes on every source.

| height | distinct sources | classes tested | hits |
|---:|---:|---:|---:|
| 30  | 11 | 275  | 0 |
| 100 | 94 | 2350 | 0 |

The source counts agree with the earlier CK sweeps.  Moreover, the complete
height-`2000` source-halving search already tested all nonzero `2`-classes
on `10,290` sources and found no divisible class.  Since doubled theta
implies divisibility of `S`, that older result logically excludes all three
special cases throughout its height box even though it did not compute their
Mumford discriminants separately.

## Decision

The doubled-theta idea survives as a realistic construction, but a larger
blind CK height box is not justified.  The next useful computation is
symbolic:

1. use the rational two-parameter CK chart;
2. compute the Mumford `u` for the three displayed classes;
3. factor `Discriminant(u)` after imposing the CK equations and removing
   `r_i=0` and `r_i^2=r_j^2`;
4. normalize the resulting curves and audit their forced p-adic boundary
   disks before any larger rational search.

This is strictly an example-finding sublocus.  Its failure would not close
the full degree-`16` source-halving cover.
