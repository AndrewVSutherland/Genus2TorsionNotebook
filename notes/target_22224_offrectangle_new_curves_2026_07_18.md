# New elliptic curves mined from the off-rectangle `A(2,2,2,8)` bank

Date: 2026-07-18

## Outcome

The 56 primitive positive off-rectangle points in `tor2228.txt` do **not**
lie on any of the four previously used rational curves, and in fact none lies
on Adam Logan's auxiliary rational first-leg cover in any ordering.  They do,
however, contain a much better structure:

* three repeated-coordinate genus-5 fibers, whose nine elliptic quotients have
  ranks between 2 and 4; and
* seven pair-scaling fibers which are themselves genus-one covers and have
  positive rank.

The seven pair-scaling fibers produce many genuinely new rational points on
the full signed `A(2,2,2,8)` cover.  A first bounded run found 160 new points
(145 in an interrupted rank-one `P1` walk and 15 in a height-capped run on
`P2`--`P7`).  Of these, 82 use positive signed branch coordinates.  No point
survived the direct good-prime test for rational 3-torsion, so there is not yet
a `[2,2,2,24]` example.

The smallest new **positive-sign** point found is

```text
[405616827844, 43646066521703, 277057639659506, 340515326975038].
```

Magma verifies that its Jacobian torsion is exactly `[2,2,2,8]` and that the
Jacobian is geometrically simple.  A strict root-power certificate is supplied
at `p=97`, with Weil polynomial

```text
9409*T^4 + 776*T^3 + 46*T^2 + 8*T + 1.
```

It is killed for the target by direct reduction (3 does not divide the
Jacobian order; the accumulated gcd bound is 128).

## 1. Adam inverse coordinates: a clean negative result

For an ordered tuple `(a,b,c,d)`, normalize

```text
x_i = (1/a_i) / sum_j(1/a_j)
```

and put

```text
A=x1+x2, B=x1+x3, C=x1+x4,
u=A/(1-A), v=B/(1-B), w^2=C*u*v/(1-C).
```

The last quantity is a rational square exactly when the symmetric ten-factor
character defining Adam's first-leg cover is a square.  It is nonsquare for
all 56 points.  Since the character is symmetric, changing the ordering does
not help.  Thus the off-rectangle bank is genuinely disjoint from this chart,
not merely missed by a poor choice of inverse coordinates.

Exact interpolation in the normalized branch chart `(b/a,c/a,d/a)` finds no
common relation of total degree 1, 2, or 3.  There are also no repeated
reciprocal pair-partition invariants.

## 2. Three shared-coordinate fibers

The only projective `3+1` coordinate-ratio coincidences are

| fixed triple | moving values | square parameter ratio |
|---|---:|---:|
| `(46,292,1679)` | `2, 722` | `1, 19^2` |
| `(71,3650,10366)` | `2116, 5476` | `1, (37/23)^2` |
| `(276,624,828)` | `13, 1573` | `1, 11^2` |

If the fixed triple is `(r,s,t)` and the first moving value is `d0`, write
`d=d0*T^2`.  Relative to the point `T=1`, the full fiber is

```text
Y_r^2 = (r+d0*T^2)/(r+d0),
Y_s^2 = (s+d0*T^2)/(s+d0),
Y_t^2 = (t+d0*T^2)/(t+d0).
```

It is a genus-5 `(Z/2)^3` cover of the `T`-line.  Each pair gives an elliptic
quotient

```text
Z^2=((r+d0*T^2)(s+d0*T^2))/((r+d0)(s+d0)).
```

The exact ranks of the three quotients on the three fibers are respectively

```text
(4,4,3), (4,3,4), (3,3,2).
```

An exact Mordell--Weil lattice walk through 61,138 quotient points recovered
only the six known full-cover points.  There was no new full-cover lift and no
3-primary survivor.  These quotients are high rank but the two residual square
conditions are very restrictive.

## 3. Seven pair-scaling elliptic fibers

There are exactly seven projective `2+2` coordinate-ratio coincidences.  After
choosing an ordering they have the form

```text
(a,b,c,d) -> (lambda*a,b,c,lambda*d).
```

The base and second parameter are:

| fiber | `(a,b,c,d)` at `lambda=1` | second `lambda` |
|---|---|---:|
| P1 | `(1,55,99,125)` | `1089/25=(33/5)^2` |
| P2 | `(20,225,304,380)` | `9` |
| P3 | `(13,276,624,828)` | `16` |
| P4 | `(41,256,800,1312)` | `6400/1681=(80/41)^2` |
| P5 | `(2,46,292,1679)` | `4` |
| P6 | `(18,158,711,1764)` | `6241/1764=(79/42)^2` |
| P7 | `(50,791,2800,5650)` | `196/25=(14/5)^2` |

The genus-one geometry is explicit.  Relative to `lambda=1`, put

```text
Ra=((a*lambda+b)(a*lambda+c))/((a+b)(a+c)),
Rb=((a*lambda+b)(d*lambda+b))/((a+b)(d+b)),
Rc=((a*lambda+c)(d*lambda+c))/((a+c)(d+c)).
```

The product `abcd` changes by `lambda^2`, so the full cover is exactly the
condition that `Ra,Rb,Rc` are squares.  Write

```text
p=Ya/Yb, q=Yc/Ya,
kp=(d+b)/(a+c), kq=(a+b)/(d+c).
```

Then

```text
lambda=(kp*c-b*p^2)/(d*p^2-kp*a)
```

and elimination gives

```text
A(p)q^2+B(p)=0,

A(p)=a(kp*c-b*p^2)+b(d*p^2-kp*a),
B(p)=-kq*d(kp*c-b*p^2)-kq*c(d*p^2-kp*a).
```

Thus

```text
z^2=-A(p)B(p),  z=q*A(p)
```

is an elliptic quotient of the full fiber.  A quotient point lifts precisely
when `Ra` is a square; the ratio equations then force `Rb` and `Rc` to be
squares as well.  The two bank points lift, and their difference has infinite
order on every one of the seven quotients.  Consequently all seven full
fibers are positive-rank genus-one curves (isogenous to the displayed
quotients).

Minimal quotient models are:

```text
P1: y^2       = x^3 + x^2 - 1621520*x + 676655700
P2: y^2+x*y   = x^3       - 30467519660*x + 1510620780362700
P3: y^2       = x^3 + x^2 - 10659112640*x + 210178158141300
P4: y^2       = x^3 - x^2 - 322452520*x - 682095494468
P5: y^2       = x^3 - x^2 - 2187523520*x + 29012594318400
P6: y^2       = x^3 + x^2 - 208814660*x + 26517397308
P7: y^2+x*y+y = x^3 + x^2 - 268718070*x + 558361642995.
```

### Exact `P1` certificate

For `P1`,

```text
A(p)=6820*p^2+396/5,
B(p)=-1375*p^2-27621/5,

z^2=9377500*p^4+37783944*p^2+10937916/25.
```

The base and second points are

```text
(p,z)=(1,34496/5), (27/125,931392/625).
```

The minimal elliptic curve has exact rank 1 and Mordell--Weil invariants

```text
[2,2,0].
```

Its generators are

```text
(955,0), (-1446,0), (4330,273600),
```

and the second bank point is exactly the free generator `(4330,273600)`.
The complete birational maps are printed reproducibly by
`code/target_22224_offrectangle_new_curves_P1_certificate.m`.

The interrupted `P1` walk wrote 147 distinct full lifts: two known points,
74 new positive-sign points, and 71 new points on the signed-negative chart.
All 147 fail the direct 3-primary filter.  The smallest new positive-sign
tuple has height about `4.89e31`:

```text
[391526449011144253734374417041,
 24364756619430415748480314825255,
 43856561914974748347264566685459,
 48940806126393031716796802130125].
```

### Bounded `P2`--`P7` run

With a height ceiling of `10^40`, a walk along the known infinite-order point
tested 984 quotient points.  Of these, 252 lifted to the full fibers, 225 were
above the height ceiling, and 27 were retained.  The retained rows comprise
12 known and 15 new curves; eight of the new curves have positive signed
coordinates.  None survives the direct 3-primary filter.

The six new curves of height below `10^18` were checked exactly.  All have
torsion `[2,2,2,8]` and geometrically simple Jacobian.  Two are on the positive
chart (`P5`); the four `P2/P7` examples use a signed-negative presentation.

## 4. Local 3-primary profiles

Direct finite-field profiles explain why blind multiplication is unproductive.
For each family there is a small prime at which no good specialization has
Jacobian order divisible by 3:

| family | prime | good lambdas | target lambdas |
|---|---:|---:|---:|
| P1 | 13 | 4 | 0 |
| P2 | 13 | 6 | 0 |
| P3 | 11 | 4 | 0 |
| P4 | 13 | 4 | 0 |
| P5 | 11 | 2 | 0 |
| P6 | 13 | 6 | 0 |
| P7 | 13 | 6 | 0 |

This does not rule out a rational target point: it forces such a point into a
boundary residue class at the indicated prime.  The correct continuation is
therefore an elliptic Mordell--Weil residue sieve imposing:

1. boundary reduction at the zero-target primes above; and
2. target or boundary residues at the remaining good primes.

That is much sharper than taking larger multiples and letting numerator and
denominator heights grow exponentially.

## 5. Artifacts

```text
code/target_22224_offrectangle_new_curves.py
code/target_22224_offrectangle_new_curves_fibers.m
code/target_22224_offrectangle_new_curves_pairfibers.m
code/target_22224_offrectangle_new_curves_pairfiber_profiles.m
code/target_22224_offrectangle_new_curves_P1_certificate.m
code/target_22224_offrectangle_new_curves_verify.m

results/target_22224_offrectangle_new_curves_analysis.txt
results/target_22224_offrectangle_new_curves_fibers.log
results/target_22224_offrectangle_new_curves_fiber_candidates.tsv
results/target_22224_offrectangle_new_curves_P1_compact.tsv
results/target_22224_offrectangle_new_curves_P1_compact.tsv.gz
results/target_22224_offrectangle_new_curves_pairfiber_small.tsv
results/target_22224_offrectangle_new_curves_pairfibers_small.log
results/target_22224_offrectangle_new_curves_verify.log
```

The uncompressed `P1` table is large because the rank-one coordinates grow
very quickly.  The `.gz` artifact is the most convenient exact archive of all
147 written rows.
