# Exact halving covers over Elkies's Clebsch--Klein surface

This note derives the universal covers on which a nonzero rational `2`-class
in Elkies's `[2,2,2,10]` family acquires a rational half.  Such a point would
give

```text
J(Q)_tors contains [2,2,2,20], of order 160.
```

The derivation is certified by

```text
code/elkies22210_halving_covers.m
```

and the independent bounded exact audit is recorded in

```text
data/elkies22210_source_halving_audit_summary.txt.
```

## 1. The source surface and its two orbits

Put

```text
S:  r1+r2+r3+r4+r5 = 0,
    r1^3+r2^3+r3^3+r4^3+r5^3 = 0
```

in `P^4_r`.  This is the Clebsch--Klein cubic surface.  Write

```text
ai = ri^2,
C_r: y^2 = x * product_i (x-ai).
```

The curve is smooth on

```text
U = S \ { Delta = 0 },
Delta = product_i ri * product_{i<j}(ri^2-rj^2).
```

There the six finite branch roots are `0,a1,...,a5`, all rational and
distinct, so `J[2](Q) = (Z/2)^4`.  The `15` nonzero `2`-classes are the
unordered branch pairs.  The natural `S5` action has two orbits:

```text
5 classes  {0,ai},
10 classes {ai,aj}, i != j.
```

It is therefore enough to derive covers for

```text
T01 = {0,a1},       T12 = {a1,a2}.
```

All the other covers follow by permuting the `ri`.

## 2. Exact Stoll--Zarhin criterion for an even split sextic

Let

```text
C: y^2 = (x-A)(x-B) product_c (x-c)
```

be monic of degree `6`, where `c` runs through the four remaining branch
roots.  The branch pair `{A,B}` represents the class `B-A` after moving
`A` to infinity.  Set

```text
t        = (x-B)/(x-A),
D        = product_c (A-c),
lambda_c = (c-B)/(c-A).
```

The exact rational change of variables

```text
X = D*(x-B)/(x-A),
Y = D^2*(A-B)^2*y/(x-A)^3
```

gives the monic odd model

```text
Y^2 = X * product_c (X-D*lambda_c).
```

Stoll's full-split Kummer map, equivalently Zarhin's division formula, now
gives the exact criterion

```text
B-A is in 2J(Q)
    iff q_c = -D*lambda_c is a square in Q for every c.
```

Since `D` contains the factor `A-c`, every `q_c` is polynomial:

```text
q_c = -D*(c-B)/(c-A).
```

This common squareclass matters.  Merely requiring all quotients
`lambda_c/lambda_d` to be squares gives a degree-`8` necessary cover, but it
does **not** imply rational divisibility.  The missing common condition is a
fourth quadratic extension.  The certification script includes an integral
negative control whose four cross-ratios are

```text
2, 8, 18, 32
```

and hence have the same squareclass, while both the exact radicands and
Magma's `IsDivisibleBy` show that the class is not divisible.  This is why the
universal rational-halving covers below have degree `16`, not `8`.

## 3. Orbit I: a class involving the distinguished root `0`

Take `A=0`, `B=a1`, with remaining roots `a2,a3,a4,a5`.  Then

```text
D = a2*a3*a4*a5 = (r2*r3*r4*r5)^2.
```

For `j=2,3,4,5`, the exact radicand is

```text
-D*(aj-a1)/aj
  = ((r2*r3*r4*r5)/rj)^2 * (r1^2-rj^2).
```

After removing the visible nonzero square factor, the exact cover is

```text
H01:
  r1+r2+r3+r4+r5 = 0,
  r1^3+r2^3+r3^3+r4^3+r5^3 = 0,
  z2^2 = r1^2-r2^2,
  z3^2 = r1^2-r3^2,
  z4^2 = r1^2-r4^2,
  z5^2 = r1^2-r5^2.
```

This is a homogeneous model in `P^8` with all coordinates of weight `1`.
It is finite etale over `U`, has dimension `2`, and has generic degree
`2^4=16`.

There is also a useful real constraint: a rational point must have

```text
|r1| > |rj| for j=2,3,4,5.
```

Thus only the branch root of largest absolute `ri` can occur in this orbit at
a real point.

## 4. Orbit II: a pair of nonzero roots

Take `A=a1`, `B=a2`, with remaining roots `0,a3,a4,a5`.  Removing the
visible square factor `r2^2` from the radicand at `0`, and `r1^2` from the
other three radicands, gives

```text
G0 = -(a1-a3)(a1-a4)(a1-a5),
G3 =  (a3-a2)(a1-a4)(a1-a5),
G4 =  (a4-a2)(a1-a3)(a1-a5),
G5 =  (a5-a2)(a1-a3)(a1-a4).
```

The exact cover is

```text
H12:
  r1+r2+r3+r4+r5 = 0,
  r1^3+r2^3+r3^3+r4^3+r5^3 = 0,
  v0^2 = G0,
  v3^2 = G3,
  v4^2 = G4,
  v5^2 = G5.
```

Each `Gi` is homogeneous of degree `6` in the `ri`.  Thus `H12` sits
naturally in

```text
P(1,1,1,1,1,3,3,3,3),
```

where the `ri` have weight `1` and the `vi` have weight `3`.  It is again a
dimension-`2`, generic degree-`16`, finite etale cover of `U`.

The apparent simplification

```text
lambda_0 = a2/a1 = (r2/r1)^2
```

only reduces the three *relative* cross-ratio equations.  The common exact
Stoll squareclass is `G0`; it is the fourth equation and cannot be omitted.

For a direct search, the three relative conditions have the useful form

```text
Gj/G0 = -(aj-a2)/(a1-aj) = wj^2,       j=3,4,5.
```

Equivalently,

```text
aj = (a2-wj^2*a1)/(1-wj^2).
```

After scaling `a1=1` and writing `a2=q^2`, the additional requirement that
`aj=rj^2` is rational becomes the genus-one quartic

```text
wj^2 = (q^2-rj^2)/(1-rj^2).
```

One must still impose the common equation `v0^2=G0` and the two CK
equations.  There is a strong exact real sieve.  Sort the five `ai` and
orient the marked pair so that `a1>a2`.  Positivity of the three displayed
ratios puts every complementary `aj` outside `[a2,a1]`, so the marked pair
must be adjacent.  Positivity of `G0` says that the number of complementary
roots above `a1` is odd.  Consequently the only eligible unordered pairs
occupy sorted positions

```text
(1,2) or (3,4).
```

Thus only `2` of the `10` orbit-II classes on any real CK source can possibly
halve.  Likewise only the unique largest `ai` is eligible in orbit I, so an
exact source audit needs to call `IsDivisibleBy` on just `3`, rather than all
`15`, nonzero `2`-classes.

## 5. Degree, connectedness, and boundary

The four squareclasses in each representative are generically independent.
This can be seen directly from valuations on the collision divisors.

For `H01`, the factor `r1-rj` occurs to odd order only in the `j`th displayed
radicand, for `j=2,3,4,5`.

For `H12`, use the six collision factors

```text
a1-a3, a1-a4, a1-a5, a2-a3, a2-a4, a2-a5.
```

The four incidence rows are

```text
G0:  1 1 1 0 0 0
G3:  0 1 1 1 0 0
G4:  1 0 1 0 1 0
G5:  1 1 0 0 0 1.
```

They have rank `4` over `F2`.  Hence both function-field extensions are
connected multiquadratic extensions of degree `16`.

The script independently specializes at Elkies's point

```text
(r1,...,r5) = (1,-8,-7,5,9).
```

For `H01`, the four values are

```text
-63, -48, -24, -80,
```

and for `H12` they are

```text
92160, -28800, -149760, 19584.
```

Both sets have squareclass rank `4`, giving a second exact degree
certificate.

Every factor appearing in the radicands vanishes only on `Delta=0`.  Thus
both covers are unramified over the smooth genus-2 locus `U`.  Their projective
closures ramify and may acquire extra components over

```text
ri=0 or ri=+-rj.
```

Those are repeated-root/singular-curve strata and must be saturated away in
any rational-point search.  Boundary points are not torsion examples.

## 6. Explicit half on either cover

The cover variables give actual square roots of the exact radicands.

For `H01`, put

```text
u_j = (r2*r3*r4*r5/rj)*z_j,  j=2,3,4,5.
```

For `H12`, put

```text
u_0 = r2*v0,
u_3 = r1*v3,  u_4 = r1*v4,  u_5 = r1*v5.
```

After the change of variables in Section 2, the curve is

```text
Y^2 = X * product_{i=1}^4 (X+u_i^2).
```

Let `s_k` be the elementary symmetric functions in the four `u_i`.  Zarhin's
half has Mumford representation

```text
q(X)     = X^2-s2*X+s4,
alpha(X) = (s1*s2-s3)*X-s1*s4.
```

The script verifies symbolically that

```text
alpha^2 - X*product_i(X+u_i^2)
```

is divisible by `q`, and verifies on an exact rational specialization that

```text
2*[q,alpha] = [(X),0].
```

The `16` independent sign choices of the four square roots give the `16`
halves of the fixed nonzero `2`-class.  Simultaneously changing all four signs
negates the half.

## 7. Direct rational chart for the first orbit

The first cover has a useful chart that avoids re-enumerating Clebsch--Klein
base points.  On `r1 != 0`, scale `r1=1` and put

```text
rj = (1-tj^2)/(1+tj^2),
zj = 2*tj/(1+tj^2),       j=2,3,4,5.
```

This parametrizes all four equations `rj^2+zj^2=1` away from their standard
points at infinity.  The two remaining equations are

```text
sum_{j=2}^5 1/(1+tj^2) = 3/2,
product_{j=2}^5 (1-tj^2) = 16.
```

Indeed the first equation is just `1+sum rj=0`.  If `e_k` denotes the
elementary symmetric functions of the four `rj`, the cubic equation on this
linear locus is `e3=-e2`.  Evaluating their monic quartic at `-1` and `0`
shows that this is equivalent to

```text
product_j(1+rj) = product_j rj,
```

which becomes the displayed product equation after substituting the `tj`.
The script certifies these identities symbolically.  Another equivalent form
of the second equation is `sum rj*zj^2=0`, since
`rj^3-rj=-rj*zj^2`.

If the four `rj` are paired with pair sums `p,q` and pair products `m,n`, the
same two conditions become

```text
p+q=-1,       p*q=p*m+q*n.
```

There is an even more economical form for exact searching.  Put `sj=tj^2`
and write `E_k=e_k(s2,s3,s4,s5)`.  Clearing denominators in the two displayed
`t`-equations and taking two invertible linear combinations gives

```text
E4=E1+5,       E3=E2-10.
```

Pair `s2,s3` and set

```text
A=s2+s3,   B=s2*s3,   C=s4+s5,   D=s4*s5.
```

Away from `(B-1)*(B-A+1)=0`, the equations determine the second pair:

```text
den = (B-1)*(B-A+1),
C   = (B^2-10*B-A^2-4*A+5)/den,
D   = ((A+6)*B-(A^2+5*A+10))/den.
```

Thus an exact direct search chooses two rational `t`-parameters, computes
`C,D`, and asks whether the roots of `T^2-C*T+D` are distinct rational
squares.  This is a two-parameter search on the cover itself, rather than a
four-parameter search or a base-surface enumeration.  There is no omitted
smooth exceptional slice: `B=1` makes `r2^2=r3^2`, while
`B-A+1=(s2-1)*(s3-1)=0` makes one branch root zero.  Both lie on `Delta=0`.

This four-parameter/two-equation chart is the natural object for a direct
rational-point search on `H01`; a larger bounded enumeration of integral
Clebsch--Klein base points does not exploit the cover geometry.

## 8. Exact computational checks and current rational-point evidence

Running

```text
magma -b code/elkies22210_halving_covers.m
```

certifies:

```text
both symbolic fractional-linear transformation identities,
both simplified orbit formula lists,
generic squareclass rank 4 and degree 16 for both covers,
Zarhin's Mumford formula and exact doubling,
agreement with Magma IsDivisibleBy on all 15 classes of Elkies's source,
an exact positive even-sextic control,
an exact negative same-cross-ratio-squareclass control.
```

The separate source enumeration constructs all `15` exact Jacobian classes
and calls `IsDivisibleBy` directly.  Its current results are:

```text
known Elkies source:  15 classes,    0 divisible,
height 30:           165 classes,    0 divisible,
height 100:         1410 classes,    0 divisible.
```

All `94` height-`100` sources retained exact torsion `[2,2,2,10]`.  This is
an independent Magma cross-check.  The faster exact integral-projective
enumerator has since completed the primitive source box

```text
max_i |ri| <= 2000:
  10,290 unique smooth CK source curves,
  0 orbit-I cover points,
  0 orbit-II cover points.
```

Searches on the covers themselves, rather than on integral base points, give

```text
orbit I:  t-parameter height 150, 94,030,041 pairs, 0 points;
orbit II: rational (r2,r3)-height 100, 148,157,584 pairs, 0 points;
full rational CK chart: (t,m)-height 100, 148,230,625 pairs, 0 points.
```

Only three orbit-I pairs even completed to rational CK tuples in the first
run.  In the second, `15,004` smooth CK completions survived, but none made
all four exact `Gi` simultaneous squares.  The full chart covers the entire
labelled smooth CK open; the stated computation bounds its two inverse
parameters.  Full run records, near-misses, and the independent height-`50`
meet-in-the-middle check are in

```text
notes/elkies22210_ck_halving_cover_direct_search_2026_07_11.md
notes/elkies22210_source_halving_local_and_search_2026_07_11.md.
```

All of these are bounded negative evidence, not a global nonexistence
result.  The two explicit dimension-`2` covers above remain the correct
geometric objects for a serious search for `[2,2,2,20]` of order `160`.
