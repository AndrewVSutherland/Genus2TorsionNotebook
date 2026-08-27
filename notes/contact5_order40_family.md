# Contact-5 plus 4-torsion family

Start from the quintic-contact 5-torsion construction

```
h(x) = 1 + a*x + b*x^2,
f(x) = h(x)^2 - (1+a+b)^2*x^5.
```

The exact search found the halving component

```
b = (a^2 - 1)/2.
```

Writing `t=a`, this gives

```
h = 1 + t*x + ((t^2 - 1)/2)*x^2,
f = h^2 - ((t + 1)^4/4)*x^5.
```

The discriminant factors as

```
(t + 1)^11*(t + 3)^4*(32*t^3 + 152*t^2 + 173*t + 37)/256.
```

For every nonsingular specialization, the contact identity gives a rational
5-torsion class and

```
H = [x^2 + 2*x/(t+1), (t+2)*x + 1]
```

satisfies `2H = [x-1,0]`, so `H` has order `4`.  Thus the family gives
rational torsion containing a point of order `20`.

The order-40 search in this family is implemented in

```
code/contact5_order40_family_search.m
```

Runs:

```
magma -b height:=35 code/contact5_order40_family_search.m
magma -b height:=50 max_hits:=10 code/contact5_order40_family_search.m
```

Both found the same order-40 specialization:

```
t = -1/3,  b = -4/9
y^2 = -324*x^5 + 1296*x^4 + 1944*x^3 - 5103*x^2 - 4374*x + 6561
J(Q)_tors = [40]
```

Magma also gives an irreducible Frobenius certificate at `p=7`.

The order-8 half is not divisible by 2 for this specialization, so this
particular curve does not give order `80`.


## Symbolic order-40 condition

The order-40 condition can be expressed by the standard 2-descent criterion:
the Mumford `u`-polynomial of

```
H = [x^2 + 2*x/(t+1), (t+2)*x + 1]
```

must be a square in the etale algebra of the branch polynomial.

The rational Weierstrass factor `x=1` gives the first necessary condition

```
(t+3)/(t+1) = r^2,        t = (3-r^2)/(r^2-1).
```

Put `s = r^2 - 1` and scale `x = s*z`.  The quartic factor becomes

```
4*z^4 + 4*(2-s)*z^3 + 4*(2-s)*z^2 + (4-s)*z + 1,
```

and the remaining condition is that `z*(z+1)` be a square modulo this quartic.

On the degree-2 square-root component

```
S = B*z^2 + C*z + D,
```

the Groebner basis gives

```
B = D*(s^2 - 3*s - 2)
2*C = D*(3*s^2 - 7*s - 4)
8*D^2 = 3*s^2 - 6*s - 1
(s - 3)*s*(s + 1) = 0.
```

Here `s=0` is the `t=infinity` boundary and `s=-1` gives the singular
member `t=-3`.  The only nonsingular point on this component is

```
s = 3,  r = +/-2,  t = -1/3.
```

This recovers the `[40]` curve above.  The full square-root cover is larger:
eliminating down to the `(D,s)` plane gives one irreducible equation with
`deg_D = 16`, `deg_s = 18`, total degree `34`, and `114` terms.  Exact
rational search on this full cover through height `50` found no other
order-40 specializations.

The symbolic helper is

```
sage --python code/contact5_order40_symbolic.py
```


## Remaining cover analysis

The remaining full cover is implemented in

```
sage --python code/contact5_order40_cover_analysis.py summary
sage --python code/contact5_order40_cover_analysis.py search 300
sage --python code/contact5_order40_cover_analysis.py finite
```

After eliminating the square-root coefficients and quotienting by the sign
symmetry `D -> -D`, the equation is an irreducible plane curve

```
Q(Y,s) = 0,        Y = D^2,
```

with

```
deg_Y = 8,  deg_s = 18,  total degree = 26,  terms = 114.
```

The known order-40 point appears in the fiber

```
s = 3,  Y = 1,  r = +/-2,  t = -1/3.
```

The affine singular ideal of `Q(Y,s)` is zero-dimensional of degree `41`,
with no rational affine singular points.  The singular `s`-coordinates are
the roots of two irreducible factors of degrees `12` and `17`; over those
fields the singularities are ordinary nodes by the Hessian test.  Thus there
are `24 + 17 = 41` affine nodes.

The Newton polygon has vertices

```
(0,0), (2,0), (8,6), (8,18), (0,10),
```

with `76` interior lattice points.  Four of the five boundary edge
polynomials are not squarefree, so the boundary is still singular and the
exact geometric genus was not obtained by the quick Sage genus routine.  In
any case, the remaining curve is high genus; it is not a hidden genus `0` or
genus `1` problem.

The eliminated-cover rational search through height `300` in `r` checked
`109588` nonboundary parameters.  It found only

```
r = +/-2,  s = 3,  Y = 1,  t = -1/3.
```

There were only `6` fibers with any rational `Y` root, and only the two above
had `Y` a rational square.

Small-prime finite residue counts with both square conditions imposed are
nonzero:

```
p 5:  1 point
p 7:  1 point
p 11: 1 point
p 13: 2 points
p 17: 3 points
p 19: 4 points
p 23: 5 points
```

So there is no simple local obstruction; the remaining difficulty is genuinely
global/high-genus.


## Independent extra 2-torsion: `[2,20]`

This was the natural next test after finding the contact-5/order-20 family.
The order-20 point already uses the rational Weierstrass factor `x=1`; an
independent rational 2-torsion point appears exactly when the residual quartic
`f/(x-1)` is reducible over `Q`.

The search script is

```
magma -b height:=300 max_hits:=100 code/contact5_order20_extra2_search.m
```

The height-300 run checked `109589` smooth specializations and found `18`
reducible residual quartics.  Every one exact-checked to

```
J(Q)_tors = [2,20]
```

and every one had an irreducible Frobenius certificate in the tested primes.

The smallest hit is

```
t = -29/5,  b = 408/25
```

with residual quartic splitting as two quadratics:

```
f/(x-1) = (x^2 - 35/48*x + 25/144)
          *(x^2 - 5/18*x + 25/576)
```

An integral model is

```
y^2 = -51840000*x^5 + 104040000*x^4 - 73950000*x^3
      + 25890625*x^2 - 4531250*x + 390625,
```

with exact torsion `[2,20]` and a simplicity certificate at `p=11`.

A second quadratic-quadratic hit is

```
t = 139/5,  b = 9648/25,
J(Q)_tors = [2,20],  simple certificate at p=29.
```

Most of the other height-300 hits have residual factorization type `[1,3]`,
which still gives one independent rational 2-torsion class and exact torsion
`[2,20]`.


## Parametrized extra-2 search and an order-80 hit

The residual quartic can be made much easier by putting

```
u = t + 1,        y = u*x.
```

Then `f/(x-1)` is equivalent to

```
u*(y^4 + 4*y^3 + 8*y^2 + 8*y + 4) - 4*y*(y+1)^2.
```

Thus the `1+3` extra-2 locus is parametrized by

```
t = -(z^4 + 4*z + 4)/(z^4 + 4*z^3 + 8*z^2 + 8*z + 4).
```

The `2+2` locus is parametrized by

```
t = -(r^6 - 2*r^5 + 2*r^4 - 4*r^3 + 4*r^2 - 8*r + 8)
     /((r^2 - 2)^2*(r^2 - 2*r + 2)).
```

The targeted search is implemented in

```
magma -b height:=1000 threshold:=80 code/contact5_extra2_param_large_search.m
```

A small validation run already finds a genuine larger example on the linear
locus:

```
z = -1/7,
t = -8233/7225,
b = 7790832/52200625,
J(Q)_tors = [2,2,20].
```

For the original odd-degree model, the residual quartic factors as

```
f/(x-1) = (x - 7225/1296)*(x - 7225/7056)
          *(x^2 - 917575/4032*x + 52200625/28224).
```

A reduced minimal model for the same curve is

```
y^2 + (x^2 + x)*y = -391671*x^6 + 1894851*x^5 + 6846924*x^4
                    - 15133525*x^3 + 3904068*x^2
                    + 2625336*x + 254016,
```

and Magma computes exact torsion `[2,2,20]`.

The higher-threshold run through parameter height `1000` checked both
parametrized extra-2 loci:

```
checked_linear = 1216767,
checked_qq     = 1216767,
unique_t       = 2027806.
```

It found `186` modular survivors with gcd bound `160`; all `186` exact torsion
computations were still `[2,20]`.  Thus no torsion larger than order `80` was
found in this parametrized search through height `1000`.


### Geometric splitting check for the `[2,2,20]` example

Sage/Lombardo endomorphism tests on the odd-degree model with

```
t = -8233/7225
```

give

```
J.geometric_endomorphism_algebra_is_field(B=100) = True
J.geometric_endomorphism_ring_is_ZZ(B=100)       = True
```

so the Jacobian is geometrically simple; in particular it is not geometrically
split as a product of elliptic curves.

A small-prime witness in the Algorithm 4.10 test is `p=71`.  The Frobenius
polynomial is

```
x^4 + 2*x^3 + 14*x^2 + 142*x + 5041,
```

and the `12`th-power transform used in the Lombardo test is irreducible.


## Further search for additional order-80 or larger examples

The double-linear residual-quartic condition was isolated explicitly.  If the
scaled residual quartic has two rational roots `z,w`, write

```
s = z + w,        p = z*w.
```

The symmetric condition is

```
(4-p)*s^2 + (-2*p^2 + 4*p + 8)*s - p^3 + p^2 + 4*p + 4 = 0.
```

Its discriminant in `s` is `4*p*(p-2)^2`, so setting `p=r^2` leaves the final
condition

```
Y^2 = (r+1)*(r^2+2*r+2)*(r^3-r^2-4*r+2).
```

The reusable Magma search is

```
magma -b height:=10000 code/contact5_order80_double_linear_search.m
```

A fast rational-height search on this genus-2 curve through height `10000`
found only

```
r = -2, -1, 0, 1/3
```

up to the sign of `Y`.  The first three are boundary/collision/degenerate
points.  The only nondegenerate point is

```
r = 1/3,
(z,w) = (-1/7, -7/9),
t = -8233/7225,
J(Q)_tors = [2,2,20].
```

Magma gives `RankBounds(Jacobian(C)) = [1,1]` and torsion `[4]` for the
auxiliary genus-2 curve, so a trivial rank-zero Chabauty certificate is not
available from the default routines.

For torsion strictly larger than `80`, the pure modular sieve

```
magma -b height:=1000 threshold:=80 prime_bound:=251 \
    code/contact5_extra2_param_modsieve.m
```

checked both parametrized extra-2 loci with `52` good primes.  It checked

```
checked_linear = 1216767,
checked_qq     = 1216767,
unique_t       = 2027806,
```

and found

```
survivors_linear = 0,
survivors_qq     = 0.
```

Thus, in the parametrized contact-5/order-20 extra-2 route, no parameter of
height at most `1000` can have rational torsion order greater than `80`, and
no additional `[2,2,20]` point was found on the exact double-linear locus
through height `10000`.


## Attempt to add rational 3-torsion to the `[2,20]` family

The natural target is to add a rational `3`-torsion class to the contact-5
order-20 family on the independent-extra-2 loci.  This would give torsion of
order divisible by `120`, e.g. `[2,60]` if the `3`-part combines with the
order-20 factor.

The necessary good-reduction test is very restrictive: for every good prime
`p != 3`, rational `3`-torsion forces

```
3 | #J(F_p).
```

The search script is

```
magma -b height:=1000 prime_bound:=251 code/contact5_extra2_plus3_search.m
```

It uses the two parametrized `[2,20]` loci:

```
linear:  t = -(z^4 + 4*z + 4)/(z^4 + 4*z^3 + 8*z^2 + 8*z + 4),
qq:      t = -(r^6 - 2*r^5 + 2*r^4 - 4*r^3 + 4*r^2 - 8*r + 8)
              /((r^2 - 2)^2*(r^2 - 2*r + 2)).
```

Finite residue counts already show forced boundary behavior.  On the linear
branch, primes `7` and `13` have no good residues with `3 | #J(F_p)`.  On the
quadratic-quadratic branch, prime `11` has no good residues with
`3 | #J(F_p)`.

A weaker run with `prime_bound=73` and `height=1000` had six smooth survivors,
but exact torsion on all six was still `[2,20]`.  The stronger run with
`prime_bound=251` checked

```
checked_linear = 1216767,
checked_qq     = 1216767,
unique_t       = 2027806,
```

and found only the singular boundary survivor

```
z = -2,  t = -3.
```

There were no smooth survivors and hence no `[2,20] + 3` example through height
`1000` on either parametrized extra-2 locus.


### Component-wise obstruction-prime boundary analysis for `[2,20] + 3`

The component-wise boundary script is

```
magma code/contact5_extra2_plus3_boundary_analysis.m
```

It decomposes the bad residues at the obstruction primes and solves the cubic
contact equations directly on the bad fibers.

For the linear branch at `p=7`:

```
z = 0,6  -> t=-1: no cubic-contact solutions
z = 3,5  -> t=-3: 2 contact solutions each
```

For the linear branch at `p=13`:

```
z = 0,12 -> t=-1: no cubic-contact solutions
z = 5,11 -> t=-3: 2 contact solutions each
z = 4,7  -> pole/infinity chart
```

The pole chart uses `s=1/t` and the square-scaled boundary fiber

```
y^2 = x^4 - x^5.
```

At `p=13` this pole fiber has `4` cubic-contact solutions, `2` with nondegenerate
`q` and `gcd(q,f)=1` on the special fiber.

For the quadratic-quadratic branch at `p=11`:

```
r = 0,1,2   -> t=-1: no cubic-contact solutions
r = 3,5,7,8 -> t=-3/discriminant component: 8 contact solutions each
```

On all finite `t=-3` components the mod-`p` contact divisor meets the branch
locus (`gcd1 = 0`), but first-order lifts modulo `p^2` exist and many move off
`t=-3`:

```
linear p=7:   each t=-3 residue has 98 p^2 lifts, 84 off t=-3
linear p=13:  each t=-3 residue has 338 p^2 lifts, 312 off t=-3
qq p=11:      t=-3 residues have p^2 off-branch lifts
```

Conclusion: the boundary analysis kills all `t=-1` components, but it does not
kill the `[2,20] + 3` strategy locally.  The only live branches are:

```
linear: p=7 on t=-3, and p=13 on either t=-3 or the pole chart;
qq:     p=11 on the t=-3/discriminant component.
```

These are exactly the branches that should be used for any further targeted
CRT/high-height search.  The previous good-prime sieve through height `1000`
and prime bound `251` found no smooth survivor on these branches.


### Targeted live-boundary CRT search for `[2,20] + 3`

Using the live branches from the component-wise boundary analysis, I added

```
code/contact5_extra2_plus3_live_boundary_search.m
```

The branches are:

```
linear_tminus3:  z = 3 or 5 mod 7,  z = 5 or 11 mod 13
linear_pole:     z = 3 or 5 mod 7,  z = 4 or 7  mod 13
qq_tminus3:      r = 3,5,7,8 mod 11
```

The script enumerates only rational parameters satisfying these CRT conditions
and then applies the same necessary good-prime condition `3 | #J(F_p)` at all
good primes up to the chosen bound.

Validation at height `1000`, `prime_bound=251`, reproduced the broad-search
result: only the singular point

```
z = -2,  t = -3
```

survived.

The higher targeted run

```
magma -b height:=5000 prime_bound:=251 max_exact:=200 \
    code/contact5_extra2_plus3_live_boundary_search.m
```

checked

```
linear_tminus3:  primitive CRT parameters = 1,085,573
linear_pole:     primitive CRT parameters = 1,085,578
qq_tminus3:      primitive CRT parameters = 10,134,217
unique_t:        8,927,496
```

and again found only the singular boundary survivor

```
linear_tminus3: z = -2, t = -3.
```

There were no smooth survivors, hence no exact torsion computations and no
`[2,20] + 3` hit through height `5000` on the live boundary branches.

## Full residual splitting on the double-linear locus

To target full rational 2-torsion while preserving the order-20 class, I next
imposed that the remaining quadratic factor in the double-linear residual
quartic also split over `Q`.

On the scaled residual quartic

```text
u*(y^4 + 4*y^3 + 8*y^2 + 8*y + 4) - 4*y*(y+1)^2,
```

write two rational roots as `z,w`, with

```text
s = z+w,   p = z*w.
```

The double-linear condition is

```text
(4-p)*s^2 + (-2*p^2 + 4*p + 8)*s - p^3 + p^2 + 4*p + 4 = 0,
```

and the discriminant in `s` is `4*p*(p-2)^2`, so `p=r^2`.  This gives the
previous genus-2 condition

```text
Y^2 = (r+1)*(r^2+2*r+2)*(r^3-r^2-4*r+2).
```

Dividing the quartic by `(y^2-s*y+p)`, the remaining quadratic has discriminant
square exactly when

```text
W^2 = (r+2)*(r^2+2*r+2)*(r^3-4*r^2-2*r+4).
```

So full residual splitting is reduced to the simultaneous two-square condition
above.  I added the exact search

```text
code/contact5_full_split_search.py
```

which checks the two degree-6 square numerators for `r=a/b` and uses small-prime
residue filters before exact square tests.

The completed run

```text
python3 code/contact5_full_split_search.py \
    --height 10000 --prime-bound 97 \
    --out data/contact5_full_split_h10000.txt
```

checked

```text
121589943 rational parameters.
```

It found no nondegenerate candidate.  The only simultaneous-square points in
this box are boundary or singular:

```text
r = -2: boundary,
r = -1: root collision z=w=-1,
r = 0: singular t=-1.
```

In particular, the known nondegenerate double-linear point

```text
r = 1/3,
(z,w) = (-1/7,-7/9),
t = -8233/7225,
J(Q)_tors = [2,2,20]
```

does not satisfy the second square condition: it has

```text
(r+2)*(r^2+2*r+2)*(r^3-4*r^2-2*r+4) = 13825/729,
```

which is not a rational square.  Thus this search found no full rational
2-torsion `[2,2,2,20]` specialization on the contact-5/order-20 double-linear
locus through height `10000` in `r`.

## Intersecting order 40 with extra rational 2-torsion

The next target after full residual splitting was `[2,40]`: combine the
order-40 cover in the contact-5/order-20 family with an independent rational
2-torsion condition, i.e. reducibility of the residual quartic `f/(x-1)`.

I added the finite-field diagnostic

```text
code/contact5_order40_extra2_finite.py
```

It enumerates the order-40 cover

```text
Q(Y,s)=0,   s+1 square,   Y square,
```

then sets `t=(2-s)/s` and checks whether the residual quartic is reducible over
`F_p`.  The run

```text
sage --python code/contact5_order40_extra2_finite.py 101 \
    > data/contact5_order40_extra2_finite_p101.txt
```

showed that the open intersection is obstructed at small primes:

```text
p=7:  cover=1, open=1, extra2=0
p=11: cover=1, open=1, extra2=0
p=13: cover=2, open=2, extra2=0
```

For larger primes there are open finite-field intersections, so this is not a
global finite-field emptiness phenomenon.  But any rational `[2,40]` candidate
on this model must reduce to the boundary at all of `7,11,13`.

I then added the boundary-filtered rational search

```text
code/contact5_order40_extra2_boundary_search.py
```

The forced boundary residues for `s` are

```text
p=7:  s in {0,1,-1}
p=11: s in {0,-1}
p=13: s in {0,-1}.
```

Here the boundary filter includes `s=0`, `s=-1`, `s=infinity`, `Y=0`, and
`Y=infinity`; at these three primes the `Y=infinity` condition adds no extra
`s` residues beyond the displayed sets.

The run

```text
sage --python code/contact5_order40_extra2_boundary_search.py 3000 \
    data/contact5_order40_extra2_boundary_h3000.txt
```

checked

```text
10944748 rational r-values,
779208 forced-boundary survivors,
0 fibers with rational Y-roots,
0 square-cover points,
0 extra2 candidates.
```

Thus this `[2,40]` route is locally forced to boundary and, through height
`3000` in the order-40 parameter `r`, the forced-boundary rational search does
not even recover a rational point on the order-40 cover, let alone one with
extra rational 2-torsion.


## Intersecting order 40 with rational 3-torsion

The next target was the clean large-torsion goal

```text
[40] + [3]  ->  torsion order divisible by 120.
```

I added a finite-field diagnostic

```text
code/contact5_order40_plus3_finite.py
```

It enumerates the order-40 cover

```text
Q(Y,s)=0,   s+1 square,   Y square,
```

sets `t=(2-s)/s`, and tests the necessary good-reduction condition

```text
3 | #Jac(C_t)(F_p).
```

The run

```text
sage --python code/contact5_order40_plus3_finite.py 101 \
    > data/contact5_order40_plus3_finite_p101.txt
```

showed that the open intersection is obstructed at small primes:

```text
p=7:  cover=1, open=1, plus3=0
p=11: cover=1, open=1, plus3=0
p=17: cover=3, open=3, plus3=0
```

There are open `plus3` residues at other primes, for example `p=13,19,23`,
so this is not a finite-field emptiness phenomenon.  But any rational
`[40] + [3]` candidate in this contact-5/order-40 model must reduce to the
boundary or bad-reduction locus at all of `7,11,17`.

I then added the boundary-filtered rational search

```text
code/contact5_order40_plus3_boundary_search.py
```

The forced boundary residues for `s` are:

```text
p=7:  s in {0,1,-1}
p=11: s in {0,-1}
p=17: s in {0,9,12,-1}
```

Here the boundary set includes `s=0`, `s=-1`, `s=infinity`, `Y=0`,
`Y=infinity`, and finite bad-reduction `s`-values of the family.

The height-1000 validation run checked

```text
1216764 rational r-values,
100770 forced-boundary survivors,
2 fibers with rational Y-roots,
0 square-cover points,
0 plus3 survivors.
```

The larger run

```text
sage --python code/contact5_order40_plus3_boundary_search.py 3000 101 \
    data/contact5_order40_plus3_boundary_h3000_p101.txt
```

checked

```text
10944748 rational r-values,
911978 forced-boundary survivors,
2 fibers with rational Y-roots,
0 square-cover points,
0 plus3 survivors.
```

Thus the contact-5 order-40 route to torsion order `120` is strongly
boundary-obstructed.  In the forced boundary classes through height `3000`, the
search does not even recover a rational square point on the order-40 cover, so
there is no candidate requiring an exact torsion computation.

