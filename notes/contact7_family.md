# Contact-7 family

This records the first direct `7`-torsion scaffold that is not obtained from
split elliptic-curve gluing.

## Construction

Work with an odd genus-2 model `C: y^2 = f(x)` and the rational point at
infinity.  If `deg(h) <= 3`, then the function

```text
x*y - h(x)
```

has pole order `7` at infinity.  Forcing

```text
h(x)^2 - x^2 f(x) = -(x - 1)^7
```

gives `div(x*y-h)=7P-7 infinity`, where `P=(1,h(1))`, provided the curve is
smooth and `h(1) != 0`.

The condition that `x^2` divide `h^2 + (x-1)^7` is enforced by

```text
h = 1 - (7/2)*x + a*x^2 + b*x^3.
```

Thus

```text
f = (h^2 + (x - 1)^7)/x^2
```

is a two-parameter family with a rational divisor class of exact order `7` in
the samples checked by Magma.

Code:

```text
code/contact7_family_search.m
code/contact7_root_even_search.m
```

## Initial broad search

Sample verification:

```text
magma -b mode:="verify" code/contact7_family_search.m
```

For `(a,b)=(0,0),(1,0),(0,1),(2,-1),(-3,2)`, Magma confirms that
`J![x-1,h(1)]` has order `7`.

Integer box:

```text
magma -b height:=25 mode:="integer" progress_interval:=5000 max_exact:=10000 \
    code/contact7_family_search.m > data/contact7_integer_h25.txt
```

Result:

```text
checked 2601
smooth 2601
verified7 2601
reducible 0
exact_tests 0
hits 0
factor types: 5 -> 2601
```

Rational height 8:

```text
magma -b height:=8 mode:="rational" progress_interval:=20000 max_exact:=10000 \
    code/contact7_family_search.m > data/contact7_rational_h8.txt
```

Result:

```text
checked 7569
smooth 7541
verified7 7541
reducible 25
exact_tests 25
hits 25
torsion counts: [14] -> 24, [28] -> 1
```

Rational height 12:

```text
magma -b height:=12 mode:="rational" progress_interval:=20000 max_exact:=10000 \
    code/contact7_family_search.m > data/contact7_rational_h12.txt
```

Result:

```text
checked 33489
smooth 33431
verified7 33431
reducible 67
exact_tests 67
hits 67
torsion counts: [14] -> 65, [28] -> 2
```

The best cyclic example found in this search is

```text
a = 11/2, b = -7/2.
```

After scaling, the curve is

```text
C: y^2 = 4*x^5 + 21*x^4 - 70*x^3 + 79*x^2 - 42*x + 9.
```

Magma verification:

```text
factorization = (x - 3/4)*(x^4 + 6*x^3 - 13*x^2 + 10*x - 3)
torsion      = [28]
marked_order = 7 for J![x-1,-1]
```

At `p=5`,

```text
#J(F_5) = 28
L_p(T) = 25*T^4 + 2*T^2 + 1
```

and this quartic is irreducible over `Q`, giving a sufficient certificate that
the Jacobian is `Q`-simple.

## Rational-root subfamily and target 56

To force a rational Weierstrass point, set

```text
r = 1 - s^2,
h(r) = eps*s^7,  eps in {+1,-1}.
```

For `r != 0`, this gives

```text
a = (eps*s^7 - 1 + (7/2)*r - b*r^3)/r^2.
```

This is implemented in `code/contact7_root_even_search.m`.

The finite-field target filter was corrected to use only the prime-to-`p` part
of the target at residue characteristic `p`.  For example, at `p=7` and target
`56`, the necessary condition is divisibility by `8`, not by `56`.

Validation at target `28`, height `8`:

```text
magma -b height:=8 target:=28 progress_interval:=20000 max_exact:=1000 \
    code/contact7_root_even_search.m > data/contact7_root_target28_h8.txt
```

Result:

```text
checked 15138
smooth 14604
verified7 14604
target_survivors 550
exact_tests 550
hits 24
torsion counts: [14] -> 526, [2,14] -> 22, [28] -> 2
```

Target `56`, height `12`, with the corrected prime-to-`p` filter:

```text
magma -b height:=12 target:=56 progress_interval:=50000 max_exact:=5000 \
    code/contact7_root_even_search.m > data/contact7_root_target56_h12_corrected.txt
```

Result:

```text
checked 66978
smooth 65866
verified7 65866
target_survivors 4
exact_tests 4
hits 0
torsion counts: [14] -> 4
```

So the contact-7 rational-root search quickly produces simple cyclic `[28]`,
but the nearby target `[56]` is strongly filtered out in the height-12 box.

## Conclusion

This is a useful new scaffold: it gives a direct two-parameter family with
rational `7`-torsion and simple examples with cyclic torsion `[28]`.  It is
more promising than imposing `7` on the earlier even-torsion surfaces, which
were repeatedly forced to bad boundary strata.

For torsion larger than the existing simple order-40 example, the next rational
step is not another blind height search.  The better next move is to derive the
halving condition for the rational Weierstrass point inside the rational-root
subfamily, turning the `[28]` condition into an explicit curve/surface, and then
test that space for a further halving to `[56]`.


## First-halving surface

Let `r = 1 - s^2` and impose `h(r)=eps*s^7`, with `eps = +/-1`.  Then
`f(r)=0`, so the rational Weierstrass class `(r,0)-infinity` gives the visible
2-torsion.

Translate `X = x-r` and write

```text
f(X+r)/X = X^4 + c3*X^3 + c2*X^2 + c1*X + c0.
```

Zarhin's formula says this class is divisible by `2` over `Q` exactly when the
quartic is

```text
prod_i (X + t_i^2)
```

with rational elementary symmetric functions `u,v,w,z` of the `t_i`.  Equating
coefficients gives

```text
c3 = u^2 - 2*v
c2 = v^2 - 2*u*w + 2*z
c1 = w^2 - 2*v*z
c0 = z^2
```

The corresponding order-4 half is represented by

```text
Q     = X^2 - v*X + z
alpha = (u*v - w)*X - u*z.
```

The equations are triangular: `c0=z^2` determines `b` linearly, then
`v=(u^2-c3)/2`, then `w=(v^2+2z-c2)/(2u)`, and only the `c1` equation remains.

Code:

```text
code/contact7_halving_surface_search.m
code/contact7_halving_surface_enum.py
code/contact7_halving_surface_finite.m
```

The known simple `[28]` example is recovered from the surface by

```text
s = 1/2, u = 4, z = -3/16, eps = -1.
```

Magma check:

```text
magma -b mode:="check_known" code/contact7_halving_surface_search.m
```

Output:

```text
ok true
a 11/2
b -7/2
r 3/4
v 7/2
torsion [28]
pass56 false
divisible8 false
simple true, pcert 5
```

## Second-halving local diagnostic

To test the route from `[28]` to `[56]`, I enumerated the first-halving surface
over small finite fields and checked whether the explicit order-4 half `H4` is
itself divisible by `2` in `J(F_p)`.

Run:

```text
magma code/contact7_halving_surface_finite.m \
    > data/contact7_halving_surface_finite.txt
```

The key obstruction is at `p=5`:

```text
p 5 checked 250 root_good 48 surface 4 exact_h4 4 h4_divisible_by_2 0
NO OPEN H4 HALVES
surface_samples:
  <s,u,z,eps,a,b,r,#J,ord(D7)>
  <2,1,2,1,3,4,2,28,7>
  <2,4,2,1,3,4,2,28,7>
  <3,1,2,4,3,4,2,28,7>
  <3,4,2,4,3,4,2,28,7>
```

Thus, on the open first-halving surface, the second halving needed for `[56]`
is obstructed modulo `5`.  Any rational `[56]` point on this contact-7 route
would have to reduce to the boundary at `5`.

For comparison, the obstruction is not universal across all small primes:
there are open second-halving residues at `p=7,11,13,17,19,...`.  The recurring
pattern is therefore again a small-prime boundary problem, now specifically
the `p=5` boundary of the contact-7 first-halving surface.

## Root-based `p=5` boundary pass

The broad Magma search on the first-halving surface was too slow because it
loops over all triples `(s,u,z)`.  I added

```text
code/contact7_halving_surface_roots.py
```

which eliminates the `u` loop.  For fixed `(s,z,eps)`, the remaining equation
is a quartic in `T = u^2`:

```text
((T-c3)^2 + 8*z - 4*c2)^2 - 64*T*(z*(T-c3) + c1) = 0.
```

The script solves this quartic over `Q`, keeps only rational square roots, and
can require the candidate to reduce to the forced `p=5` boundary.

The first completed boundary-filtered run was

```text
python3 code/contact7_halving_surface_roots.py \
    --s-height 8 --z-height 20 --require-boundary \
    --out data/contact7_halving_surface_roots_boundary5_s8_z20.txt
```

with output summary

```text
checked=85848 quartics=85848 root_t=8 square_t=2 hits=0 boundary_hits=0.
```

For comparison, the same box without the boundary filter found exactly four
rows, all the known `[28]` example and sign symmetries:

```text
s = +/-1/2, u = +/-4, z = -3/16, eps = -sign(s),
a = 11/2, b = -7/2, r = 3/4, v = 7/2, w = +/-1.
```

Magma verifies the representative

```text
s = 1/2, u = 4, z = -3/16, eps = -1
```

as

```text
torsion [28]
pass56 false
divisible8 false.
```

So the first serious `p=5` boundary pass found no route from the direct
contact-7 `[28]` scaffold to `[56]` in the small rational chart.  The remaining
way to continue this line would be a genuine local chart expansion at the
specific `p=5` boundary components, rather than another undirected height
increase.

A larger completed boundary-filtered pass was

```text
python3 code/contact7_halving_surface_roots.py \
    --s-height 8 --z-height 40 --require-boundary \
    --out data/contact7_halving_surface_roots_boundary5_s8_z40.txt
```

with output summary

```text
checked=329112 quartics=329112 root_t=8 square_t=2 hits=0 boundary_hits=0.
```

Thus increasing the `z` height from `20` to `40` did not find any candidate on
the forced `p=5` boundary.  The only rational-square `u^2` roots encountered
in this larger box are again the open `[28]` roots, which cannot lift to `[56]`
by the `p=5` finite-field obstruction.


## Local `p=5` open-neighborhood check for `[56]`

I then made the local check around the four open first-halving surface points
over `F_5`:

```text
<s,u,z,eps> = <2,1,2,1>, <2,4,2,1>, <3,1,2,-1>, <3,4,2,-1>.
```

Code:

```text
code/contact7_halving_surface_local5.py
```

Run:

```text
python3 code/contact7_halving_surface_local5.py --max-power 3 \
    > data/contact7_halving_surface_local5.txt
```

The first-halving surface itself lifts smoothly at these points:

```text
mod 5:   1 lift per base point
mod 25:  25 lifts per base point
mod 125: 625 lifts per base point
```

But the second-halving count is zero throughout these good open
neighborhoods.  The reason is already visible modulo `5`: the finite-field
diagnostic gives

```text
p 5 checked 250 root_good 48 surface 4 exact_h4 4 h4_divisible_by_2 0
```

Since the curve has good reduction and `2` is prime to `5`, any `5`-adic half
of the order-4 class would reduce to a half in `J(F_5)`.  No such half exists.
So these four open neighborhoods cannot contain a rational `[56]` example.

Conclusion: the genuine remaining `[56]` work is not around the four open
surface points.  It must be on excluded boundary charts of the first-halving
parametrization, for example where `s=0`, `s^2=1`, `u=0`, `eps*s=-1`, or a
denominator becomes nonintegral at `5`.


I then checked the most obvious excluded finite chart, `u=0`, without using
the solved formula for `w`.  The script

```text
code/contact7_halving_boundary5_finite.m
```

enumerates the unsolved symmetric equations over `F_5` in the open root chart.
Run:

```text
magma -b p:=5 code/contact7_halving_boundary5_finite.m \
    > data/contact7_halving_boundary5_finite_p5.txt
```

Result:

```text
root_checked 20
root_good 12
surface 4
surface_u0 0
exact_h4 4
exact_h4_u0 0
h4_divisible_by_2 0
h4_divisible_by_2_u0 0
```

So `u=0` contributes no first-halving points at all in the good open root
chart over `F_5`.

Finally I removed the `s` parametrization entirely and enumerated finite
contact-7 residues `(a,b,r)` with `f(r)=0`, then imposed the unsolved
first-halving equations:

```text
code/contact7_root_halving_all_finite.m
```

Run:

```text
magma -b p:=5 code/contact7_root_halving_all_finite.m \
    > data/contact7_root_halving_all_finite_p5.txt
```

Result:

```text
curve_checked 25
curve_good 18
root_good 10
surface 2
exact_h4 2
h4_divisible_by_2 0
TAG_ROOTS: open_s 6, r=0 4
TAG_SURFACE: open_s 2
TAG_EXACT_H4: open_s 2
TAG_DIV2: empty
```

Thus no finite good-reduction boundary residue over `F_5` supports the second
halving.  The `r=0` finite boundary has roots but no first-halving surface
points, while `s=0` is absent after requiring `h(1) != 0` and good reduction.
The remaining `[56]` possibility, if any, must therefore involve genuinely
bad or nonintegral reduction at `5`, not just a missing finite affine chart.


## `5`-adic valuation scan for the remaining `[56]` boundary

After closing the finite good-reduction boundary over `F_5`, I added a
valuation scan for the root-eliminated first-halving equation.  It classifies
both input quartics and actual rational first-halving candidates by the
`5`-adic valuations of

```text
s, u, z, a, b, r = 1-s^2, eps*s + 1.
```

Code:

```text
code/contact7_halving_valuation5_scan.py
```

Validation box:

```text
python3 code/contact7_halving_valuation5_scan.py \
    --s-height 8 --z-height 20 \
    --out data/contact7_halving_valuation5_s8_z20.txt
```

Result:

```text
checked 85848
quartics 85848
root_t 8
zero_t 2
square_t 2
hits 4
boundary_hits 0
```

The four hits are exactly the known open `[28]` point and sign variants:

```text
s = +/-1/2, u = +/-4, z = -3/16, a = 11/2, b = -7/2, r = 3/4.
```

Larger valuation pass:

```text
python3 code/contact7_halving_valuation5_scan.py \
    --s-height 12 --z-height 40 --progress 100000 \
    --out data/contact7_halving_valuation5_s12_z40.txt
```

Result:

```text
checked 705240
quartics 705240
root_t 8
zero_t 2
square_t 2
hits 4
boundary_hits 0
```

The scan explicitly included many boundary-looking input signatures, including

```text
v5(s) = 1          // s == 0 mod 5
v5(s) = -1         // s nonintegral at 5
v5(r) = 1          // r = 1-s^2 == 0 mod 5
v5(eps*s+1) = 1    // denominator boundary
```

and also nonintegral `z` cases in the smaller validation box.  None produced a
rational square `u^2` and hence none produced a rational first-halving point.
All rational square roots occurred in the open-integral valuation class.

Conclusion: within these valuation boxes, the missing `[56]` boundary is not
showing up even as a first-halving point.  Together with the finite `F_5`
checks, this makes the contact-7 `[28] -> [56]` route look essentially dead
unless there is a much more singular nonintegral model outside this chart.


## Rational-root subfamily plus 3-torsion

I also tested whether the contact-7 rational-root family can be combined with an independent rational `3`-torsion condition.  This would give torsion containing the prime-to-characteristic target `42`, and possibly `84` if a further `2`-part appears.

Code:

```text
code/contact7_root_plus3_search.m
```

The finite-field diagnostic enumerates the `(s,b,eps)` chart and asks whether the prime-to-`p` part of `42` divides `#J(F_p)`:

```text
magma -b mode:=finite code/contact7_root_plus3_search.m > data/contact7_root_plus3_finite.txt
```

There is local room at every odd prime tested after `p=3`; for example:

```text
p 5  good 12   pass_target 2
p 7  good 40   pass_target 4
p 11 good 142  pass_target 48
p 13 good 208  pass_target 68
p 17 good 414  pass_target 144
```

The chart has no good finite points modulo `3`, so `p=3` is again a bad boundary prime for this parametrization rather than a useful open-chart filter.

Rational height `12`:

```text
magma -b mode:=search height:=12 progress_interval:=20000 max_exact:=1000 code/contact7_root_plus3_search.m > data/contact7_root_plus3_h12.txt
```

Result:

```text
checked 66978
smooth 65866
verified7 65866
target_survivors 0
exact_tests 0
hits 0
```

Rational height `20`:

```text
magma -b mode:=search height:=20 progress_interval:=50000 max_exact:=1000 code/contact7_root_plus3_search.m > data/contact7_root_plus3_h20.txt
```

Result:

```text
checked 522242
smooth 519156
verified7 519156
target_survivors 0
exact_tests 0
hits 0
```

The first-kill distribution at height `20` was spread across the small prime filters, led by `p=5` and `p=7`:

```text
p 5:  114232
p 7:  145274
p 11: 100260
p 13: 67782
p 17: 42156
```

Conclusion: this route is locally plausible in finite-field samples, but the low-height rational chart is very sparse.  In the tested box there were no curves even satisfying the necessary reduction conditions for `42`, so this is currently weaker than the direct contact-7 `[28]` scaffold and the first-halving `[56]` analysis.


### `p=3` boundary diagnostic for the rational-root plus 3 route

The finite root chart has no good open points modulo `3`, so I separated the
height search by `3`-adic boundary branch.

Code:

```text
code/contact7_root_plus3_boundary3.m
```

Finite `p=3` summary:

```text
magma -b mode:=finite code/contact7_root_plus3_boundary3.m \
    > data/contact7_root_plus3_boundary3_finite.txt
```

There are two relevant finite-looking branches:

```text
s = 0, r = 1:
  all b in F_3 have h(1)=0 and discriminant 0.

r = 0 cancellation branch, eps*s = 1, a -> 1:
  b = 0: h(1)=0, singular
  b = 1: h(1)=1, singular
  b = 2: smooth, #J(F_3)=14
```

Thus the only good `p=3` boundary branch is

```text
s = eps mod 3, b = -1 mod 3, a = 1 mod 3.
```

This branch is compatible with the prime-to-`3` target: `#J(F_3)=14`.
So the `p=3` boundary itself does not force failure of target `42`.

The branch-tagged height `20` run was

```text
magma -b mode:=search height:=20 progress_interval:=50000 \
    code/contact7_root_plus3_boundary3.m \
    > data/contact7_root_plus3_boundary3_h20.txt
```

Result:

```text
checked 522242
smooth 519156
survivors_away3 0
```

The good `p=3` branch occurred `36180` times in the height-20 box, split
equally between the two signs of `s`.  Every one was killed by an away-from-3
prime.  Aggregate first-kill counts on this good branch:

```text
p 5:  7488
p 7:  9162
p 11: 7162
p 13: 4644
p 17: 3504
p 19: 1814
p 23: 1352
p 29: 578
p 31: 238
p 37: 160
p 41: 36
p 43: 24
p 47: 14
p 53: 2
p 59: 2
```

Conclusion: the bad `p=3` behavior is mostly a chart/boundary issue, not a
decisive local obstruction.  The only viable `p=3` branch is explicit and
small, but it still gives no height-20 candidate satisfying the away-from-3
necessary conditions for rational `3`-torsion.
