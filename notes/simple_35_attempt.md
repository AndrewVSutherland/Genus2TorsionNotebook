# Simple 35-torsion attempt

Goal: find a genus-2 curve over `Q` with `Q`-simple Jacobian and rational
torsion containing `35`.

The known easy construction of `35`/`70` torsion uses elliptic-curve gluing and
is not suitable here because the Jacobian is isogenous to a product.  The tests
below only use contact families that can have simple Jacobians.

## Contact-7 plus 5

Start from the direct contact-7 family

```text
h = 1 - (7/2)*x + a*x^2 + b*x^3
f = (h^2 + (x - 1)^7)/x^2.
```

This gives a rational order-7 divisor class.  A rational 5-torsion class would
force `5 | #J(F_p)` at every good prime `p != 5`.

Code:

```text
code/contact7_plus5_search.m
```

Finite-field diagnostic:

```text
magma -b mode:="finite" code/contact7_plus5_search.m \
    > data/contact7_plus5_finite.txt
```

Key result:

```text
p 3 total 9 good 5 pass5 0
```

Thus the open contact-7 surface is obstructed modulo `3`: any rational
`35`-torsion example on this route must have bad reduction at `3`.

Rational search:

```text
magma -b mode:="search" height:=20 progress_interval:=50000 max_exact:=500 \
    code/contact7_plus5_search.m > data/contact7_plus5_h20.txt
```

Result:

```text
checked 261121
smooth 260959
five_survivors 0
exact_tests 0
hits 0
```

First-kill counts:

```text
3  83761
7  89622
11 55710
13 18813
17 8264
19 3311
23 1044
29 314
31 93
37 24
41 3
```

So no height-20 point even survives the necessary 5-divisibility filters.

## Contact-5 plus 7

The dual direction starts from the quintic-contact 5-torsion family

```text
h = 1 + a*x + b*x^2
f = h^2 - (1+a+b)^2*x^5.
```

A rational 7-torsion class would force `7 | #J(F_p)` at every good prime
`p != 7`.

Code:

```text
code/contact5_plus7_search.m
```

Finite-field diagnostic:

```text
magma -b mode:="finite" code/contact5_plus7_search.m \
    > data/contact5_plus7_finite.txt
```

Key result:

```text
p 3 total 9  good 3  pass7 0
p 5 total 25 good 15 pass7 0
```

Thus the open contact-5 surface is obstructed at both `3` and `5`.  Any
rational `35`-torsion example on this route would have to reduce to the bad
boundary at both primes.

Rational search:

```text
magma -b mode:="search" height:=12 progress_interval:=20000 max_exact:=500 \
    code/contact5_plus7_search.m > data/contact5_plus7_h12.txt
```

Result:

```text
checked 33489
smooth 33300
seven_survivors 0
exact_tests 0
hits 0
```

## Simultaneous contact equations

I also tried a more structured subroute: keep the contact-7 family and impose
an independent point-contact 5 condition

```text
q(x)^2 - f(x) = const*(x-r)^5,    q = c0 + c1*x + c2*x^2.
```

The derivative equations give a three-parameter exact search.  The top
equation is

```text
c2^2 - b^2 = 5*r - 7.
```

Writing `c2-b=d`, `c2+b=e`, gives

```text
r  = (d*e + 7)/5
b  = (e-d)/2
c2 = (e+d)/2.
```

The next two derivative equations solve for `a` and `c0`; the remaining two are
checked exactly.

Code:

```text
code/contact7_contact5_point_enum.py
code/contact7_contact5_point_verify.m
```

Run:

```text
python3 code/contact7_contact5_point_enum.py --height 6 \
    --out data/contact7_contact5_point_candidates_h6.txt
```

Result:

```text
checked 95128
hits 0
```

A height-10 run was interrupted after `1,500,000` exact checks, also with
`hits=0`.

## Conclusion

No simple `35` example was found in these first direct searches.  More
importantly, both natural simple-contact surfaces have immediate small-prime
open obstructions:

```text
contact-7 plus 5: obstructed on the open surface modulo 3
contact-5 plus 7: obstructed on the open surface modulo 3 and modulo 5
```

The contact-7 direction remains the less bad of the two, since it only forces
one small-prime boundary.  The next serious step would be a `p=3` boundary
analysis of the contact-7 family, using the discriminant factorization

```text
Disc(f) = (2*a + 2*b - 5)^7 * Q5(a,b) / 256.
```

The visible boundary component `2*a+2*b-5=0` should be analyzed first.


## Contact-7 `p=3` boundary residues

For the contact-7 family modulo `3`, the bad affine residue classes are:

```text
(a,b) = (0,1): disc=0, h(1)=0, visible component 2a+2b-5=0
(a,b) = (1,0): disc=0, h(1)=0, visible component 2a+2b-5=0, also Q5=0
(a,b) = (1,1): disc=0, h(1)!=0, Q5=0
(a,b) = (2,2): disc=0, h(1)=0, visible component 2a+2b-5=0
```

The five good affine residue classes modulo `3` all fail the necessary
`5 | #J(F_3)` condition.  Thus the next boundary analysis should begin with
these four residue classes, with `(1,1)` being the only one where the marked
7-contact point does not itself degenerate modulo `3`.


## Boundary-at-3 analysis for contact-7 plus 5

Code:

```text
code/contact7_plus5_boundary3_search.m
code/contact7_plus5_boundary3_components.m
```

First I separated the `p=3` boundary classes and reran the necessary
`5 | #J(F_p)` filters away from `p=3`.

Full boundary search at height `20`:

```text
magma -b height:=20 progress_interval:=50000 max_exact:=200 \
    code/contact7_plus5_boundary3_search.m \
    > data/contact7_plus5_boundary3_h20.txt
```

Result:

```text
checked 261121
smooth 260959
boundary 177198
survivors 0
exact_tests 0
hits 0
```

The finite boundary classes at `p=3` had counts:

```text
finite:a=0,b=1,h1=0,L=0       16297
finite:a=1,b=0,h1=0,L=0,Q5=0  16297
finite:a=1,b=1,Q5=0           18224
finite:a=2,b=2,h1=0,L=0       18184
```

There are also infinity boundary classes from denominators divisible by `3`;
none produced a survivor either.  For every class, the first kills are mostly
at `7` and `11`, with smaller tails at later primes.

### Special fibers

The special-fiber diagnostic is in
`data/contact7_plus5_boundary3_components.txt`.  The finite boundary classes
behave as follows:

```text
(a,b)=(0,1): fbar=(x+2)^2*(irreducible cubic), generalized_order=10
(a,b)=(1,1): fbar=x^2*(irreducible cubic),       generalized_order=14
(a,b)=(2,2): fbar=(x+2)^2*(irreducible cubic), generalized_order=10
```

The `(1,0)` class is more singular:

```text
(a,b)=(1,0): fbar=x*(x+2)^4
```

Thus the clean `(1,1)` branch is locally dead on the boundary stratum itself:
its limiting generalized Jacobian has order `14`, hence no 5-part.  The only
finite boundary strata with a visible 5-part on the limiting generalized
Jacobian are `(0,1)` and `(2,2)`, but in both of those the marked contact-7
point also degenerates modulo `3` (`h(1)=0`).

### Targeted height-30 boundary branches

The clean `(1,1)` branch:

```text
magma -b height:=30 only_tag:="finite:a=1,b=1,Q5=0" \
    progress_interval:=100000 max_exact:=200 \
    code/contact7_plus5_boundary3_search.m \
    > data/contact7_plus5_boundary3_a1b1_h30.txt
```

Result:

```text
checked 1234321
smooth 1233959
boundary 74528
survivors 0
exact_tests 0
hits 0
```

The two finite boundary branches with a 5-part in the limiting generalized
Jacobian were also tested to height `30`:

```text
magma -b height:=30 only_tag:="finite:a=0,b=1,h1=0,L=0" \
    progress_interval:=200000 max_exact:=200 \
    code/contact7_plus5_boundary3_search.m \
    > data/contact7_plus5_boundary3_a0b1_h30.txt

magma -b height:=30 only_tag:="finite:a=2,b=2,h1=0,L=0" \
    progress_interval:=200000 max_exact:=200 \
    code/contact7_plus5_boundary3_search.m \
    > data/contact7_plus5_boundary3_a2b2_h30.txt
```

Results:

```text
(a,b)=(0,1) branch: boundary 77166, survivors 0, hits 0
(a,b)=(2,2) branch: boundary 74442, survivors 0, hits 0
```

### Boundary conclusion

The `p=3` boundary analysis did not produce any simple `35` candidate.  More
importantly:

```text
(1,1) boundary: limiting generalized Jacobian has no 5-part.
(0,1),(2,2): limiting generalized Jacobian has a 5-part, but the marked
              contact-7 point degenerates and height-30 searches have no
              survivors away from p=3.
(1,0): highly singular x*(x+2)^4 boundary; not a promising first branch.
```

At this point the direct contact-7 plus 5 route looks much less promising than
it did before the boundary check.  Any remaining chance would require a more
careful blow-up analysis of the degenerate `h(1)=0` branches `(0,1)` and
`(2,2)`, not another broad height increase.
