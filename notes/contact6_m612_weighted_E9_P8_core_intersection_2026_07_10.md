# The rational `P8` endpoint family against the `b=0` contact core

Date: 2026-07-10.

## Result

The exact rational `P8` family of `R3` halves was searched against the
independent cubic-contact `[3,3]` core through parameter height `30`.
No rational contact-open core point was found.

The completed exact run gave

```text
rational tau values enumerated                    1111
distinct nonzero e fibers tested                   858
duplicate e fibers                                 251
finite parameter poles                               0
e=0 parameter-boundary values                        2
dimension > 0 exceptional core fibers                0
unit-ideal (dimension -1) open core fibers            0
saturation failures                                  0
Variety failures                                     0
rational points on the saturated M=L^2 quotient      0
rational-square M lifts                               0
verified independent [3,3] cores                      0
dual torsion audits / [6,12] hits                    0 / 0
runtime                                      239.94 seconds
```

This is a bounded parameter-height result on the rational `P8`
normalization, not a proof that its fiber product with the core has no
rational point.


## The `P8` parameter

Write the normalization parameter as `tau`.  The conic coordinate is

```text
t = 4*(tau^2+tau-6)/(tau^2+6),
y = 2*(tau^2-24*tau-6)/(tau^2+6),
y^2 = 100-6*t^2.
```

The exact functions `e(tau), mu(tau), nu(tau)` are generated and verified in

```text
code/contact6_m612_weighted_R3_geometry_p8_param.m.
```

The scan enumerates every reduced rational `tau=n/d` with

```text
|n| <= 30,  1 <= d <= 30,
```

and also audits the projective value `tau=infinity`.  Fibers are deduplicated
by `e`, since the contact core depends only on `e`.

The two finite boundary parameters are

```text
tau=-3, 2,
```

where `(e,mu,nu)=(0,0,0)` and the square-quartic chart is not open.  There
were no finite denominator poles in the height-30 set.

At projective infinity,

```text
(e,mu,nu)=(-200/409,
           -36320/167281,
           112592/501843).
```

Its `e`-fiber is the same as the finite `P8` point at `tau=12`, so it was
already tested after deduplication.  The latter point is

```text
(e,mu,nu)=(-200/409,
           -36320/167281,
           38136/167281).
```


## Exact fixed-core test

For every nonzero `e`, put `a=1/e`, `b=0` and solve the exact three-equation
core in `(M,U,v)`, saturated by

```text
M*v*(U^2-4*v^2).
```

Every rational point would then have to pass

1. `M=L^2` over `Q`;
2. the exact cubic-contact identity;
3. orders `6` and `3` for the two marked classes;
4. independence of those classes;
5. exact dual torsion and direct `R3` divisibility on both Richelot twists.

The height-30 run did not reach step 1: every saturated fiber had an empty
rational `Variety`.

The known `P8` fiber `e=-200/409` illustrates why the saturation is
essential.  Before removing the repeated-root boundary, its only rational
point with `M != 0` is

```text
(M,U,v)=(1,-2,1).
```

Indeed this is the universal exact section

```text
M=1, U=-2, v=1, N=0, R=1/e,
```

valid for every `e`.  But

```text
x^2+U*x+v^2=(x-1)^2,
U^2-4*v^2=0,
```

so it is not an independent rational `3`-direction.  The contact-open
saturation removes it, and the known `P8` fiber is then rationally empty.


## Generic-degree warning

The later function-field calculation works over `Q(e)`, quotients
`L -> -L` by writing `M=L^2`, eliminates the two linear variables `N,R`,
and removes the obvious component

```text
M=0, U=2*v.
```

The resulting generic support fiber has length `40` and exact factorization

```text
1 + 12 + 27.
```

The degree-`1` factor is the universal repeated-support boundary, while the
degree-`12` orthogonal factor and degree-`27` nonorthogonal factor are
contact-open.  Thus there is no hidden degree-`1` contact-open component.
The rational-point problem still cannot treat `e` alone as a parameter: the
relevant support cover has degree `12`, and its signed lift has degree `24`.


## Reproduction

From `torsion_jac`:

```text
magma -b height:=30 progress_interval:=200 \
  code/contact6_m612_weighted_E9_P8_core_intersection.m

magma -b code/contact6_m612_weighted_E9_core_generic_degree.m
magma -b code/contact6_m612_relative3_exact_reconstruct.m
```

The first script classifies Magma dimension `-1` as an exact empty ideal,
not as an exceptional or failed fiber.  Only positive-dimensional fibers or
actual saturation/`Variety` failures are reported as exceptional.


## Next step

The rational `P8` lane remains globally possible but did not meet the
independent core in the tested range.  The proof-level continuation has since
isolated the relevant characteristic-zero component: an irreducible
degree-`12` orthogonal support factor and its connected irreducible
degree-`24` signed lift.  The height-`1000` modular sieve has no nonboundary
survivor, while `p=7` and `p=17` both have smooth target branches in some
forced disks.  The signed curve over the `e`-line has exact genus `10`, so
its rational points, Jacobian, and `S3` quotients are the first global target;
the exact P8 support/signed genera are `73/145`, and their global points and
the remaining forced primes `19,23,41` are the later tasks, not another
undirected height increase.  See
`notes/contact6_m612_p8_relative3_status_2026_07_11.md` for current caveats.
