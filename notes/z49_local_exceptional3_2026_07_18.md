# Exceptional `p=3` charts for cyclic `[49]` in the contact-7 family

This note treats the two mod-`3` discriminant residues not covered by the
ordinary `h(1)=0` analysis, plus a small parameter-pole probe.  The family is

```text
h = 1-(7/2)x+a*x^2+b*x^3,
f = (h^2+(x-1)^7)/x^2,
D7 = [(x-1),h(1)].
```

Reproducible files:

```text
code/z49_local_q5_special_fiber.m
data/z49_local_q5_special_fiber.txt
code/z49_local_exceptional3_charts.m
data/z49_local_exceptional3_h12_p43.txt
data/z49_local_exceptional3_h12_p43_time.txt
```

Every away-from-`3` filter uses the exact necessary condition

```text
D7_bar in 7*J(F_p),
```

computed in the invariant-factor coordinates of `AbelianGroup(J(F_p))`.
It does not use mere divisibility of `#J(F_p)` or of the group exponent.
The primes are

```text
5,11,13,17,19,23,29,31,37,41,43.
```

## The four singular affine residues modulo 3

The discriminant identity is

```text
256*Disc(f) = H^7*Q5(a,b),       H=2*a+2*b-5.
```

The ordinary points `(0,1)` and `(2,2)` lie only on `H=0` and are handled
elsewhere.  The remaining residues are

```text
(1,1): Q5=0, H!=0;
(1,0): Q5=H=0.
```

At `(1,1)`, `dQ5/da=2 mod 3`, so the `Q5` component is smooth.  At `(1,0)`
both derivatives vanish.  If `a=1+A`, `b=B`, its tangent cone is

```text
Q5 = -(A^2-A*B+B^2) + terms of degree >=3
   = -(A+B)^2              mod 3,
```

which is the doubled tangent of `H=0`.  Thus the intersection is genuinely
nontransverse.

## The `Q5`-only ordinary node at `(1,1)`

The exact special fiber is

```text
f = x^2*(x^3+2*x+1).
```

It has a split ordinary node at `x=0`; its two normalization points are
`(0,+1)` and `(0,-1)`.  The normalization is

```text
E: Y^2=x^3+2*x+1,
E(F_3) = [7].
```

Since `y=xY`, the marked point maps to `(1,1) in E(F_3)`, which has exact
order `7`.  The split torus has order `2`, so the Neron identity group has
order `14`.  If the node thickness `n` is prime to `7`, its component group
adds no 7-part and the marked generator cannot lie in seven times the special
fiber group.  A necessary condition is therefore

```text
7 | n.
```

Because `H` is a unit and `Q5` is a smooth local equation here,

```text
n = v_3(Q5(a,b)).
```

The first possible layer is `v_3(Q5)=7`.  For every small 3-integral
`b=1 mod 3`, the script Hensel-lifts the unique root `a0(b)` modulo `3^7`
and writes

```text
a = a0(b) + 3^7*u.
```

With local-coordinate height `12`, the search had

```text
b rows                              43
u rows                             135
pairs checked                    5,805
smooth points with v3(Q5)=7      3,878
away-prime survivors                 0
```

The first-kill counts were

```text
p=5: 1602, p=11: 1312, p=13: 721, p=17: 188,
p=19: 43, p=23: 8, p=29: 4.
```

Thus the first locally possible thickness layer contains no candidate in
this bounded Hensel-coordinate box.

## Tangential intersection at `(1,0)`

Here

```text
f = x*(x-1)^4,
```

so this is a nonordinary four-root cluster, not a one-node fiber.  Put
`H=2*a+2*b-5` and `z=x-1`.  The contact polynomial expands exactly as

```text
h(1+z) = H/2 + (H+b+3/2)z
         + (H/2+5/2+2b)z^2 + b*z^3.
```

At the residue, the constant and linear coefficients vanish while the
quadratic coefficient is a unit.  Consequently the ordinary-node component
argument does not determine the full component group.

As a bounded first probe, the script tests the natural contact-depth slice

```text
a = 1+3*u,
b = 5/2-a+3^7*s,
h(1)=3^7*s,             v3(s)=0.
```

This is the analogue of the first potentially 7-divisible contact layer, but
it is not claimed to exhaust every relative valuation in the four-root
cluster.  At local-coordinate height `12`:

```text
pairs checked        11,610
smooth               11,610
away-prime survivors      0
```

The final two points were killed at `p=31`; no exact rational division test
was needed.

## Small nonintegral parameter charts

For a common parameter pole write

```text
a=u/3^e, b=v/3^e,
min(v3(u),v3(v))=0.
```

After multiplying the equation by `3^(2e)`, the first special polynomial is

```text
x^2*(u+v*x)^2,
```

so this chart is also highly nonordinary and the first scaling alone gives no
component-group sufficiency statement.  The bounded probe used `e=1,2,3`
and 3-integral local coordinates of height `12`:

```text
checked               47,472
smooth                47,458
away-prime survivors       0
```

Again all candidates were removed by exact marked finite-group masks; the
last seven were killed at `p=31`.

## Resources and conclusion

The combined Magma run completed with exit status zero in `26.03` seconds and
used `82,688 KB` maximum resident memory under a hard `200 MB` cap.  The
separate special-fiber certificate also ran under `200 MB`.

The strongest conclusion is on the `Q5`-only branch: its ordinary-node
geometry rigorously forces `7 | v3(Q5)`, and the first possible layer
`v3(Q5)=7` has no candidate in the searched local box.  The intersection and
infinity computations are useful negative bounded probes, not complete local
resolutions.  In particular, deeper `Q5` thicknesses `14,21,...`, other
relative valuations at the four-root intersection, and deeper/asymmetric
parameter-pole charts remain outside the claim.
