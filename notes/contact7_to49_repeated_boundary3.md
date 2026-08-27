# The repeated-root locus on the forced 3-adic boundary

This continues the specialization

```
Q = 2P-K,                 U=(x-p)^2,
14*(P-infinity) = +/-D7.
```

No smooth contact-7 fiber over `F_3` contains a point satisfying the
displayed relation.  A rational solution must therefore meet a bad
curve fiber or have nonintegral parameters/section at 3.

## Ordinary `h(1)=0` branches

Write `H=2*a+2*b-5=2*h(1)`.  At the residues `(0,1)` and `(2,2)` the
special fiber has one ordinary node.  Its component group has order

```
n = 7*v3(h(1)).
```

The marked class has order 7 in this cyclic component group.  If a
generator is oriented so that `D7=7*g` when `n=49`, the repeated-point
condition gives

```
14*t = +/-7 (mod 49), hence t = +/-4 (mod 7).
```

More generally the component equation can be solvable only if
`49 | n`; equivalently `7 | v3(h(1))`.  The first feasible blow-up is

```
a = a0+3*u,
b = 5/2-a+3^7*s,       v3(s)=0,
a0 in {0,2}.
```

The congruence on `t` is an additional necessary condition on the
section component.  The bounded computation conservatively retained
all section components and imposed the exact repeated-point relation
at every usable away prime.

## The `Q5=0` branch

At `(a,b)=(1,1)` the special fiber has a split ordinary node at zero.
Its normalization has `E(F_3)=Z/7`, and the identity group has order
14.  When the node thickness is prime to 7, multiplication by 14
kills the only 7-primary subgroup and cannot produce the marked
generator.  Thus

```
7 | v3(Q5(a,b))
```

is necessary.  Since `Q5=0` is smooth at this residue, the first chart
is obtained by Hensel lifting its root and writing

```
a = a0(b)+3^7*u,       b=1 (mod 3),
v3(Q5(a,b))=7.
```

This is necessary, not a proof of local sufficiency: the extension of
the identity and component groups can impose a further condition on
the section.

## Intersection and pole charts

At `(1,0)`, `f=x*(x-1)^4`; the `H` and `Q5` components are tangent and
the fiber is a four-root cluster.  The ordinary-node calculation does
not resolve its component group.  The computation tests the natural
first contact-depth slice

```
a=1+3*u, b=5/2-a+3^7*s, v3(s)=0,
```

but does not claim that this exhausts all relative valuations.

For a common parameter pole `a=u/3^e, b=v/3^e`, the first scaled
special equation is `x^2*(u+v*x)^2`, again nonordinary.  The run tests
`e=1`.  For the point section it tests both 3-integral `x(P)` and the
first pole layer `x(P)=w/3`.

## Exact bounded search

For every curve in these charts and every tested abscissa, the script
checks whether `f(x(P))` is an exact rational square.  Each resulting
rational point is then screened by the exact condition

```
14*[P-infinity] = +/-D7 in J(F_p)
```

for all usable primes through 43.  This is stronger than testing
`D7 in 7*J(F_p)`.

At local-coordinate height 6 the run examined 3,622 smooth curves,
206,454 abscissas and 7,656 rational points.  No point survived the
finite masks, so no rational exact survivor remained to check.  Peak
memory was 24 MB.

The result excludes this bounded collection of first boundary charts,
not the whole expected repeated-root curve.  Deeper contact or Q5
thicknesses, a full resolution of the `(1,0)` cluster, asymmetric
parameter poles, and deeper section-pole layers remain open.

Files:

```
code/contact7_to49_repeated_boundary3.m
data/contact7_to49_repeated_boundary3_h6.txt
```
