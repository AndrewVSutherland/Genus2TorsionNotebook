# Z/35 contact-7: first blow-up at the degenerate `h(1)=0` branches

This note continues the `contact-7 plus 5` route from
`simple_35_attempt.md`, but does not repeat the broad boundary-at-3 height
search.  I only looked locally at the visible `p=3` boundary component

```text
2*a + 2*b - 5 = 0,
```

at the degenerate finite residue branches `(a,b)=(0,1)` and `(2,2)`.  These
are the two branches where the naive limiting generalized Jacobian has a
5-part, but the marked contact-7 point reduces to the node.

## Local chart

The contact-7 family is

```text
h = 1 - (7/2)*x + a*x^2 + b*x^3,
f = (h^2 + (x - 1)^7)/x^2.
```

Since

```text
h(1) = a + b - 5/2,
```

the visible boundary is `h(1)=0`.  For the two degenerate branches use the
3-adic chart

```text
a = a0 + 3*u,
b = 5/2 - a + 3*s,
h(1) = 3*s,
```

where `a0=0` gives `(a,b)=(0,1) mod 3`, and `a0=2` gives `(a,b)=(2,2) mod 3`.
The first blow-up layer is `s` a 3-adic unit.

Near the node put `z=x-1` and `Y=x*y`.  Then

```text
Y^2 = H(z)^2 + z^7,       H(z) = h(1+z),
```

and with `R=Y-H(z)`:

```text
R*(R + 2*H(z)) = z^7.
```

On the chart,

```text
H(z) = 3*s + lambda*z + O(z^2),
lambda = 4 - a0 - 3*u + 9*s,
```

and `lambda` is a 3-adic unit on both branches.  Also

```text
Q5(0,1) = 2 mod 3,
Q5(2,2) = 1 mod 3,
```

so away from the visible factor the discriminant is a 3-adic unit.  Thus

```text
Disc(f) = (2*a + 2*b - 5)^7 * Q5(a,b) / 256
        = (6*s)^7 * unit / 256.
```

For `v3(s)=r >= 0`, this gives local node thickness

```text
n = 7*(r+1).
```

In particular, the first blow-up layer has thickness `7`.  The contact divisor
does not disappear: the function `Y-H` still has the contact relation
`(Y-H)*(Y+H)=z^7`, and the order-7 class moves into the thickness-7 component
group.  The probe script also checks sample fibers directly in Magma; all
tested samples had exact order `7`.

## Probe

Code:

```text
code/agent_Z35_contact7_boundary_probe.m
```

Primary branch, first layer:

```text
magma -b branch_a0:=0 height:=18 first_layer_only:=1 \
    progress_interval:=20000 max_contact_checks:=5 \
    code/agent_Z35_contact7_boundary_probe.m
```

Summary:

```text
branch <0, 1>
checked 59994
smooth 59994
survivors 0
contact_checks 5
contact_failures 0
LAYER_COUNTS
0 59994
FIRST_KILL_5_AWAY_FROM_3
7 24846
11 19757
13 8329
17 4137
19 2250
23 506
29 117
31 34
37 15
41 3
```

Sibling branch, first layer:

```text
magma -b branch_a0:=2 height:=16 first_layer_only:=1 \
    progress_interval:=20000 max_contact_checks:=3 \
    code/agent_Z35_contact7_boundary_probe.m
```

Summary:

```text
branch <2, 2>
checked 36498
smooth 36498
survivors 0
contact_checks 3
contact_failures 0
LAYER_COUNTS
0 36498
FIRST_KILL_5_AWAY_FROM_3
7 14145
11 11867
13 5556
17 3542
19 1080
23 218
29 60
31 21
37 9
```

I also allowed one deeper local layer on the primary branch:

```text
magma -b branch_a0:=0 height:=14 first_layer_only:=0 \
    max_layer:=1 progress_interval:=20000 max_contact_checks:=3 \
    code/agent_Z35_contact7_boundary_probe.m
```

Summary:

```text
branch <0, 1>
checked 33930
smooth 33930
survivors 0
contact_checks 3
contact_failures 0
LAYER_COUNTS
0 26130
1 7800
FIRST_KILL_5_AWAY_FROM_3
7 13349
11 11095
13 5057
17 3224
19 936
23 201
29 45
31 11
37 11
41 1
```

## Outcome

The first blow-up does **not** kill the contact-7 class.  On the first layer
`v3(s)=0`, the singularity has thickness `7`, and the rational contact class
is naturally accounted for in the component group.  The direct Magma checks of
sample fibers confirm exact order `7`.

However, the branch looks dead for the `35` search at low local order.  The
necessary condition `5 | #J(F_p)` for good primes `p != 3,5` killed every
tested local lift on `(0,1)`, on `(2,2)`, and on the one-step-deeper
`v3(s)=1` primary chart.  Combined with the existing height-30 broad boundary
search, I would not treat these `h(1)=0` branches as genuinely viable unless a
later argument specifically forces very high 3-adic depth or very large
blow-up coordinates.
