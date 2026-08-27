# The bad and nonintegral `3`-adic boundary of the cyclic `[49]` incidence

## Outcome

The two open incidence points modulo `3` admit a clean weighted blowup. Its
exceptional divisor is a single affine line on each sign branch, and the
pullback of the ordinary-node discriminant has ramification index `7` along
that line. Exact modular Hensel calculations resolve six successive residue
layers, through discriminant depth `42`.

A one-dimensional rational reconstruction search in the resulting local
parameter found no rational point. At local height `100` it performed
`18,206` branch lifts; none of the other four coordinates reconstructed with
height at most `5000`, hence there was no exact open hit.

The balanced first-pole chart was also blown up three times. It remains a
nonreduced, higher-dimensional initial divisor, so it does not yet justify a
one-dimensional search. No five-variable box was run.

Reproducible artifacts are

```text
code/z49_incidence_3boundary.py
data/z49_incidence_3boundary_h100.txt
data/z49_incidence_3boundary_h100_time.txt
```

## The two open special points

Eliminating the four coefficients of `A` from

```text
A(x)^2-B(x)^2*f(x) = (x-1)*(x-r)^7,
B=u+v*x,
```

leaves four integral equations `E0,...,E3` in `(a,b,u,v,r)`. The only open
incidence points over `F_3` are

```text
(a,b,u,v,r)=(1,1,epsilon,0,0),  epsilon=+1,-1.
```

Their Jacobian matrices modulo `3` are

```text
epsilon=+1:                 epsilon=-1:
[1 0 0 0 0]                [1 0 0 0 0]
[1 0 1 1 0]                [1 0 2 2 0]
[2 0 2 1 1]                [2 0 1 2 1]
[1 0 2 2 1]                [1 0 1 1 1].
```

Thus the `(a,u,v,r)` block is invertible, the `b` column is zero, and the
tangent direction is exactly `db`.

## First weighted blowup

Use

```text
b = 1 + 3*T,
a = 1 + 9*A,
u = epsilon + 9*U,
v = 9*V,
r = 9*R.
```

Every `Ei` is divisible by `9`. After division and reduction modulo `3`,

```text
e0 =  A + T^2 + T + 1,
e1 =  A + T^2 + epsilon*U + epsilon*V - 1,
e2 = -A + R - T^2 - epsilon*U + epsilon*V + 1,
e3 =  A + R + T^2 - T - epsilon*U - epsilon*V + 1.
```

They have exactly the three `F_3` points on the exceptional line

```text
A = -(T^2+T+1),
U = epsilon*T,
V = -epsilon,
R = -1.
```

This is only an exceptional-divisor equation. Substitution of these first
approximations as rational formulas is not an exact component over `Q`;
higher Hensel corrections are essential.

## Sevenfold discriminant pullback

The family discriminant is

```text
256*Disc(f) = H^7*Q5(a,b),       H=2*a+2*b-5.
```

Here `H` is a `3`-adic unit, while

```text
Q5(1+9*A,1+3*T)/9 = -A-T^2-T-1  (mod 3),
```

the negative of `e0`. Thus the incidence exceptional divisor forces extra
cancellation. Exact modular Hensel lifting gives the same equations on both
sign branches:

```text
T = W:            Q5/3^7  =  1-W       (mod 3),
T = 1+3*W:        Q5/3^14 = -1-W       (mod 3),
T = 7+9*W:        Q5/3^21 =  1-W       (mod 3),
T = 16+27*W:      Q5/3^28 = -1-W       (mod 3),
T = 70+81*W:      Q5/3^35 =  1-W       (mod 3),
T = 151+243*W:    Q5/3^42 =  1-W       (mod 3).
```

At each stage exactly one residue proceeds to the next layer. The unresolved
deeper disc after six stages is

```text
T = 394 (mod 729).
```

Away from that nested discriminant point the observed node thickness is
`7,14,21,...`, as required by the Neron-model condition `7 | v3(Q5)` for
the split node at `(a,b)=(1,1)`. The computation certifies the six displayed
finite layers; it is not an all-depth analytic identity.

## One-dimensional reconstruction search

The search enumerates every reduced

```text
T=m/n, |m|<=100, 1<=n<=100, gcd(m,n)=1, 3 does not divide n,
```

sets `b=1+3*T`, and lifts `(a,u,v,r)` on both branches to `3^18`. Each is
rationally reconstructed with numerator and denominator at most `5000`.
Since `3^18 > 2*5000^2`, reconstruction is unique in that box. Any candidate
is substituted into all four exact equations and the exact open conditions
`r!=1`, `H!=0`, `B(1)!=0`, `B(r)!=0`, and `Q5!=0`.

```text
rational T values                    9,103
sign branches                            2
Hensel lifts                         18,206
v3(Q5)=7                             12,176
v3(Q5)=14                             4,018
v3(Q5)>=18                            2,012
all four coordinates reconstructed       0
exact open hits                           0
```

The last valuation row combines depth `21` and deeper discs because the
search precision is `18`; the separate layer calculation resolves their
first six exceptional equations. This is not redundant with the earlier
height-`200` search in `b`: it bounds the natural coordinate `T=(b-1)/3`,
so some tested `b` values lie outside that naive-height box.

## Balanced first-pole chart

Consider only

```text
v3(b)=v3(v)=-1,  v3(a),v3(u),v3(r)>=0,
a=A, b=B/3, u=U, v=V/3, r=R,
B,V units.
```

The four lowest initial forms have valuations `-16,-14,-12,-10`. The last
is `V^8*(V^2-B^2)`, and the others are compatible with it, forcing

```text
V=delta*B,  delta=+1,-1.
```

Put `v=delta*B/3+W`. The next fourth initial equation is
`B^9*(-B-delta*W)=0`, hence `W=-delta*B`. Now write

```text
v=delta*B/3-delta*B+3*X.
```

All four third-stage initial forms are divisible over `F_3` by

```text
delta*B*X-R+1.
```

After imposing it, `A,U,R` remain free. Thus the chart is still thick and
nonreduced and has not become the one-dimensional normalization. Searching
`(A,U,R)` would be another blind box, so none was performed.

This covers only the balanced first pole with `a,u,r` integral. Asymmetric
poles, depths `v3(b)<-1`, and charts where `a,u,r` also have poles remain
unresolved.

## Resources and conclusion

The full height-`100` run took `54.37` seconds and used `85,296 KB` maximum
resident memory under a hard `200 MB` address-space cap. The derivation-only
run used `85,292 KB`. No large elimination or multidimensional box was used.

The bad integral boundary is now a controlled pair of one-parameter branches,
with visibly sevenfold discriminant pullback through depth `42` and no
rational reconstruction in the stated natural box. This is meaningful
negative evidence, not a proof that the incidence curve has no rational
point. The main unresolved local issue is normalization of the nonintegral
charts, beginning with the thick balanced pole divisor above.
