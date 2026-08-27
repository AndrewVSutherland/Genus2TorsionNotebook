# Geometry of the first orbit-12 radicand

This note studies the first square condition on the orbit-12 cover of the
Elkies Clebsch--Klein family.  The fixed marked class is
`{r1^2,r2^2}`, and

```text
G0 = -(r1^2-r3^2)(r1^2-r4^2)(r1^2-r5^2).
```

The exact calculations are certified by

```text
code/elkies22210_orbit12_first_radicand_geometry.m.
```

The main outcome is slightly subtler than the hoped-for elliptic
fibration.  The raw first-radicand cover has genus-5 fibers in either
coordinate of the complete rational CK chart.  It has a very useful
elliptic **S3 quotient**, but that elliptic surface has no non-boundary
rational section.  Moreover, the one rational genus-drop fiber selected
by the branch calculation, `q=-2`, has no smooth rational point at all.

## 1. Pullback to the complete CK chart

Use the complete labelled chart

```text
R1 = 1+t(t+2)m,
R2 = tm(m-t-2),
R3 = -1+m+t(t+1)m^2,
R4 = 1+t-m-tm^2,
R5 = -(1+t)(1+tm^2).
```

It satisfies `sum Ri=sum Ri^3=0`.  Define

```text
A = t^2*m^2-t^2*m+t*m^2-2*t*m+m-2,
B = t^2*m+t*m^2+2*t*m-t+m,
C = t^2*m^2+t^2*m+t*m^2+2*t*m+t+2.
```

Exact factorization gives

```text
G0(t,m) = -m(m-1)t(t+1)(t-m+2)(tm+1)
            *(tm+m-1)(tm+t+1)*A*B*C.
```

It is squarefree of degree `12` in `m` over `Q(t)` and degree `12` in
`t` over `Q(m)`.  Thus both coordinate projections give generic
hyperelliptic fibers of genus `5`.  In particular, merely rewriting the
complete chart does not produce a genus-one curve.

## 2. The symmetric quotient

On `r1 != 0`, scale `r1=1` and put

```text
q = r2/r1,
s = (r3*r4+r3*r5+r4*r5)/r1^2.
```

If `x3,x4,x5` are the normalized complementary roots, the two CK
equations give

```text
e1(x3,x4,x5) = -(1+q),
e3(x3,x4,x5) = (1+q)(q-s).
```

Consequently the complementary roots are the roots of

```text
p_q,s(U) = U^3+(1+q)U^2+sU-(1+q)(q-s),
```

and a direct symmetric-function calculation gives the squareclass
identity

```text
G0/r1^6 = q(s-q)(2-q^2+(q+2)s).
```

Since `r1^6` is a square, the quotient of the first-radicand cover by
permutations of `r3,r4,r5` is

```text
Ebar: v^2 = q(s-q)(2-q^2+(q+2)s).
```

For fixed `s` this is a genus-one quartic.  A point on the original
labelled cover maps to a point of `Ebar`; lifting back additionally
requires `p_q,s` to split completely over `Q`.  This splitting condition
must not be dropped: the elliptic surface is a degree-six quotient of the
labelled problem, not a birational model of it.

## 3. Weierstrass model and all rational sections

Put

```text
a = 2s(s+1),
X = a/q,
Y = a*v/q^2.
```

Then the quotient has the Weierstrass equation

```text
E_s: Y^2 = X^3+(s^2-2s-2)X^2-4s^2(s+1)X+4s^2(s+1)^2
           = (X-2(s+1))
             *(X^2+s^2*X-2s^2(s+1)).
```

Its discriminant and `c4` are

```text
Delta = 256*s^2*(s+1)^4*(s^2+8s+8),
c4    = 16*(s^4+8s^3+12s^2+8s+4).
```

The geometric singular fibers are

```text
I2 at s=0,
I4 at s=-1,
I4 at s=infinity,
I1 at each root of s^2+8s+8.
```

This is a rational elliptic surface.  Shioda--Tate gives

```text
rank E(Qbar(s)) = 10-2-(1+3+3) = 1.
```

The quartic has the geometric non-torsion section

```text
q=-1,  v=i(s+1),
```

or, on the Weierstrass model,

```text
(X,Y)=(-2s(s+1), 2i*s*(s+1)^2).
```

Its specialization at `s=2` has infinite order over `Q(i)`, as certified
by Magma.  Complex conjugation sends this section to its negative.  Since
the geometric free rank is one, the invariant free part is zero:

```text
rank E(Q(s)) = 0.
```

There is a visible point of order four

```text
P4 = (0,2s(s+1)),
2P4 = (2(s+1),0).
```

Generic torsion injects into the good specialization `s=1`, whose exact
torsion subgroup is `Z/4`.  Hence

```text
E(Q(s)) = Z/4.
```

In quartic coordinates these four sections are `q=0`, `q=s`, and the two
points with `q=infinity`.  They all lie on the CK/discriminant boundary.
Therefore:

```text
the elliptic quotient has no non-boundary rational section over Q(s).
```

This completely exhausts the proposed search for a generic rational
section.  A record point would have to arise from a genuine multisection
or from an isolated rational point on a specialized fiber.

## 4. The rational genus-drop fiber q=-2

Up to a nonzero constant, the branch resultant in the fixed-`q` model is

```text
q^11*(q-1)^2*(q+1)^4*(q+2)
    *(q^2-2q-4)*(q^2+3q+8/3)^2.
```

The values `q=0,+-1` are boundary, and both displayed quadratic factors
are irreducible over `Q`.  Thus `q=-2` is the only non-boundary rational
value at which the relevant branch configuration degenerates.  Set

```text
r1=1, r2=-2, r3=x.
```

After eliminating `r4,r5`, the first-radicand square condition becomes

```text
C0: u^2 = -x(x-1)(x+1)(x-2).
```

The condition that the remaining quadratic split over `Q` is

```text
Csplit: w^2 = (x-1)(x^3+x^2-x-9).
```

Both quotients have genus one.  Their third `V4` quotient is

```text
z^2 = -x(x+1)(x-2)(x^3+x^2-x-9),
```

of genus two, so the simultaneous fiber product has genus
`1+1+2=4`.

In fact the first equation alone kills the fiber.  Magma maps `C0` to

```text
E_-2: y^2 = x^3+x^2-4x-4
             = (x+1)(x-2)(x+2).
```

It certifies

```text
rank E_-2(Q) = 0,
E_-2(Q)_tors = (Z/2)^2.
```

Pulling back the complete Mordell--Weil group gives exactly

```text
x = -1, 0, 1, 2.
```

All four values make the first-radicand quartic zero, hence lie on the
repeated-root CK boundary.  Thus there is no smooth rational point on the
`q=-2` first-radicand fiber, before imposing the split-completion equation
or any of the other three orbit-12 square conditions.

## 5. Consequence for the search

The first-radicand geometry gives two rigorous negative conclusions:

1. the elliptic quotient has no non-boundary `Q(s)`-section;
2. the only rational genus-drop fiber `q=-2` has no smooth rational point.

The useful surviving route is therefore not a larger blind `(t,m)` box and
not generic section enumeration.  It is a search for low-degree rational
multisections of the elliptic quotient, or fiberwise Mordell--Weil searches
at rational `s` constrained by the existing compatible CRT disks, followed
by exact splitting of `p_q,s` and the remaining three radicand tests.
