# An explicit exact cyclic `[60]` HLP control

The Howe--Leprévost--Poonen equation (3), specialized at `t=1/3`, is

```text
(-5/27)y^2 = (2u-1)(4u^2-2u-1).
```

It has the small nondegenerate rational point `(u,y)=(-1,9)`.  Applying HLP
Proposition 4 to the corresponding universal elliptic curves with rational
10- and 12-torsion gives the genus-2 curve

```text
C: y^2 = -185*(125*x^2-1728)
              *(2000*x^4-48525*x^2+294912).
```

A smaller reduced Weierstrass model is

```text
y^2 + x^2*y = -106560*x^6 + 704600*x^4 - 1550855*x^2 + 1136640.
```

Magma certifies

```text
J(C)(Q)_tors = [60].
```

One generator in Mumford coordinates is

```text
(x^2 + 74/15*x + 128/25,
 -1081510/9*x - 1226624/3,
 2),
```

and has exact order `60`.

This is an exact cyclic-`[60]` positive control, but it is **not** the desired
geometrically simple example.  The displayed sextic is even, so
`(x,y) -> (-x,y)` is a non-hyperelliptic involution and gives elliptic
quotients; equivalently, by the HLP construction its Jacobian is
`(2,2)`-isogenous to the product of the two input elliptic curves.  The simple
search should therefore continue on the full quadratic-`B` and unrestricted
Hermite loci.

Reproduce with `code/hlp_z60_explicit_verify.m`.
