# Contact-7 two-root surface plus 3: geometry stop gate

## Two-root surface

Put

```text
A(u)=(2u^5+4u^4+6u^3+8u^2+10u+5)/(2(u+1)^2).
```

The signed parameters absorb the two sign choices.  With

```text
r=1-s^2,  w=1-t^2,
h(r)=s^7, h(w)=t^7,
```

the contact-7 coefficients are

```text
b=(A(t)-A(s))/(s^2-t^2),
a=A(s)-b(1-s^2).
```

Thus the two-root locus is a rational surface.  Translating the first root
to zero makes the cubic-contact constant equation `E^2=V^3`.  On its open
normalization use `V=j^2,E=j^3`; alternatively retain `(V,E)` and eliminate
high-to-low.  For

```text
q=x^2+U*x+V,
H=x^3+A2*x^2+B2*x+E,
H^2-q^3=M*f,
```

the coefficients of `x^5,x^4,x^3` reconstruct `A2,B2,E`, leaving three
equations.  The exact rational 3-torsion cover further adjoins `L^2=M`.

The generic surface equations already have total degrees `19,20,30` and
`328,385,462` terms in the smaller direct reconstruction.  Fixed `s=2`
equations have degrees `12,12,16`; their direct four-variable elimination
is substantially larger than the gatekeeper below and was stopped.

## The removable `s=1` gatekeeper

The apparent `s=1` pole is removable and must not be discarded.  Here

```text
a=35/8,
b=-(8t^4+24t^3+48t^2+45t+15)/(8(t+1)^3).
```

The resulting quintic has roots `0` and `1-t^2`.  Up to a nonzero rational
constant, its discriminant is

```text
t^14 (t-1)^2
*(t^2+9t/8+3/8)^7
*(t^3+4t^2+3t+3/4)^2
*(t^3+29t^2/16+5t/4+5/16)^2
*S6(t)
/(t+1)^36,
```

where

```text
S6=t^6+55t^5/27+115t^4/108-83t^3/288
   -83t^2/288+5t/96+5/96
```

is irreducible.  Its six roots are therefore six distinct simple
discriminant zeros.  Since the quintic is monic, each is a one-node
Picard--Lefschetz degeneration.

Exact elimination over `F_5(t)`, after saturating the trivial `M=0`
component, gives a triangular basis with an irreducible degree-40
polynomial in `M`.  Adjoining `L^2=M` is also irreducible.  This certifies
connected degree-40 and degree-80 covers at the good prime 5.  The equations
have denominators supported only at 2 and `t+1`, the generic quintic remains
smooth modulo 5, and the degrees do not drop.  Consequently the usual
primitive-polynomial/good-reduction argument lifts irreducibility to
`Q(t)`: a characteristic-zero factorization would reduce to a nontrivial
factorization of the displayed modular elimination polynomial.

The exact transcript is in
`data/contact7_two_root_plus3_geometry_s1_f5.txt`.

## Rigorous genus lower bounds

A transvection on the 40 lines of `F_3^4` fixes the 13 lines contained in
its orthogonal hyperplane and permutes the other 27 lines in nine 3-cycles.
Its Riemann--Hurwitz contribution is therefore `9*(3-1)=18`.

The six simple nodal fibers alone give

```text
2g_40-2 >= -2*40 + 6*18 = 28,
g_40 >= 15.
```

On the 80 nonzero vectors, the same transvection fixes 26 vectors and has
18 three-cycles, contributing 36.  Hence

```text
2g_80-2 >= -2*80 + 6*36 = 56,
g_80 >= 29.
```

All omitted boundary fibers can only increase these genera.  In particular,
the only natural quotient of the exact cover, the sign quotient, already
has genus at least 15; it is not genus 0, 1, or 2.

## Stop verdict

This route fails the agreed stop rule.  Even its simplest removable
one-parameter branch has a connected exact 3-torsion cover of genus at
least 29, and its natural quotient has genus at least 15.  There is no
natural elliptic or genus-2 quotient exposed by the contact construction.
Further generic Groebner elimination or normalization is not justified for
the record search.

This is a stop verdict for the natural one-parameter gatekeeper, not a proof
that the entire two-parameter surface has no rational point on its 3-torsion
cover.  On the full surface the result remains the complete bounded
height-100 search recorded separately.

Code:

```text
code/contact7_two_root_plus3_geometry.m
code/contact7_two_root_plus3_geometry_s1.m
code/contact7_two_root_plus3_geometry_s1_functionfield.m
```
