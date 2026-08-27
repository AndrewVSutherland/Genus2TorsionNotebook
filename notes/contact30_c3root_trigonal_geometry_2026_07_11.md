# Geometry of the order-30 `C3`-root cover

## Outcome

The genuine rational-root cover is not a low-genus trigonal curve.  Its
normalization has genus `12`.  It has one rational involution, and quotienting
by that involution gives a trigonal curve of genus `6`.  The genus-6 quotient
has no further nontrivial automorphism over `Q`.

More decisively for the proposed `[2,60]` Richelot construction, the
dual-halving squareclasses descend to the genus-6 quotient.  Their three
quadratic covers have genera

```text
squareclass       genus       relative ramification degree
disc(A0)            17                    12
disc(B0)            18                    14
disc(A0)disc(B0)     24                    26
```

The individual Richelot gates `U` and `V` use the genus-18 and genus-17 covers,
respectively.  The third gate `W` requires **both** discriminants to split, so
its actual cover is the `V4` fiber product, not the genus-24 product-square
quotient.  Its genus is

```text
17+18+24-2*6 = 47.
```

Thus even the weakest necessary correlated condition has genus `24`, while
the actual simultaneous gate has genus `47`.  This route fails the agreed
geometry stop rule.

The exact, assertion-driven computation is

```text
code/contact30_c3root_trigonal_geometry.m
```

## The degree-three function field

For either sign `eps=+/-1`, use the order-30 parametrization

```text
t = (5R^2-20R+19)/(R^2-5),
Y = -2(5R^2-22R+25)/(R^2-5),
u = t^3,
s = t^5+t^4+(5/2)t^3+(1/2)t
    + eps*t*(t-1/2)*(t+1)*Y.
```

Recover `q`, and hence `A,B,C`, as in the order-30 construction.  The root
cover is the cubic extension of `Q(R)` defined by

```text
2*rho^3 + (A-3)*rho^2 + (B+3)*rho + (C-1) = 0.
```

The cubic is irreducible over `Q(R)` for both signs.  The two sign models are
isomorphic under

```text
R -> (11R-25)/(5R-11),     rho -> rho.
```

It therefore suffices to discuss `eps=-1`.

## Exact ramification and genus 12

The raw polynomial discriminant contains large repeated factors caused by a
nonmaximal affine order.  Factoring the different of the normalized function
field removes these apparent branches.  For `eps=-1`, the actual ramification
is:

| base place | degree | ramification | different contribution |
|---|---:|---:|---:|
| `R=2` | 1 | `e=2` | 1 |
| `R=3` | 1 | `e=2` | 1 |
| `R=7/3` | 1 | `e=3` | 2 |
| `31R^4-220R^3+574R^2-660R+291=0` | 4 | `e=2` | 4 |
| irreducible `P20(R)=0` | 20 | `e=2` | 20 |

The degree-20 polynomial is printed exactly by the script.  Its leading
normalization begins

```text
P20 = R^20 -(24788818/675105)R^19
      +(1289032208/2025315)R^18 - ... .
```

The different has total degree

```text
1+1+2+4+20 = 28.
```

Riemann--Hurwitz for the degree-three map to `P1_R` gives

```text
2g-2 = -6+28 = 22,
g = 12.
```

For `eps=+1`, the only change in the rational ramification is that the totally
ramified place is `R=1`; the total different degree and genus are unchanged.

## The unique involution and genus-6 quotient

The coefficients of the cubic are fixed term by term by

```text
iota: R -> (5R-7)/(3R-5),     rho -> rho.
```

An invariant is

```text
z = (3R^2-7)/(3R-5),
```

and the quadratic relation over the invariant field is

```text
R^2-zR+(5z-7)/3 = 0.
```

Writing `A,B,C` in the basis `[1,R]` over `Q(z)` gives zero `R` coefficient
identically.  Hence the cubic descends exactly to `Q(z)`.  For `eps=-1`, its
normalized ramification is:

| base place | degree | ramification | different contribution |
|---|---:|---:|---:|
| `z=5` | 1 | `e=2` | 1 |
| `z=2` | 1 | `e=2` | 1 |
| `z=14/3` | 1 | `e=3` | 2 |
| `31z^2-220z+352=0` | 2 | `e=2` | 2 |
| `P10(z)=0` | 10 | `e=2` | 10 |

Here

```text
P10 = 675105*z^10 - 24788818*z^9 + 404115156*z^8
      - 3847306032*z^7 + 23658817344*z^6
      - 98072909568*z^5 + 277203253248*z^4
      - 526965620736*z^3 + 644275372032*z^2
      - 457272459264*z + 143107555328.
```

Thus the different has degree `16`, and

```text
2g-2 = -6+16 = 10,
g = 6.
```

The sign-changing isomorphism descends to

```text
z -> (47z-220)/(10z-47).
```

## Automorphisms

A genus-6 trigonal curve has a unique trigonal pencil: two independent
degree-three maps to `P1` would contradict the Castelnuovo--Severi bound
`g <= (3-1)(3-1)=4`.  Therefore every `Q`-automorphism of the quotient must
descend to a `Q`-automorphism of `P1_z` preserving the ramification data.

The unique `e=3` rational branch is `z=14/3`, while the rational `e=2`
branches are `z=2,5`.  Besides the identity, the only possible base map fixes
`14/3` and swaps `2,5`:

```text
z -> (106z-532)/(21z-106).
```

Direct substitution shows that it does not preserve the unique quadratic
branch place `31z^2-220z+352`.  Hence it cannot lift.  The cubic extension is
not Galois, so it has no nontrivial automorphism over `Q(z)`.  Consequently

```text
Aut_Q(genus-6 quotient) = 1.
```

The same argument on the genus-12 curve leaves only the identity and `iota`,
both of which lift.  Thus

```text
Aut_Q(genus-12 root cover) = <iota> ~= C2.
```

There is no further natural automorphism quotient to exploit.

## The actual Richelot-halving gates

After choosing the rational root `rho`, write

```text
A0 = Q2,
B0 = C3/(x-rho).
```

The discriminant squareclasses are generated by

```text
DA = disc(A0)
   = (B-3)^2 - 4(A+3)(C+1),

DB = disc(B0)
   = (A-3)^2 - 4(A-3)rho - 12rho^2 - 8(B+3),

DA*DB.
```

All three expressions are invariant under `iota`, so they descend to the
genus-6 quotient.  Exact quadratic function-field normalization gives:

```text
genus(D(sqrt(DA)))       = 17,
genus(D(sqrt(DB)))       = 18,
genus(D(sqrt(DA*DB)))    = 24.
```

The exact gates from the Richelot calculation are:

```text
U: DB is a square,                         genus 18,
V: DA is a square,                         genus 17,
W: DA and DB are both squares.
```

For `W`, adjoining both square roots gives a `V4` cover of the genus-6 base.
Its three index-two quotients are precisely the genus `17`, `18`, and `24`
curves above.  The standard `V4` genus relation therefore gives

```text
g(W) = 17+18+24-2*g(D) = 17+18+24-12 = 47.
```

For comparison, the first two covers before quotienting have genera `35` and
`36`.  The genus-24 product-square quotient is already a stop certificate;
the actual simultaneous `W` gate is substantially worse at genus `47`.

## Computational search context

The independent projective sieve through height `25000` found only
`R=1,2,3,7/3`, all boundary degeneracies; see

```text
notes/contact30_c3root_projective_search_2026_07_11.md
```

That is a bounded exclusion rather than a global point theorem.  Combined
with the geometry above, it is enough to recommend parking this route rather
than beginning another high-genus descent or Chabauty computation.
