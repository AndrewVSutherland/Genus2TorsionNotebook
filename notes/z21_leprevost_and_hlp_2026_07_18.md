# Explicit infinite cyclic `[21]` families

## Outcome

There are two viable routes, and the stronger one was already in the primary
literature but had not been transcribed into this workspace.

1. Leprévost (1991) gives an explicit one-parameter genus-2 family with a
   marked rational class of order `21`.  The local verifier now checks its
   order exactly over `Q(t)` and finds a geometrically simple specialization
   whose complete rational torsion is exactly `[21]`.
2. HLP Lemma 11 and Proposition 4 give a very compact one-parameter subfamily
   of their two-parameter elliptic-gluing construction.  It also has a
   specialization with exact torsion `[21]`, but every member is geometrically
   split by construction.

The direct contact-7-plus-3 searches already in the repository were attempted
ansätze for order `42`, not an established order-`21` family, and found no
hits in their bounded searches.

## Leprévost's one-parameter family

Put

```text
p2 = t^14+4*t^13+19*t^12+32*t^11+113*t^10+188*t^9
     +379*t^8+448*t^7+379*t^6+188*t^5+113*t^4+32*t^3
     +19*t^2+4*t+1,

p1 = t^10+4*t^9+17*t^8+24*t^7+54*t^6+56*t^5+54*t^4
     +24*t^3+17*t^2+4*t+1,

p0 = t^6+4*t^5+15*t^4+16*t^3+15*t^2+4*t+1,

A21(x) = p2*x^2 - 2*(t^2+1)^2*p1*x + (t^2+1)^4*p0,

k21 = 64*t^4*(t+1)^2*(t^2+1)^3
      *(t^4+2*t^3+6*t^2+2*t+1)^3.
```

Then

```text
C_t: y^2 = f21(x) = A21(x)^2-k21*x^3*(x-1)^2
```

is generically a smooth quintic genus-2 model with one rational point at
infinity.  The point `(0,A21(0))` defines

```text
D21 = [(0,A21(0))-infinity].
```

Magma arithmetic over `Q(t)` verifies

```text
21*D21 = 0,
3*D21 != 0,
7*D21 != 0.
```

Thus `D21` has exact order `21` generically, without relying only on a
specialization.

In fact the complete generic rational torsion is exactly `[21]`.  The family
has good reduction at `t=1`, and torsion over `Q(t)` injects into that
characteristic-zero fiber.  Since the fiber has full rational torsion `[21]`
and the generic marked class already supplies `[21]`, there is no room for
additional generic rational torsion.

### Explicit exceptional set

The primary paper leaves the bad specializations implicit.  Direct exact
factorization gives

```text
Disc_x(f21) = -2^44
 * t^33*(t+1)^30*(t^2+1)^26*(t^2+t+1)
 * (t^4+2*t^3+6*t^2+2*t+1)^8
 * (t^6+4*t^5+15*t^4+16*t^3+15*t^2+4*t+1)^3
 * R24(t),
```

where

```text
R24 = t^24+5*t^23+37*t^22+81*t^21+104*t^20-1211*t^19
      -6479*t^18-23783*t^17-59305*t^16-117982*t^15
      -191446*t^14-249302*t^13-286592*t^12-249302*t^11
      -191446*t^10-117982*t^9-59305*t^8-23783*t^7
      -6479*t^6-1211*t^5+104*t^4+81*t^3+37*t^2+5*t+1.
```

All displayed non-linear factors are irreducible over `Q`.  Hence for rational
affine parameters the only singular values are `t=0,-1`; over `Qbar`, exclude
the roots of every displayed factor.  The compactified parameter also has a
degenerate fiber at `t=infinity`.

### Small exact specialization

At `t=1`, divide `y` and `A21` by `128`.  Since

```text
A21/128 = 15*x^2-16*x+7,
```

one obtains

```text
C_1: y^2 = -216*x^5+657*x^4-696*x^3+466*x^2-224*x+49,
D_1 = [(0,7)-infinity].
```

Exact Magma output is

```text
Order(D_1) = 21,
J(C_1)(Q)_tors = [21].
```

At the good prime `5`, `#J(F_5)=21`; the checks also give orders `126`,
`252`, `336`, and `483` at `11,13,17,19`.

The independent Sage/Lombardo check gives

```text
geometric_endomorphism_algebra_is_field(B=100) = True
geometric_endomorphism_ring_is_ZZ(B=100) = True
```

Thus this exact `[21]` specialization is geometrically simple with geometric
endomorphism ring `ZZ`.  In particular, the generic Leprévost Jacobian is
geometrically simple; this is a local computational certificate, not a claim
made in the 1991 note.

The 1991 paper proves that the family contains infinitely many pairwise
non-`Qbar`-isomorphic curves.

## Bridge to the local contact-7 plus cubic-contact model

The Leprévost family can be put exactly on the direct cover that motivated the
local search.  From

```text
div(y-A21) = 3*D21+2*P7
```

one gets `P7=9*D21`, hence an order-7 class, while `7*D21` has order 3.
The script `code/contact7_plus3_leprevost_bridge.m` constructs a rational
change of variable over `Q(t)` and verifies

```text
h7(z)^2-z^2*g(z) = -(z-1)^7,
H3(z)^2-q3(z)^3 = L3^2*g(z).
```

The first identity marks the order-7 class at `z=1`; the second gives the
independent cubic-contact order-3 class.  Thus the direct contact-7-plus-3
cover itself contains a one-parameter `[21]` curve, rather than merely isolated
search hits.

At `t=1` the normalized data are

```text
g  = z^5-(227/36)*z^4+(1271/81)*z^3-(27733/1458)*z^2
       +(298/27)*z-257/108,
h7 = 1-(7/2)*z+(86/27)*z^2-(5/6)*z^3,
q3 = z^2-2*z+73/81,
H3 = z^3-(53/18)*z^2+(8/3)*z-997/1458,
L3 = 1/3.
```

The script checks both identities over `Q(t)`, proves the marked classes have
orders 7 and 3 generically, and verifies full torsion `[21]` after specializing
at `t=1`.  It used about `83 MB` and `0.19` seconds under a hard `280 MB` cap.

## Compact HLP gluing family

For comparison, start with the Tate normal form on `X_1(7)`

```text
E7_t: y^2+(1-c)xy-b*y=x^3-b*x^2,
b=t^3-t^2, c=t^2-t,
P7=(0,0).
```

Write its short model as

```text
y^2=x^3+A*x+B.
```

With the standard Weierstrass quantities

```text
a1=1-c, b2=a1^2-4*b, b4=-a1*b, b6=b^2,
A=-(b2^2-24*b4)/48,
B=(b2^3-36*b2*b4+216*b6)/864,
```

HLP Lemma 11 supplies

```text
q=B^2/A^3,
E3_t: y^2=x^3+(x-q)^2,
P3=(0,q).
```

The Galois-equivariant correspondence between the nonzero 2-torsion points is

```text
beta = -alpha^2/A+(B/A^2)*alpha-1,
alpha^3+A*alpha+B=0.
```

Indeed `beta^3+(beta-q)^2=0`.  Applying HLP Proposition 4 and removing only a
rational square factor gives

```text
C_t^HLP: y^2 = A^7-3*A^5*B*x^2
                  +(A^6+3*A^3*B^2)*x^4-A*B^3*x^6.
```

Its discriminant is

```text
64*A^52*B^3*(4*A^3+27*B^2)^2.
```

The quotient map has 2-primary kernel, so the rational point `(P7,P3)` of
order `21` retains exact order `21` in the Jacobian.  At `t=2`, Magma gives
the full rational torsion `[21]`.

This HLP family is an even sextic and has elliptic quotients; its Jacobian is
`(2,2)`-isogenous to `E7_t x E3_t`.  It is therefore an explicit positive
control, not a geometrically simple family.

## Prior local status

Before this transcription, the local inventory had no `[21]` row.  The tracked
version of `notes/internet_infinite_families_sweep.md` only said that
Leprévost's order-`21` family existed and that its equations had not been
extracted; that note has now been updated.

The contact-7 rational-root plus 3-torsion attempt is recorded in
`notes/contact7_family.md`, lines 607-762, with outputs

```text
data/contact7_root_plus3_h12.txt
data/contact7_root_plus3_h20.txt
data/contact7_root_plus3_boundary3_h20.txt.
```

Those computations target the even group `[42]` and report zero hits, so they
are attempted ansätze rather than order-`21` examples or families.

## Artifacts and resources

```text
code/z21_leprevost_family_verify.m
code/leprevost_z21_geom_check.sage
code/contact7_plus3_leprevost_bridge.m
code/hlp_z21_family_verify.m
data/z21_leprevost_family_verify.txt
data/z21_leprevost_family_verify_time.txt
```

The generic and specialized Leprévost verification used `84,408 KB` peak RSS
and `0.26` seconds under a hard `200 MB` address-space limit.  The exact
discriminant factorization used `58,916 KB`.  The Sage check was run under a
hard `2 GB` address-space limit on a host with substantially more available
memory.

Reproduction commands:

```bash
magma -b MemMB:=180 code/z21_leprevost_family_verify.m
magma -b code/contact7_plus3_leprevost_bridge.m
magma -b MemMB:=180 code/hlp_z21_family_verify.m
prlimit --as=2000000000 sage code/leprevost_z21_geom_check.sage
```

Primary references:

- F. Leprévost, *Familles de courbes de genre 2 munies d'une classe de
  diviseurs rationnels d'ordre 15, 17, 19 ou 21*, C. R. Acad. Sci. Paris
  Sér. I Math. 313 (1991), 771-774.  Official scan:
  <https://gallica.bnf.fr/ark:/12148/bpt6k57325582/f775.image>.
- E. Howe, F. Leprévost, B. Poonen, *Large torsion subgroups of split
  Jacobians of curves of genus two or three*, Forum Math. 12 (2000), 315-364:
  <https://math.mit.edu/~poonen/papers/large.pdf>.
