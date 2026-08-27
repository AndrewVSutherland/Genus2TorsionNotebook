# Exact marked section of the contact-9 `[3,18]` q-cover

This independent certificate extracts the rational section of the unsigned
order-3 q-cover coming from the known class `3*D9`, verifies its contact
identity over `Q(t)`, and removes it scheme-theoretically from the degree-40
cover.

## Files and commands

```text
code/audit_contact9_318_marked_section.sage
data/audit_contact9_318_marked_section_Q.txt
data/audit_contact9_318_marked_section_strong_Q.txt
```

The recorded files were produced by:

```bash
timeout 60s sage code/audit_contact9_318_marked_section.sage \
  --section-only > data/audit_contact9_318_marked_section_Q.txt

timeout 180s sage code/audit_contact9_318_marked_section.sage \
  --strong-check > data/audit_contact9_318_marked_section_strong_Q.txt
```

Both runs completed with exit status zero.

## Symbolic extraction

Let `D9=[(1,h(1))-infinity]` on the rational-root contact-9 family.  Sage
computes `3*D9` over `Q(t)` by exact Cantor arithmetic.  Its monic Mumford
polynomial is

```text
q = x^2 + u0*x + v0,

u0 = (3*t^5 + 15/4*t^4 - 3*t^3 - 147/16*t^2
      - 27/4*t - 27/16)/(t+1)^4,

v0 = (-2*t^5 - 47/16*t^4 + 1/4*t^3 + 7/2*t^2
      + 11/4*t + 11/16)/(t+1)^4.
```

Substitution of `u0,v0` into the three triangular q-cover equations gives a
common monic gcd of degree one in `z`.  Its root is

```text
z0 = (4*t^10 + 14*t^9 + 81/4*t^8 + 27/2*t^7 + 9/8*t^6
      - 9/2*t^5 - 183/64*t^4 - 3/8*t^3 + 9/32*t^2
      + 1/8*t + 1/64)/(t+1)^8.
```

It is visibly a square:

```text
z0 = s^2,
s = 2*(t+1/2)^2*(t^3 + 3/4*t^2 - 1/4)/(t+1)^4.
```

The full monic cubic `G` is printed in the data file.  The script verifies
coefficient by coefficient that

```text
G^2 - z0*f - q^3 = 0.
```

All three eliminated q-cover residuals are also exactly zero.  With
`H=G/s`, reduction modulo `q` gives the negative of the Mumford `v`
polynomial of `3*D9`.  Thus this is precisely the unsigned pair
`{3*D9,-3*D9}`, with the displayed normalization choosing `-3*D9`.

The script additionally checks exactly that

```text
z0 != 0,  disc(q) != 0,  Res(q,f) != 0,
```

so the section is on the generic open q-cover.

## Specialization at `t=4 mod 7`

The section specializes to

```text
(u,v,z) = (1,4,2),
q = x^2+x+4.
```

This is exactly the q-pair previously identified in the generalized
Jacobian calculation as the marked `3D9` branch.  The other two smooth
q-pairs in that calculation remain in the degree-39 residual cover.

## Exact removal and degree 39

As in the prior degree-40 computation, two quotients by the joint boundary
ideal `(z,u^2-4v)` change the raw ideal first from dimension one to dimension
zero.  The resulting finite algebra has degree 40.  Let

```text
m = (u-u0, v-v0, z-z0).
```

Every generator of the degree-40 ideal reduces to zero modulo `m`, and `m`
has degree one.  The exact ideal quotient

```text
Ires = I40 : m
```

has vector-space degree 39.  The strong run further certifies

```text
m + Ires = (1),
m intersect Ires = I40.
```

Hence this is a genuine reduced, comaximal degree-1 plus degree-39
decomposition, not just subtraction of point counts.

The strong run deliberately uses the already certified double-boundary
quotient of degree 40 and skips the more expensive post-saturations by
`disc(q)` and `Res(q,f)`.  The marked section itself is checked to be open;
the earlier q-cover computation recorded that the postfilters do not alter
the target degree.

## Expected `1+12+27` structure

The degree pattern has a group-theoretic explanation.  The unsigned cover
parametrizes the 40 pairs of nonzero elements of `J[3]` modulo sign.  If
`T=3D9`, then:

```text
{+/-T}:                                      1 pair,
(T^perp - {0,+/-T})/{+/-1}:                 12 pairs,
(J[3] - T^perp)/{+/-1}:                     27 pairs.
```

Thus removal of the marked section must precede any search for a `12+27`
split.  This audit certifies the `1+39` decomposition exactly.  It does not
claim an exact global `12+27` ideal decomposition; that requires adjoining
an explicit Weil-pairing condition or completing the residual factorization.
