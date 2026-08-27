# Packaged simple `[6,6]` contact-6 example

The specialization

```text
a = 133/39,  b = -7/13,
L = 29/16,  U = -9/4,  v = 5/2
```

on the contact-6 core gives the PARI-minimal, reduced integral model

```text
C: y^2 = 1872*x^5 - 3000*x^4 + 6969*x^3 - 1691*x^2 + 4875*x.
```

This is obtained from the earlier integral equation by dividing its right
side by `78^2`.  It factors over `Q` as

```text
1872*x*(x^2 - 23*x/13 + 125/39)*(x^2 + x/6 + 13/16).
```

## Exact torsion certificate

On the minimal model, use the Mumford classes

```text
D = (x - 1, 95),
E = (x^2 - 9*x/4 + 25/4, 25*x - 300),
T = (x, 0),
W = E + T.
```

Exact Jacobian arithmetic gives orders `6`, `3`, `2`, and `6`, respectively.
The 36 classes `i*D+j*W`, for `0 <= i,j < 6`, are distinct.  Moreover,

```text
#J(F_7)  = 36,
#J(F_11) = 144.
```

Both primes are of good reduction.  Reduction is injective on prime-to-`p`
torsion; using the two different residue characteristics shows that the order
of `J(Q)_tors` divides `gcd(36,144)=36`.  The explicit subgroup therefore
proves

```text
J(Q)_tors = Z/6Z x Z/6Z.
```

The marked classes come directly from the identities checked in the Sage
script:

```text
H6^2 - F = (39/2)^2*(x-1)^6,
H3^2 - F = (97344/841)*(x^2 - 9*x/4 + 25/4)^3.
```

The original Magma search independently printed exact invariants `[6,6]` and
orders `ord(D)=6`, `ord(E)=3` in
`data/contact6_m36_core_slice_h14_simple.txt`.  The dedicated verifier is
`code/contact6_m36_66_package.m`.  It was prepared but could not be rerun on
2026-07-10 because Magma was not installed on this system.  The Sage certificate above is an independent exact proof, not a
replacement assertion based on the old output.

## Geometric simplicity

Sage 10.7 gives

```text
geometric_endomorphism_algebra_is_field(B=100): True
geometric_endomorphism_ring_is_ZZ(B=100): True
```

There is also a standalone good-reduction certificate.  At `p=37`,

```text
P_37(X) = X^4 - 10*X^3 + 54*X^2 - 370*X + 1369,
```

and the polynomial whose roots are the 12th powers of the roots of `P_37` is

```text
V^4 + 5679639644*V^3 + 16271908860833814246*V^2
    + 37388795186918383904326279964*V
    + 43335257111193343900365036083324748961.
```

It is irreducible over `Q`.  This is the Frobenius-power test in Lombardo's
Algorithm 4.10 and proves geometric simplicity.  In addition, the minimal
model discriminant factors as

```text
2^8 * 3^5 * 5^12 * 13^4 * 17^3 * 19^6 * 29,
```

so the exponent-one odd factor used by the same implementation certifies that
the geometric endomorphism ring is `ZZ`.

Run and capture the complete reproducible certificate with

```bash
sage code/contact6_m36_66_package.sage \
  > data/contact6_m36_66_package_sage.txt
```
