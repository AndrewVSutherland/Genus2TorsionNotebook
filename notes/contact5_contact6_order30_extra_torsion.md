# Extra Torsion Conditions on the Contact-5/Contact-6 Order-30 Family

This records algebraic conditions for asking whether the order-30 family has
extra independent `2`- or `3`-torsion.

## Extra Independent 2-Torsion

For every member of the family,

```text
f = h6^2 - (x-1)^6
  = (h6-(x-1)^3)*(h6+(x-1)^3)
  = Q2*C3,
```

where

```text
Q2 = (A+3)*x^2 + (B-3)*x + (C+1),
C3 = 2*x^3 + (A-3)*x^2 + (B+3)*x + (C-1).
```

This `2+3` factorization gives the built-in rational `2`-torsion class.  An
additional independent rational `2`-torsion class occurs exactly when the
factorization refines over `Q`, i.e.

```text
1. Q2 splits over Q, or
2. C3 has a rational root rho.
```

For the explicit parameter `R` and branch `eps=+/-1`, the quadratic splitting
condition is the hyperelliptic square condition below.

For `eps=-1`:

```text
W^2 = -1095*(R^2 - 4*R + 13/3)
      *(R^6 - 4044/365*R^5 + 18249/365*R^4
        - 42664/365*R^3 + 54039/365*R^2
        - 34764/365*R + 1751/73).
```

For `eps=+1`:

```text
W^2 = -2135*(R^2 - 32/7*R + 37/7)
      *(R^6 - 3504/305*R^5 + 16761/305*R^4
        - 42784/305*R^3 + 61611/305*R^2
        - 47664/305*R + 3119/61).
```

The cubic-root condition is simply

```text
2*rho^3 + (A-3)*rho^2 + (B+3)*rho + (C-1) = 0,
```

after substituting `A,B,C` from the `R,eps` family.  After clearing
denominators this is irreducible in the tested factorization; it has degree
`20` in `R` and degree `3` in `rho`.

A direct exact scan through rational parameter height `120` found no extra
`2`-torsion specialization on either branch.

## Extra Independent 3-Torsion

For a genus-2 curve `y^2=f(x)=sum c_i*x^i`, a rational `3`-torsion class with
quadratic Mumford support is equivalent to a cubic contact:

```text
q3   = x^2 + U*x + V,
ell3 = M*x^3 + N*x^2 + P*x + S,
ell3^2 - f = M^2*q3^3.
```

Equating coefficients gives:

```text
2*M*N - c5 - 3*M^2*U = 0,
N^2 + 2*M*P - c4 - M^2*(3*U^2+3*V) = 0,
2*M*S + 2*N*P - c3 - M^2*(U^3+6*U*V) = 0,
P^2 + 2*N*S - c2 - M^2*(3*U^2*V+3*V^2) = 0,
2*P*S - c1 - 3*M^2*U*V^2 = 0,
S^2 - c0 - M^2*V^3 = 0.
```

For this family, substitute

```text
f = (x^3 + A*x^2 + B*x + C)^2 - (x-1)^6.
```

The existing `3`-torsion is the non-reduced cubic contact

```text
q3 = (x-1)^2, so U=-2, V=1, ell3=+/-h6.
```

Thus an additional independent `3`-torsion class must satisfy the contact
equations away from this branch, and then be checked for independence in
`J[3]`.

## Finite-Prime Diagnostic

At every good prime `p != 3`, a second independent rational `3`-torsion class
forces

```text
9 | #J(F_p).
```

The finite scan found:

```text
p=11: eps=-1 good=3 allowed9=0; eps=+1 good=3 allowed9=0
p=13: eps=-1 good=4 allowed9=0; eps=+1 good=4 allowed9=0
```

So any rational specialization with extra independent `3`-torsion must reduce
to the bad/boundary locus at both `p=11` and `p=13`.  This does not prove
nonexistence, but it means a rational search should be boundary/CRT-guided,
not a blind height scan.

## Files

```text
code/contact5_contact6_order30_extra_conditions.py
code/contact5_contact6_order30_extra3_residue_scan.m
data/contact5_contact6_order30_extra_conditions.txt
data/contact5_contact6_order30_extra3_residue_scan.txt
```
