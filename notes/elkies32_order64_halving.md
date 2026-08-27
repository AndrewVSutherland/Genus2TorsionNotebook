# Trying to Halve the Elkies `[32]` Point

Goal: in the reconstructed Elkies `[32]` family, look for a specialization where

```text
D = (0,1) - infinity
```

is divisible by `2` in the Jacobian.  Since the reconstructed family has
`16D=(r,0)-infinity`, such a specialization would give a rational point of
order `64`, provided the curve is nonsingular and the original point has exact
order `32`.

## Algebraic Halving Conditions

For an odd quintic model

```text
y^2 = F(x) = F0 + F1*x + ... + F5*x^5
```

with `F(0)=1`, the divisor `(0,1)-infinity` is divisible by `2` if and only if
there are rational `m,n` and a sign `k=+/-1` such that

```text
F(x) - (m*x^2+n*x+k)^2 = F5*x*(x^2 + alpha*x + beta)^2,
```

where

```text
alpha = (F4 - m^2)/(2*F5),
beta  = (F3 - 2*m*n - F5*alpha^2)/(2*F5).
```

Equivalently, the two equations are

```text
K2 = F2 - n^2 - 2*m*k - 2*F5*alpha*beta = 0,
K1 = F1 - 2*n*k - F5*beta^2 = 0.
```

These are the family-level order-64 conditions.  No coefficient search is used
before this point.

## Substitution into the Reconstructed `[32]` Family

In the reconstructed contact model

```text
y^2 = (a*x^3 + b*x^2 + c*x + 1)^2 - a^2*x^5*(x+1),
```

write `p=a-b+c-1` and `z=a/p`.  The first `[32]` rank condition gives

```text
a = z*p,
c = (-z^4*p + 4*z^3*p - 8*z^2*p + 8*z^2 + 8*z*p - 16)/(4*z^2),
b = (-z^4*p + 8*z^3*p - 12*z^2*p + 4*z^2 + 8*z*p - 16)/(4*z^2).
```

The full `[32]` component is the genus-0 bidegree `(4,4)` equation `C32(z,r)=0`
from `data/elkies32_reconstruct_conditions.txt`, together with the rational
common-root formula

```text
p = 4*Pnum(z,r)/Pden(z,r).
```

Thus the order-64 cover is explicitly

```text
C32(z,r) = 0,
p = 4*Pnum(z,r)/Pden(z,r),
K1(z,p,m,n) = K2(z,p,m,n) = 0
```

for one of the two signs `k=+/-1`, away from boundary denominators.

After clearing denominators, for each sign the two halving numerators have
sizes

```text
K2: degrees [24,6,6,2] in (z,p,m,n), 241 terms,
K1: degrees [32,8,8,2] in (z,p,m,n), 511 terms.
```

## Finite-Field Sanity Check

The cover has points modulo many small primes, so there is no immediate simple
local obstruction.  For example the affine checks gave:

```text
q: #C32_affine  #usable_p_formula  #usable_with_halving
7:   4   1   1
11: 11   4   3
13:  6   4   2
17: 13  10   5
23: 21  19  14
31: 27  20   8
43: 40  37  19
73: 66  64  29
```

This means a quick congruence obstruction is not available.

## Search Performed

I used Magma's degree-8 genus-0 parametrization of the reconstructed `[32]`
component.  For rational parameter `t=U/V`, the search first applied residue
filters from the exact halving equations modulo

```text
11, 13, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59.
```

Then it ran the exact rational halving equations only on survivors.

Result for height `max(|U|,|V|) <= 5000`:

```text
total reduced parameters scanned: 30,401,831
passed residue sieve: 82
nonsingular exact checks: 82
hits: 0
```

Every exact survivor was rejected by the two rational halving equations.

The saved search log is

```text
data/elkies32_order64_param_search_B5000.txt
```

## Correction: Scale of `p`

A later audit found that the copied factorized formula for `p` missed the rational unit from Sage factorization.  The correct formula is `p = 4*Pnum(z,r)/Pden(z,r)`.  The scripts and data have been regenerated with this correction.  In the corrected finite-field table, `p=11` gives no good-residue halving points, so any order-64 candidate would have to lie over a bad `11`-adic boundary residue.

## Caveat

This is a targeted parameter search, not a proof of nonexistence.  The sieve is
performed in Magma's chosen degree-8 parameter; at some primes that
parametrization has bad reduction or collapses.  The script avoids using the
obviously collapsed prime `17` in the parameter sieve, but the search should
still be viewed as evidence, not a theorem.

## Files

- `code/elkies32_order64_conditions.py` derives the halving equations and finite-field sanity table.
- `data/elkies32_order64_conditions.txt` records the derived equations and counts.
- `code/elkies32_order64_param_search.m` performs the parameter residue sieve and exact survivor checks.
- `data/elkies32_order64_param_search_B2000.txt` is the corrected smaller run.
- `data/elkies32_order64_param_search_B5000.txt` is the larger corrected run.

## Conclusion

We now have the correct algebraic order-64 cover for the reconstructed Elkies
`[32]` family.  The first targeted search did not find a rational specialization
where `D` halves.  Since the finite-field cover has plenty of points, the next
serious step is geometric: eliminate the halving variables over the genus-0
base and compute the resulting cover's genus/components, rather than simply
pushing the height bound.
