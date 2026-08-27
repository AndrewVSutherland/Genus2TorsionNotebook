# Boundary correction: the removable root-at-zero branch

## Correction

The raw formulas with denominator `r^2*w^2*(w-r)` made `u=+1` look
like a chart boundary.  It is actually removable.  The cancelled
parametrization

```text
A(z) = (2z^5+4z^4+6z^3+8z^2+10z+5)/(2(z+1)^2),
b = (A(v)-A(u))/(u^2-v^2),
a = A(u)-b*(1-u^2)
```

extends across `u=+1` and `v=+1`.  The finite-field masks and the full
height-100 surface search were rerun with this corrected model.

## Exact classification of the excluded affine branches

- `u=+1` is removable and generically smooth.  It gives the valid
  root `r=0` and fixes `a=35/8`.
- `u=-1` is incompatible: it asks for `h(0)=(-1)^7=-1`, while the
  contact-7 normalization has `h(0)=1`.
- `u=0` puts the forced root at the marked contact point.  Both
  `h(1)` and the curve discriminant vanish, so this is singular and
  does not retain the contact-7 class.
- `u=v` makes the two forced roots coincide and the curve discriminant
  vanish.
- `u=-v`, away from zero, asks for opposite values of `h` at the same
  root and is a genuine parameter pole.  At zero it is already in the
  singular contact-collision branch.

The same statements hold after swapping `u` and `v`.  The remaining
finite discriminant factors `H`, `R(u,v)`, `R(v,u)`, and `K` are all
singular fibers, not omitted smooth components.

The exact projective parameter at infinity is distinct from a
denominator residue disk and was not part of the finite `n/d`
enumeration.  It is nevertheless not a missed smooth fiber.  On the
valid `u=1` branch, as `t` tends to infinity,

```text
h/t   -> -x^3,
f/t^2 -> x^4.
```

Thus after the natural scaling `y -> y/t`, the limiting curve is
`y^2=x^4`, which is singular; the second forced Weierstrass point has
merged with infinity.  More generally, with one full-surface parameter
`u` tending to infinity and the other `v` finite,

```text
A(u)=u^3+2u+O(1/u),
b=-u+O(1/u),
a=(1-v^2)u+O(1),
f/u^2 -> x^2(x-(1-v^2))^2,
```

again a square and hence singular.  Directions with both parameters at
infinity have the same square-limit degeneration; their exceptional
diagonal and antidiagonal directions are already the coincident-root
and incompatible-pole loci.  By contrast, *residue* denominator disks
modulo each sieve prime were retained conservatively as unknown.

These assertions are checked symbolically by
`code/contact7_two_root_plus3_u1_boundary.m` and by the full surface
discriminant factorization in `code/contact7_two_root_boundary.m`.

## The valid one-parameter family

Set `u=1` and write the other signed parameter as `t`.  Then

```text
r = 0,                 w = 1-t^2,
a = 35/8,
b = -(8t^4+24t^3+48t^2+45t+15)/(8(t+1)^3),
h(1) = -t^2(8t^2+9t+3)/(8(t+1)^3).
```

The finite rational parameter exclusions visible in these formulas are
`t=-1` (pole), `t=0` (contact collision), and `t=1` (coincident root).
The projective point `t=infinity` has the singular limit just described.
The remaining factors of the exact discriminant describe singular
algebraic fibers.  For example, `t=2` is smooth and Magma verifies

```text
visible divisor orders = 7,2,2,2
J(Q)_tors              = [2,14]
```

with an ordinary `D4` geometric-simplicity certificate at `p=7`.

## Local masks and height-10,000 search

The corrected finite masks have a sharper obstruction at `p=5`:

```text
full surface: 4 good ordered residues, all fail target 84;
u=1 branch:   t=2,3 are good, and both have #J(F_5)=28.
```

Therefore any global record hit must lie in one of the `5`-adic disks
`t=0,1,-1,infinity`.  Those disks were retained conservatively.

The complete projective-height convention was

```text
t=n/d, gcd(n,d)=1, |n|<=10000, 1<=d<=10000.
```

The search checked

```text
121,589,943 rational parameters
121,589,940 valid characteristic-zero parameters
```

using conservative masks through `p=83`.  They left 86 candidates.
Exact reconstruction and independent Magma reductions through `p=251`
killed every candidate:

```text
89:42, 97:31, 101:7, 103:4, 107:1, 127:1.
```

Thus the formerly omitted smooth branch contains no target-84 example
through height 10,000.  This is a bounded negative result, not a proof
that the branch has no rational 3-class.

## Files

```text
code/contact7_two_root_plus3_sieve.cpp
code/contact7_two_root_plus3_verify.m
code/contact7_two_root_plus3_u1_boundary.m
data/contact7_two_root_plus3_u1_h10000.txt
data/contact7_two_root_plus3_u1_h10000_survivors.m
data/contact7_two_root_plus3_u1_h10000_exact_summary.txt
```
