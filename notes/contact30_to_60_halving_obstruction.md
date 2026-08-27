# The contact-(5,6) order-30 family cannot yield order 60 by halving

## Construction tested

For the simultaneous contact-5/contact-6 family in
`code/contact5_contact6_order30_family.m`, write

```text
f = Q2*C3,
Q2 = h6-(x-1)^3,
C3 = h6+(x-1)^3.
```

After making `Q2` monic, the distinguished rational 2-class is

```text
T2 = [Q2,0] = 3*D6 = 15*(D5+D6).
```

Here `D5`, `D6`, and `D5+D6` have orders 5, 6, and 30.  If `2*H=T2`,
then `H` has order 4 and

```text
Q60 = H + 16*(D5+D6)
```

has order 60.  Thus this was a direct route from the exact order-30 family
to the target.

## Exact halving cover

Put

```text
u = Monic(Q2) = x^2+u1*x+u0,
g = f/u = g3*x^3+g2*x^2+g1*x+g0.
```

A half exists precisely when there are `M != 0`, `N` and a monic quadratic
`G` such that

```text
u*(M*x+N)^2-g = M^2*G^2.
```

With `z=N/M` and `k=1/M^2`, coefficient comparison reduces this to two
equations:

```text
r = z+(u1-g3*k)/2,
s = ((u1+g3*k)*z+u0-g2*k-(u1-g3*k)^2/4)/2,
F1 = 2*r*s-(u1*z^2+2*u0*z-g1*k) = 0,
F0 = s^2-(u0*z^2-g0*k) = 0,
```

where `k` must be a nonzero rational square.  The resultant in `z` has
degree 8 in `k`.  On either contact-family branch, its cleared global form
has bidegree `(352,8)` in `(R,k)`; the `branch=-1` form has 3,177 terms.

The fiberwise implementation is
`code/contact30_to_60_halving_search.m`.  It reconstructs every candidate
half, verifies `2*H=T2` in the Jacobian, computes the exact torsion subgroup,
and applies the full D4/root-power geometric-simplicity test to any hit.
The compact cover and Magma's independent `IsDivisibleBy` agreed on every
smooth height-2 test fiber and at `R=5/2` on both branches.

## Projective finite-field obstruction

The local masks must include the projective parameter `R=infinity`; a
rational parameter whose denominator is divisible by `p` reduces there.
After this correction, the `p=11` table on each branch is

```text
P^1(F_11): 12 residues
good:       4
halvable:   0
boundary:   8
```

The four good residues are `R=0,8,9,infinity`.  At all four, the order-6
class `D6` is not divisible by 2 in `J(F_11)`.  The boundary residues are

```text
R=1,6:      f reduces to zero,
R=2,3:      f = x^2*(x+1)*(x+4)*(x+9),
R=4,5,7,10: coefficient-pole charts.
```

The same class-specific gate is much sharper than checking only whether the
finite Jacobian exponent is divisible by 60.  For example, at `p=13` there
are good fibers with exponent divisible by 60, but none in which `T2` is a
double.

## Complete Q_11 obstruction, including the boundary

Center at the contact-6 point by putting `T=x-1` and

```text
h = T^3+Ac*T^2+Bc*T+Cc,
f = (h-T^3)*(h+T^3).
```

The 2-descent image of `D6=[T,Cc]` is the class of `T` in the two etale
factors.  If `D6` is a double, `T` must be a square in each factor.  Taking
norms gives the necessary conditions

```text
Nminus = Cc/Ac,
Nplus  = -Cc/2.
```

They resolve six of the eight boundary disks:

- In the disks `R=2,3`, `Nminus` reduces to `6`, a nonsquare in `F_11`.
  The reduced four-equation halving cover independently has no point.
- In every punctured disk around `R=1` and `R=7/3` (the latter reduces to
  `6`), `Nplus` has even valuation and nonsquare leading unit.  The centers
  themselves are degenerate.

For the four pole disks, the first norm leaves deeper directions.  The full
quadratic Kummer component closes them.  After normalizing `h-T^3`, its
reduction on every pole annulus is

```text
lambda*(T+1)^2.
```

Consequently the quadratic algebra is split or ramified with residue field
`F_11`, and the image of `T` reduces to `-1`.  Since `-1` is a nonsquare in
`F_11`, `T` cannot be a square in that factor.  The first and second blow-up
charts, including assertions for this double-root Kummer reduction, are in
`code/contact30_to_60_p11_boundary.m`.

Therefore neither branch has a nondegenerate parameter
`R in P^1(Q_11)` for which `D6`, `T2`, or `D5+D6` is divisible by 2.

## Conclusion

The contact-(5,6) order-30 family cannot produce a rational point of order
60 by halving its rational 2-class.  This is a complete local impossibility
result, not a bounded-search failure.  As an independent check, the
projective finite sieve found no survivor among 48,927 rational parameters
of height at most 200 on either branch.

The order-60 search should instead move to a construction carrying an order-4
class from the outset, such as the full two-parameter `M(12)` surface plus
5-torsion, or to the order-20 plus 3-torsion locus.
