# Contact-9 rational-root scout for `[3,18]`

## Scope

The rational-root contact-9 family is

```text
r = 1 - s^2,
h(r) = eps*s^9,  eps = +/-1,
a = (eps*s^9 - h0(r))/r^4,
h0(x) = 1 - (9/2)x + (63/8)x^2 - (105/16)x^3,
f = (h^2 + (x - 1)^9)/x^4.
```

It has the marked order-9 contact class and a rational Weierstrass point,
so its rational torsion generically contains `[18]`.  This scout asks
whether it is locally possible to add an independent rational 3-torsion
class, giving `[3,18]`.

## Correct finite subgroup test

Divisibility of `#J(F_p)` by the target order is not a subgroup test.  At
a good prime `p`, the prime-to-`p` part of rational torsion injects into
`J(F_p)`.  For the target `[3,18]`, the required finite subgroup is

```text
p = 2:  [3,9]
p = 3:  [2]
p != 2,3: [3,18].
```

For each prime `ell` and exponent `k`, the implementation compares

```text
#{target invariant factors divisible by ell^k}
```

with the corresponding count for the full invariant factors of `J(F_p)`.
These inequalities for all `ell,k` are equivalent to finite abelian
subgroup containment.

Code:

```text
code/contact9_318_finite_scout.m
code/contact9_318_finite_scout.sage
```

The Magma version uses `AbelianGroup(J)`.  The independent Sage version
enumerates the reduced Mumford representatives, checks their number against
the order obtained from `#C(F_p)` and `#C(F_{p^2})`, and recovers the
Sylow invariant factors from the kernels of multiplication by `ell^k`.

Run:

```text
sage code/contact9_318_finite_scout.sage --prime-bound 7 \
    > data/contact9_318_finite_p7.txt
```

## Result

```text
p 3 required [2]    good 0
p 5 required [3,18] good 0
p 7 required [3,18] good 2 pass_embedding 0
```

At `p=7`, the only good affine root-chart residues are

```text
(s,eps) = (2,+1), (5,-1).
```

Both give `a=1`, `r=4`, and the same curve

```text
y^2 = x^5 + 6*x^4 + x^3 + 5*x + 2
```

over `F_7`.  Its finite Jacobian has full invariants

```text
[2,36].
```

Thus its 3-primary subgroup is cyclic of order 9, whereas `[3,18]`
requires `[3,9]` in the 3-primary part.  There is no subgroup embedding.
In this particular fiber even the weaker order test fails, since
`#J(F_7)=72` is not divisible by 54.

## Boundary classification at 7

The other 12 affine `(s,eps)` residues are boundary:

```text
r = 0:
  s = 1 or 6, either eps                         (4 residues)

singular:
  (s,eps,a,r) =
  (0,+/-1,0,1),
  (2,-1,4,4),
  (3,-1,6,6), (3,+1,4,6),
  (4,-1,4,6), (4,+1,6,6),
  (5,+1,4,4)                                    (8 residues)
```

The projective parameter residue `s=infinity` is not part of this affine
chart.  Therefore the finite result says:

* no rational `[3,18]` specialization can have 7-integral `s` and good
  reduction in the open root chart;
* any possible specialization must escape through one of the listed
  `r=0` or singular residue classes, or through `s=infinity`.

Resolving those cases would require a 7-adic blowup/stable-model analysis.
A further affine height search would only revisit residue classes already
excluded by the good-reduction injection.

## Decision

The proposed independent cubic-contact 3-torsion cover was not derived or
searched.  The prerequisite of good-open residues at every tested small
prime already fails at `p=7`.  This route should remain stopped unless a
separate project explicitly analyzes the 7-adic boundary charts.
