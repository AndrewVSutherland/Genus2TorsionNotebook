# Contact-7 with two rational roots plus a 3-class

This is the bounded-search and local-sieve audit of the proposed
`[2,42]` record lane.  It complements the cubic-contact geometry in
`code/contact7_two_root_plus3_geometry.m`.

## Signed two-root coordinates

The four apparent sign choices are redundant.  Put

```text
u = eps*s,       v = delta*t.
```

Then

```text
r = 1-u^2,       w = 1-v^2,
h(r) = u^7,      h(w) = v^7.
```

Thus the search space is an unordered rational `(u,v)`-surface, not a
fourfold sign cover.  The cancelled parametrization is essential.  Put

```text
A(z) = (2z^5+4z^4+6z^3+8z^2+10z+5)/(2(z+1)^2).
```

Then

```text
b = (A(v)-A(u))/(u^2-v^2),
a = A(u)-b*(1-u^2).
```

After cancellation, `u=+1` and `v=+1` are removable and can give
smooth curves with a root at zero.  The genuine parameter poles are
`u=-1`, `v=-1`, and `u=-v`; `u=0`, `v=0`, and `u=v` are singular.

For

```text
h = 1-(7/2)x+a*x^2+b*x^3,
f = (h^2+(x-1)^7)/x^2,
```

the curve `y^2=f` has two independent rational Weierstrass classes and
the marked contact-7 class.  Exact Magma checks give visible orders
`7,2,2,2`, hence the generic visible subgroup is
`Z/2 x Z/14`, of order 28.  An additional rational 3-class would give
`Z/2 x Z/42`, of order 84.

For independent implementation checks, the degree-five polynomial is

```text
f = x^5 + (b^2-7)x^4 + (2ab+21)x^3
      + (a^2-7b-35)x^2 + (2b-7a+35)x + (2a-35/4).
```

## Conservative finite-field masks

At a good prime `p`, the prime-to-`p` part of 84 must divide
`#J(F_p)`.  The C++ sieve computes `#C(F_p)` and `#C(F_{p^2})`, then
recovers `#J(F_p)`.  A residue is discarded only when:

1. both parameters lie in the displayed affine chart;
2. the two-root/contact conditions remain nondegenerate;
3. `f mod p` is squarefree of degree five; and
4. the necessary target divisibility fails.

All denominator, chart-boundary, and bad-reduction residues are retained.
Magma independently reproduced the mask counts at `p=5,7,11,13` and the
individual Jacobian orders.

The most informative local fact is at `p=5`:

```text
good affine residues = 4
target-84 passes     = 0
#J(F_5)              = 28 for all four
```

Consequently a global record hit on this surface would have to reduce
into a 5-adic parameter-boundary disk.  This is not by itself a global
obstruction: the sieve deliberately retains all such disks.

At `p=7` there are 16 good residues, two of which pass the target
condition.  The earlier apparent `p=7` obstruction came from using the
uncancelled raw formula and incorrectly omitting `u=+1` and `v=+1`.

## Complete height-100 search

Height means `u=n/d` in lowest terms with
`|n| <= 100` and `1 <= d <= 100`, and similarly for `v`.  Root exchange
was removed by enumerating unordered pairs.

```text
rational parameters       12,175
unordered raw pairs    74,109,225
exact chart boundaries     30,433
open pairs tested       74,078,792
```

Conservative masks at

```text
5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83
```

left 167 candidates.  Exact rational reconstruction followed by
independent Magma reduction at good primes from 89 through 251 killed
all 167.  The first-kill distribution was

```text
89:74, 97:60, 101:14, 103:8, 107:6,
109:1, 113:2, 127:2.
```

Thus there is no target-84 candidate in the complete signed unordered
height-100 box on the open two-root chart.

Two intermediate height-30 survivors were also sent through full
`TorsionSubgroup` and geometric-simplicity checks:

```text
(u,v)=(-20/9,21/26): torsion [2,14], order 28, D4 at p=17;
(u,v)=(-30/23,-4/5): torsion [2,14], order 28, D4 at p=31.
```

Both are simple false positives, not record examples.

## Interpretation

The bounded experiment is decisively negative, but it does not prove
that the whole surface has no rational 3-class.  The arithmetic forces
any global hit into a 5-adic boundary disk.  The removable `u=+1` branch
has now been searched separately through height 10,000; see the boundary
addendum.  The remaining exact affine exclusions are genuinely
singular or incompatible, not omitted smooth families.

## Files

```text
code/contact7_two_root_plus3_sieve.cpp
code/contact7_two_root_plus3_verify.m
data/contact7_two_root_plus3_h100_corrected.txt
data/contact7_two_root_plus3_h100_corrected_survivors.m
data/contact7_two_root_plus3_h100_exact_summary.txt
notes/contact7_two_root_plus3_boundary_addendum_2026_07_11.md
```
