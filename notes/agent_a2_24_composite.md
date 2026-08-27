# [2,24] via coprime composition 24 = 8 x 3 on the A(8) chart

Date: 2026-07-03.  After the A(12)-halving, Pell-scan, and backward-CF
routes all confirmed [2,24] is high-height/construction-only, this route
sidesteps BOTH walls.

## The idea

The A(8) chart carries a visible order-8 class D8 for FREE.  Any curve
on it that ALSO has a rational 3-torsion class has, automatically,
Z/8 + Z/3 = Z/24 in its torsion (coprime orders always compose) -- no
order-8->16 halving (the [4,16] wall) and no degree-24 Pell unit (the
CF-construction wall).  For the [2,24] target, additionally impose
2-rank 2 (the extra [2] factor).

## Search (agent_a2_24_composite8x3.m)

Sweep the A(8) chart (r,p,t); skip even sextics (bielliptic => split);
prefilter 3 | #J(F_q) at the strong prime set {7,11,13,17,19,23,29,31,
37,41} (necessary for rational 3-torsion, no false negatives); exact
TorsionSubgroup on survivors; 2-rank filter for the [2,24] variant.

## VALIDATED: geometrically simple cyclic Z/24 (H=12)

Two geometrically simple Z/24 curves (2-rank 1), independently verified
(torsion [24]; irreducible Frobenius with D4 Galois group at many primes):

```text
(r,p,t) = (5, -5/2, -9/2):
  f = 701706240 x^6 - 2463436800 x^5 - 1709890560 x^4
      - 58438195200 x^3 + 367508378880 x^2 - 494073043200 x - 6779350080
  torsion [24]; SIMPLE (chi_17 = x^4-5x^3+16x^2-85x+289, irred; D4 at 8 primes)

(r,p,t) = (1/3, -1/9, -1):
  f = 1549681956 x^6 - 2066242608 x^5 + 1607077584 x^4 - 1836660096 x^3
      + 2040733440 x^2 - 1088391168 x + 241864704
  torsion [24]; SIMPLE (chi_13 = x^4+2x^3-6x^2+26x+169, irred; D4 at 8 primes)
```

(Z/24 itself is already known -- in LMFDB -- so these validate the route
rather than being new; the target is [2,24].)

## [2,24] hunt: in progress

Strong 3-prefilter + 2-rank>=2 filter, height 20 (agent_a2_24_composite_h20_*).
The bottleneck is that 2-rank-2 AND rational-3-torsion are each rare on
A(8), so their intersection (the [2,24] locus) is thin -- a codim-2
1-dim family in (r,p,t).  If the blind-scan intersection stays thin, the
denser route is to parametrize the 2-rank-2 A(8) sublocus (W = Q^2+q
splits into two rational quadratics, quadratic in p for fixed (r,t,beta)
-- from the [4,16] work) and apply only the 3-prefilter on it.

## [2,24] hunt result: 3-torsion is the wall (2026-07-03)

Dense W-split + 3-torsion search (agent_a2_24_wsplit_3tors.m):

```text
part 0 (of 3): 87.7M (r,t,beta) tested, 991,275 curves 2-rank-2 (~90%!),
1709 passed the 14-prime 3-prefilter, ALL torsion-computed:
  histogram [2,8]:1708, 0 order-24, 0 [2,24].
```

So the 2-rank-2 order-8 structure ([2,8], the 2-part of [2,24]) is
ABUNDANT and dense, but rational 3-torsion on the 2-rank-2 A(8) locus is
absent through ~1M samples.

**Structural reading.**  The composite route decomposed [2,24] into:
 - order-8 (2-adic exponent): free on A(8);
 - extra Z/2 (2-rank 2): free via W-split (dense, [2,8] everywhere);
 - Z/3 (3-torsion): the wall.
The two simple Z/24 curves found earlier had order-8 + rational 3-torsion
but 2-rank 1.  The W-split ([2,8]) curves have order-8 + 2-rank-2 but no
3-torsion.  The two sub-loci (3-torsion; W-split) barely intersect: the
[2,24] locus is codim-2 (a 1-dim family in (r,p,t)) that appears
rational-point-poor at accessible height -- rational 3-torsion is
anti-correlated with 2-rank 2 here.

No local obstruction (F_13..F_29 all [2,24]-compatible), so simple
[2,24] exists; it is just a sparse global point.

## Milestone achieved and options

ACHIEVED: two geometrically simple cyclic Z/24 genus-2 Jacobians
(validated, D4 certificates); dense simple [2,8]; no local obstruction
to [2,24].

Next options for [2,24] (all real work):
 1. Parametrize the 3-torsion sublocus of the W-split family (cubic
    contact h3^2 - f = kappa*q3^3 on f = q*W1*W2) and study the geometry
    of the resulting 1-dim [2,24] family: if rational, [2,24] falls; if
    high genus, it is genuinely hard.
 2. Large distributed brute run on the W-split+3-torsion search.
 3. Literature [2,2n] family (Leprevost/Kulesz) -- the source of the
    known [2,26],[2,28].
