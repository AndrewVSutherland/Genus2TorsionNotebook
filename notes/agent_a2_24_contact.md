# Cubic-contact 3-torsion machinery (validated) — for [2,24] derivation

## Contact form (VALIDATED on both known simple Z/24 curves)

A rational order-3 class D=(u,v) (Mumford; u monic deg 2, v deg<=1) on
J(y^2=f), deg f = 6, gives the cubic-contact identity

        f  =  h3^2  +  kappa * q3^3

with q3 := u (the "contact conic"), h3 the "contact cubic" (deg 3),
kappa a constant, and the leading-coefficient relation

        lead(f) = lead(h3)^2 + kappa.

Geometric origin: div(y - h3) = (zeros of f - h3^2) - 3*inf+ - 3*inf-.
For order-3 we need f - h3^2 = kappa*q3^3, i.e. 3*(P1+P2) - 3*(inf+ +inf-)
where P1,P2 lie over the roots of q3 with y = h3.

## Extraction (code/agent_a2_24_contact_extract.m)

Given f and D=(u,v): u | (f - v^2) (since (root(u),v) is on the curve),
f - v^2 = u*g4.  Write h3 = v + u*(s x + w); then
   g4 - 2 v (s x + w) - u (s x + w)^2 = kappa * u^2   (deg-4 identity),
5 coeff eqns in (s,w,kappa); solve via Variety.  Verified: residual = 0.

## Validated data (the two known simple Z/24 curves)

Curve A  (r=5, p=-5/2, t=-9/2;  chi_17 simple):
   q3 = x^2 - 435/73 x + 2529/292   (disc = 4608/5329 = (48/73)^2 * 2,
        IRREDUCIBLE over Q -- splits only over Q(sqrt2));
   h3 = 1205280 x^3 - 10769040 x^2 + 31705560 x - 30713580;
   kappa = -1451998172160;  lead(h3)^2 + kappa = 701706240 = lead(f).  OK.

Curve B  (r=1/3, p=-1/9, t=-1;  chi_13 simple):
   A(8) data: d=0 (so Q = x^2), q = (4/9)(x^2 - 4/3 x + 16/27);
   q3 = x^2 - 2/3 x = x*(x - 2/3)  -- SPLITS RATIONALLY;
   h3 = 13122 x^3 + 26244 x^2 - 34992 x + 15552;  kappa = 1377495072.
   => order-3 divisor = P1 + P2 with P1=(0, 64/243), P2=(2/3, 32/243)
      BOTH RATIONAL points.

## Two flavours of rational 3-torsion divisor

 - q3 irreducible (curve A): the two contact points are conjugate over a
   real quadratic field; the class is still rational.
 - q3 split (curve B): the 3-torsion is a sum of two rational points.
This split/non-split dichotomy is a natural way to stratify the search:
the "two rational points of order-3-sum" sublocus may be more tractable.

## Plan

Harvest many Z/24 (r,p,t) samples (code/agent_a2_24_ztors_sample.m, no
2-rank filter), extract contact (q3,h3,kappa) for each, and fit how the
3-torsion condition sits on the A(8) chart.  Then intersect with the
W-split (2-rank-2) locus to parametrize [2,24].
