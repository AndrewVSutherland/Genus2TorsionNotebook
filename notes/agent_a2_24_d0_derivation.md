# Contact-3 derivation on the d=0 slice (rigorous, negative for [2,24])

Goal: use the validated cubic-contact 3-torsion (f = h3^2 + kappa*q3^3,
see notes/agent_a2_24_contact.md) to parametrize simple Z/24 curves on
the A(8) chart and intersect with the W-split (2-rank-2) locus to reach
[2,24].  Carried out on the tractable d=0 slice; the outcome rigorously
explains why [2,24] resists this route.

## The d=0 slice
d = e+2p-r^2 = 0  <=>  p = r(r+t)/2.  Then Q=x^2, f = q*(x^4+q), and x=0
is a free rational point (f(0)=q(0)^2=C^2).  Curve B (r=1/3,p=-1/9,t=-1),
one of the two known simple Z/24, lives here.

## Split-through-0 contact and its spurious branch
Ansatz q3 = x*(x-beta) (a 3-torsion point at x=0): then
q3^3 = x^6 - 3b x^5 + 3b^2 x^4 - b^3 x^3 has only degrees 3..6, so the low
coefficients of f = h3^2 + kappa q3^3 pin h3 rationally
(l=C, j=f1/2l, n=(f2-j^2)/2l), leaving (m,kappa,beta) and 4 equations;
kappa=f6-m^2 (x^6), beta from x^5, and x^4,x^3 give Eq4(m),Eq3(m).

Resultant_m(Eq4,Eq3)=Phi(r,t) factors; the deg-6 component through B is a
rank-1 elliptic curve y^2=x^3+x^2-4x -- BUT this is SPURIOUS: it is the
kappa=0 locus (f a perfect square), giving torsion [8], not [24].  Every
generated point except B had torsion [8].  (code/agent_a2_24_d0_derive.m,
d0_family, d0_harvest.)

## Rigorous statement (code/agent_a2_24_d0_saturate.m)
Enforcing kappa<>0 (Rabinowitsch w*kappa-1), the ideal
<Eq4,Eq3,w*kappa-1> over Q(r,t) collapses to (1): genuine rational
3-torsion is NOT cut out by any hypersurface in (r,t).  It is a covering
condition -- the genuine Z/24 are rational points of the cover
   X = { (r,t,m) : Eq4=Eq3=0, kappa=f6-m^2 <> 0 },
a space curve, with B at (r,t,m)=(1/3,-1,2/9) (m0=2/9, kappa=32/81,
beta=2/3, divisor of order exactly 3 verified).

## Genus and rank of the genuine cover X
Eliminating r (or t) gives the genuine component: a degree-12 plane curve
of GENUS 1 (both projections agree).  (code/agent_a2_24_d0_cover_rank.m.)
Magma's EllipticCurve/Jacobian on the deg-12 model is very slow, so the
rank was settled by direct rational-point search
(code/agent_a2_24_d0_cover_ptsearch.m):

  rational points on the (r,m) model, height <= 40:  ONLY 4
    (0,0),(1,0),(2,0)  -- m=0, degenerate (h3 not a genuine cubic)
    (1/3,2/9)          -- B, the sole genuine point.

A bounded 4-point set => RANK 0.  Independently, the complete d=0 (r,t)
scan of 6.3M curves (code/agent_a2_24_d0_fast.m, H=45) found ONLY B.

## Conclusion
On the d=0 slice the genuine simple-Z/24 cover X is a genus-1, rank-0
curve: B is essentially the only genuine simple Z/24 there, and there is
NO positive-dimensional / infinite family to intersect with W-split.  So
the contact route does not reach [2,24] via d=0.  More broadly this
pins down WHY [2,24] is hard: rational 3-torsion is a covering (thin-set)
condition, and even the most tractable slice's genuine cover has rank 0.

## Open
The general (non-d=0) contact (q3=x^2+Ux+V, irreducible allowed, as in
curve A) gives a higher-dimensional genuine cover whose rank/rationality
is not settled here.  That is the remaining door for the contact route.
