# [2,24] incompatibility localized: the cubic-coordinate 5-adic/2-adic obstruction — 2026-07-21

Continuation of notes/claude_prod_06_224.md (ranked lane #2). Goal: explain why the
dense mech-A [2,12] family (6,420 exact points) and the dense halvable locus never
intersect ([2,24] = both). Scripts: code/claude_224_kernel_scan.m,
code/claude_224_cubic_probe.m, code/claude_224_cubic_probe2.m (+ logs in results/).

## Findings (empirical, large-sample, zero exceptions)

Setup: odd quintic f5 from the M(12) halving chart (z = (t^2+1)/(2t), a = (1-z^2)/(4(r+1)),
w = 2(r+1)/(1+z)); D = [x+1/w, L r(r+1)/w^3] the order-12 class; mech-A points give a
second rational root k of f5. Halving criterion (validated Cassels/Schaefer x-T):
coordinates c(xP - rho0), c(xP - k) [rational roots] and eps = c(xP - theta) in the
cubic algebra L must all be squares (c = lc(f5), xP = -1/w).

1. s1 = c(xP - rho0) is a square at ALL 6,420 mech-A points. PROVEN symbolically:
   its odd square-class part is z^2 - 1 = ((t^2-1)/(2t))^2 on the t-parametrization —
   identically square. (The rho0-coordinate is trivial by chart construction.)
2. s2 = c(xP - k) is a square at 389/6,420 points; kernels random-looking (no fixed
   class): no obstruction in the second rational coordinate.
3. At ALL 389 both-square points the torsion is still exactly [2,12] — the halving
   failure is ENTIRELY in the cubic coordinate eps.
4. Probing eps at 40 such points: Norm(eps) always a square (product relation),
   eps totally positive at every real embedding (both signatures (1,1) and (3,0)) —
   the obstruction is invisible to norm and signature.
5. MECHANISM (12/12 deeper probes): the ideal (eps)O_L has ODD valuation at a prime
   above 5 (usually the inertia-degree-2 prime of norm 25; once a pair of degree-1
   5-primes), AND eps is a non-square unit at a 2-adic place in every case. Either
   fact alone forces eps notin (L*)^2, i.e. D not halvable.

## Conjecture (theorem-shaped, next to prove)

On the mech-A locus of the M(12) halving surface, with s2 a square, the 5-part of
(eps) is never an ideal square (odd valuation above 5), hence D is never 2-divisible
and [2,24] does NOT occur through this component. Proof route: 5-adic boundary
analysis of the chart (f5 mod 5, position of xP, splitting of the cubic mod 5) per
the component-boundary-analysis methodology; the 2-adic unit obstruction is a
parallel handle. This would upgrade the 6,426-point empirical emptiness to a
structural closure, and redirect the [2,24] hunt to OTHER components (the M(12)
surface's other mechanisms, or entirely different charts).

## Status of other overnight lanes (same session)
- Z/96, Z/160 via Elkies-32: CLOSED (notes/claude_z96_z160_closure.md).
- (8,8) fleet: 550+ bases, gate-passing members appear (~1.5%), all exact J1
  torsion [4,4] so far; continuing.
- 3-gate scan extension to u-height 400: 0 passes through height 300+.
- Multigrade deep sieve: complete, still exactly 3 orbits (X2<=250k, X3<=1e7).

## Final overnight-fleet tallies (2026-07-21, run stopped to free the machine)

(8,8) J1-hunt on Nicholls Lambda_334 (code/claude_88_worker.m + claude_88_drive.sh):
775 bases (m,n) by height, 589 with fibration rank >= 1, 1,020 double-stage-1
members built, 32 passed the (Z/8)^2-embedding gate at 6 primes, all 32 exact
TorsionSubgroup(J1) = [4,4] — the lift wall holds across a 250x larger base
sample than the 2026-07-18 evidence. Summary: data/claude_88_hunt_summary.txt.
Next escalation for (8,8) is the SYMBOLIC lift layer (Kummer coordinates of D_i
over Q(s,t,v), then the specialize-and-profile treatment), not more scanning.

Z/96 sporadic exclusion extended: 3-divisibility gate scan to u-height 400:
0/194,711 members pass (results/claude_z96_height_scan_h400.log).
