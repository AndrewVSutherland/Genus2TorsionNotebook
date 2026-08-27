# Sibling hunt summary — 2026-07-18 evening: SECOND (2,2,2,12) curve found

Four-lane hunt for more (2,2,2,12) curves (notes/claude_sib_{A,B,C,D}_*.md). Wall ~66 min.

## Addendum 2026-07-20: curve #3 found on the ABC surface

The two homogeneous square conditions supplied after this campaign factor as
\[
F=(A^2-B^2)(B^2-C^2)D,\qquad G=F+D^3,\qquad D=A^2-B^2+C^2.
\]
A 16 MB modular block sieve exhaustively searched primitive ABC-height
through 20000.  It found the first new orbit at
\[
(A,B,C)=(4165,19661,5364),
\]
with canonical key
\([16660,21456,78793]\mid[27593,78644]\).  A smaller curve-parameter chart is
\[
(s,m,n)=(1500518600,253638081,138777800).
\]
Magma verifies torsion exactly [2,2,2,12], four strict geometric-simplicity
certificates at \(p=37,127,131,179\), and G2-invariants distinct from curves
#1 and #2.  All 12 symmetry charts have the same G2-invariants.

The same calculation identifies the birational main surface as a nodal
\((2,4)\) complete intersection in \(\mathbf P^4\), hence a surface of
general type with \(\kappa=2,K^2=8,p_g=5,q=0\).  Full details:
notes/22212_abc_surface_search.md; curve certificate:
data/22212_abc_curve3.txt; verifier: code/verify_22212_abc_hit.m.

### Addendum 2026-07-31: curve #3 independently re-verified (fresh session)

Curve #3 is now at the same verification standard as curve #2:

- `code/claude_sib_curve3_verify.m` rebuilds the curve by a DIFFERENT code
  path — direct chart \((A,B,C)=(4165,19661,5364)\) → symmetric-model point
  \([x:y:z:u:v]=[16660,21456,78793,78644,27593]\) (both symmetric equations
  verified; this is exactly the canonical orbit representative) → sextic
  \(W^2=T\prod(T-\alpha^2)\).  Its reduced minimal model matches the recorded
  one in data/22212_abc_curve3.txt verbatim; torsion exactly [2,2,2,12],
  order 96; three strict root-power certificates at FRESH primes
  \(p=193,197,211\) (disjoint from the original 37,127,131,179); G2-invariants
  equal to the quintic discovery model and distinct from curves #1, #2.
  Log: results/claude_sib_curve3_verify.log.
- `code/claude_sib_curve3_lombardo.sage`: Sage certifies
  `geometric_endomorphism_algebra_is_field(B=100) = True` and
  `geometric_endomorphism_ring_is_ZZ(B=100) = True`, i.e.
  \(\mathrm{End}(J_{\overline{\mathbf Q}})=\mathbf Z\) — an endomorphism-based
  simplicity certificate independent of the Frobenius criterion.
  Log: results/claude_sib_curve3_lombardo.log.

## HIT: curve #2, independently found by two lanes, triple-verified

(s,m,n) = (2208, -8303, -7200) on M(2,2,2,6); from (u,g) = (-23/75, 459/23), i.e.
rho' = -95^2/3519 (lane A) = law-sweep point (u,rho') = (-23/75, -9025/3519) (lane C).
  y^2 = prod(A_i + B_i x), A = [1,1,1,2,2],
  B = [25648128, -36568896, -52466496, -59781600, 23309856]
Minimal model: y^2+(x^2+x)y = 36x^6+36750x^5-462983772x^4-301623595823x^3
  +1518598238654317x^2+397058962729817115x-1282993930035013443975
- TorsionSubgroup EXACTLY [2,2,2,12]; simplicity certificates p = 71, 103, 127, 137
  (chi and chi^12 irreducible deg 4); verified on multiple representations by lanes A and C
  and INDEPENDENTLY RE-VERIFIED by the orchestrator (code/claude_sib_curve2_verify.m).
- G2-invariants distinct from curve #1; far smaller parameters ((2208,...) vs (336396,...)).
- 12 pencil representations known (data/claude_sib_t5surf_hitclassify.txt).

## Structural theorems (the hunt is now a theory)
- (u,g) chart of the hit surface S: rho' = (2gu+(u-1)(g^2+1))/(2g*q(u)); on S only THREE
  independent conditions (V'3 dependent, proven), each a QUADRATIC in g. (Lane A)
- PROVEN LAW (lane B/C): W1*W2 == rho'(rho'-1) mod squares => DEAD-member certificate;
  live members are exactly rho' = q^2/(q^2-p^2) (up to sigma) — proves the +-square
  numerator law, kills the "odd" and "240" parts of the old empirical law.
- PROVEN involutions (lane C): tau: (u,rho') -> ((4rho'u-3rho'-1)/(4rho'(u-1)), rho') and
  sigma: (u,rho') -> ((4u-3)/(4u-4), 1-rho'); <sigma,tau> = (Z/2)^2; on S sigma is
  g -> 1/g-type swap. Fixed loci of tau KILLED 3-adically (no infinite family there).
- g-RIGIDITY (lane A): all 24 known hit points lie on FOUR g-fibers +-725/288 (curve 1),
  +-459/23 (curve 2): one special fiber pair = one curve; exactly 6 points per fiber to
  u-height 1e5; u-numerators form 3-term APs ({97,133,169} diff 36).
- Rank signal: special fibers have E2(g) rank 5-6 vs 1-2 generic; necessary-flavored but
  not sufficient (control fiber +-67/22 rank>=5 is hitless to u-height 30000).
- Lane D scans: known curve unique among all members |rn|,rd<=300 to u-height 3000,
  +-odd^2/240 members to numerator 999^2 at u-height 2e5, quintic chart d<=2100
  (unconditional), MW lattices to exponential heights on tasked members. Curve #1 now has
  >=4 representations (two at u-heights 7105, 13872 found by the new C_rho' genus-1 scan).
- Isogeny neighborhood of curve #1: all 17 two-power-isogenous Jacobians have torsion
  order exactly 24 ([2,12] x12, [2,2,6] x5) — 96->24 uniform degrade, no new group;
  twists |D|<=30 all exactly [2,2,2,2]. (Lane D)

## Where curve #3 will come from (ranked)
1. Lane A g-height extension (./t5surf 150 16000 6, then 32000): every new special |g|
   is a new curve; u-height 150 suffices for first contact empirically.
2. Triple-rank filter over g (E1,E2,E4 ranks at candidate g, GRH bounds) to pre-locate
   special fibers cheaply.
3. Lane B smallmem pipeline over fresh live members rho' = q^2/(q^2-p^2), small (p,q).
4. Lane C rational-curve hunt on the triple-conic cover of the (v,t)-plane — the
  remaining infinite-family route (tau-fixed loci are dead; this one is open).
