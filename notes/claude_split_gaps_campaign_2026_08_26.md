# Split-torsion gap campaign (2026-08-26 session)

AVS directive: make a serious attempt at the unrealized geometrically split
torsion groups.  Explicitly requested: (a) the x -> x^2 bielliptic sweep over
the genus-0 torsion families (not just LMFDB seeds), (b) systematic gluing
checks (e.g. 5 (x) 7 for [35] — already realized by HLP2000, row 7436.a.2),
(c) **Q-simple geometrically-split constructions now count**: gluings need
not be defined over Q, factors may live over quadratic fields (conjugate
pairs, Weil restrictions, isogeny classes of E1 x E2 / Res_K(E)).
Early-abort mod-p checks before every exact TorsionSubgroup (split_lab
Funnel does this).  AWS spot use pre-authorized for scale-out.

## Target inventory (as of this session's start)

30 Mazur-product gaps (product/split_torsion_table.md §5):
- (3,3)-route: [56], [2,28], [2,40], [2,56], [10,10], [2,2,20], [2,2,40],
  [2,8,8], [2,2,8,8]
- 2-glue + gain: [2,60], [4,24], [9,9], [12,12], [2,2,2,12], [2,2,2,24],
  [2,2,6,6] (last: sigma-chart exhausted 08-14, thin)
- no 2-/3-route: [42], [72], [84], [90], [2,36], [2,42], [2,72], [3,18],
  [3,36], [6,18], [2,2,30], [2,4,12], [2,4,24], [2,6,12]

NEW target class unlocked by the Q-simple relaxation: groups realized
simple-only that quadratic-conjugate gluings could realize split.  Odd
torsion of E(K), K quadratic, injects into the Jacobian glued from
(E, E^sigma) — Kenku–Momose allows 11- and 13-torsion over quadratic
fields, so **[11] and [13]** (both simple-only: Ogg 353.a, Flynn/Leprevost
349.a) are candidate NEW split groups; also [22]/[26]-adjacent shapes if the
2-part cooperates.  Blocker to derive: E[2] =~ E^sigma[2] as Gal_K-modules
(a twisted 2-congruence condition on the quadratic point of X1(11)/X1(13));
alternatively glue along 3.  NOT yet implemented.

## Engines built and validated this session

1. **Contravariant pencil = X_E^-(3) (the anti-3-congruent family), with
   correct twists built in.**  For U := GenusOneModel(3, E), (P,Q) :=
   Contravariants(U) (Fisher's g1-model machinery in Magma):
   E2(lam) := Jacobian(P + lam*Q) is 3-congruent to E (0/43 trace
   mismatches, every sampled lam) and Genus2Elliptic3(E, E2(lam)) returns a
   gluing for EVERY member; the covariant Hesse pencil U + lam*H is the
   symplectic family and NEVER glues (0/3).  This makes (3,3)-gluing
   searches pure P^1-sweeps — no Frengley surface, no twist alignment, no
   parametrization-of-fibers problems.  (The 2026-08-12 z3 lane failed
   because blind Z(3,2) surface sweeps are torsion-thin: two torsion
   conditions + twist pinning are codim ~3 on rational points.  Fibering
   Z(3,2) over j1 also dead-ends: the genus-0 fibers have essentially no
   small rational points in the (u,v)-chart.)
2. **lane_c33.m** (rewritten): per family member E1(t) (families 7, 8, 10,
   24=[2,4], 28=HLP E_{2,8}), sweep the contravariant pencil, gate E2(lam)
   by exact divisibility ordm | #E2(F_p) at ~6 primes (no twist slack),
   exact elliptic torsion on gate-passers, Genus2Elliptic3 + split_lab
   Funnel on target hits.  ~1500-1900 members/s single-core.
3. **lane_biell_fam.m**: the x -> x^2 sweep over families
   {5,6,7,8,9,10,12,[2,4],[2,6],[2,8]} with shift s: C: y^2 = f_t(x^2+s),
   J ~ E1(t) (2,2)-glued with the reversed-cubic partner E2(t,s).  Gates:
   E2 torsion >= 5 (mod-p gcd then exact), or 2-primary gain possible
   (v2(#E1 #E2 (F_p)) >= v2(|T1_2||T2_2|)+1 at all test primes).  Funnel
   with OddInvs = T1odd x T2odd pinned.  Magma pitfalls hit: EllipticCurve
   REQUIRES monic cubics (non-monic input stalls/errors — always monicize
   y^2=ax^3+bx^2+cx+d as Y^2=X^3+bX^2+acX+a^2d); X1(12) Kubert variant is
   (b,c) = (+,-) signs.
4. **x9join.gp**: deep 2-congruence hash-join on X1(9) (any 2-congruent
   pair of 9-torsion curves 2-glues to J >= [9,9] — the 2-division cubic is
   forcibly irreducible).  ARTIFACT: the diamond <2> acts on the t-line by
   t -> (t-1)/t (order 3), so every field bucket contains the whole
   diamond orbit of each curve — genuine pairs must be j-distinct.  H=80:
   7861 curves, 3931 fields, 0 genuine pairs.  Square-disc (cyclic-cubic)
   members would even SELF-glue to [9,9]: the square-disc locus of X1(9) is
   y^2 = t(t-1)(t^2-t+1)(t^3-6t^2+3t+1), genus 3, no nondegenerate points
   to height 1e4 (X1(7) analogue: genus 2, same emptiness).

## Running (2026-08-26 evening, local box)

- x9join H=800 (bg task btvib41i7) -> product/logs/x9join_H800.log
- biellfam all 10 families TH=SH=24 sequential (bg bo2uwh8hg) ->
  product/logs/biellfam_F*_TH24.log
- c33 all 9 (3,3)-targets TH=16 SH=60 sequential (bg b554dizzx) ->
  product/logs/c33_*_TH16.log

## Honest odds + structured next stages

The (3,3) target loci and the 9x9 2-congruence locus are (twists of)
modular curves of substantial level — likely high genus, so blind height
sweeps are luck-plays (cheap, so run anyway).  The structured follow-ups,
in value order:
1. **[11]/[13] split via quadratic-conjugate gluing** — derivation lane:
   parametrize quadratic points of X1(11) (genus 1) resp. X1(13) (genus 2,
   hyperelliptic fibers), work out the twisted 2- (or 3-) congruence
   condition E[2] =~ E^sigma[2], find the right normal form (palindromic
   sextics y^2 = a x^6 + b x^5 + ... + d^3 a are the shape), or compute the
   glue over K and descend.  A hit = a brand-new split group.
2. **Component analysis of the (3,3) correspondence** for one flagship
   target ([10,10] or [2,2,20]): j(E2(lam)) = j(E_{N2}(m)) fiber product
   over the (t,lam)-surface; look for genus <= 1 components (the
   winning-pattern: parametrize the condition, not the box).
3. **[12,12]/[2,60]/[4,24]/[2,2,2,12]/[2,2,2,24] gain systems** (plan note
   claude_split_2266_1212_plan_2026_08_13.md) at heights >> the H=250
   already swept — C presieve + spot instance.
4. Isogeny-class walks (Richelot + (3,3)) around any new interesting glue.

## Session artifacts

- product/code/lane_c33.m, lane_biell_fam.m, x9join.gp (new lanes)
- this note; logs in product/logs/
- x9join H=800 RESULT: 778997 curves, 389499 fields, 0 genuine pairs, 0 cyclic-cubic members (log product/logs/x9join_H800.log).  The [9,9] 2-glue supply is empty to t-height 800; route parked pending a structural idea.

## RESULT (2026-08-26 late): first split realization of [11]

**New group for the split table: [11]** — previously simple-only (Ogg 353.a).
Witness (minimal model):

    C: y^2 = 1204142x^6 - 5109634x^5 + 31412066x^4 - 65405928x^3
             + 99564424x^2 + 94188680x + 9471400,   J(Q)_tors = [11] exactly.

Construction (lane_qglue.m, fully automated end-to-end, 1 s to rediscover):
X1(11) raw model r^2 - r(s^3-3s^2+4s) + s = 0 ((b,c) = (rs(r-1), s(r-1)));
the fiber s0 = 4/5 gives E over K = Q(sqrt 11) with E(K)_tors = [11], j
irrational, non-CM, E NOT isogenous to E^sigma (17/28 trace mismatches, no
Phi_d relation d <= 13).  Among 1109 quadratic points to s-height 30, s0=4/5
is the ONLY one with E[2] ~ E^sigma[2]; the 2-division field is S3 over K, so
the equivariant iso psi is UNIQUE, hence the sigma-descent condition
psi^sigma = psi^{-1} is automatic and graph(psi) is Gal_Q-stable.  The
HLP/BHLS glue over the splitting field gives an even sextic over K on which
sigma acts by x -> 1/x; matrix Hilbert-90 (N = A + M A^sigma) + scalar
rescale descends it to Q; the twist is pinned by 11 | #J_d(F_p) (d = 2 works,
d = -5 is the twin descent with the 11-torsion in the minus part, torsion
[]).  Verification: product/code/verify_split11.m (asserts exact torsion,
Weil-restriction count identities at 18 split + 22 inert primes, Q-simplicity
via trace mismatch + j irrational; 0.4 s).  Provenance scripts:
`code/split_campaign_2026_08_26/` (probe -> glue -> diag -> check -> descend
-> final -> validate, with a narrative README); the Mestre/IgusaClebsch
reconstruction route FAILS here (the curve
is bielliptic, Aut = V4, so Q-forms are more than quadratic twists and Mestre
returns the wrong form family) — honest Weil descent is required.

## [13]: structural finding (negative for the 2-glue)

On X1(13) (raw model r^3 - r^2 s^4 + 5r^2 s^3 - 9r^2 s^2 + 4r^2 s - 2r^2
- rs^3 + 6rs^2 - 3rs + r - s^3, genus 2, hyperelliptic model
y^2 = 9x^6 - 180x^5 + ... + 12753), EVERY hyperelliptic-fiber quadratic point
has E^sigma isomorphic to a quadratic TWIST of E (the hyperelliptic
involution is a diamond operator; x0 = 2 gives the classical Q(sqrt 17)
13-torsion point).  The unique equivariant psi is then the restriction of the
degree-1 twist isomorphism — Kani-degenerate, so the (2,2)-glue is
decomposable and NO genus-2 Jacobian arises: the 2-glue route to [13] is
dead on the whole fibral family.  Live [13] routes: (i) exceptional
(non-fibral) quadratic points of X1(13), if any; (ii) the chi_d-twisted
(3,3)-glue: an ANTI-isometry E[3] -> E^{(d)}[3] is never twist-induced
(the twist iso is isometric of degree 1), so it would glue nondegenerately;
necessary condition a_p(E) = 0 mod 3 at every chi_d-inert prime — a cheap
strong scan filter over the fibral family.  Not yet implemented.

## Depth of the [11] uniqueness (post-hit scans)

lane_qglue N=11 at s-height 250: **76,093 quadratic points of X1(11), still
exactly one sigma-2-congruent point** (s0 = 4/5, the witness).  The twisted
2-congruence correspondence over X1(11) appears to have essentially a single
accessible rational point — the [11]-split witness is sporadic, not the
smallest member of a visible family.  (Logs qglue11_SH70/SH250.)

## [13]: the chi_d-twisted (3,3)-route is also empty at accessible height

lane_q13c3 SH=120: 17,540 X1(13) fibral quadratic points, ZERO survivors of
the necessary condition (a_p = 0 mod 3 at >= 8 chi_d-inert primes).  Both
structured routes to a split [13] are now closed on the hyperelliptic-fiber
family: the 2-glue is Kani-degenerate identically, and the twisted 3-glue
locus has no points to height 120.  Remaining ideas: exceptional
(non-fibral) quadratic points of X1(13) if any exist (literature: quadratic
points of X1(13) — check); chi_d-twisted glues along N >= 5 (thinner);
entirely different isogeny-extra mechanisms.  [13]-split stays OPEN.

## LMFDB cross-check of the [11] configuration (AVS question 2026-08-26)

The witness's E is NOT in the LMFDB: over 2.2.44.1 = Q(sqrt 11) the DB has
2378 curves reaching conductor norm 500 only, and Norm(cond E) = 6050 =
2*5^2*11^2 (conductor p2 * p5 p5' * p11^2).  The LMFDB's five known
quadratic-11-torsion conjugate pairs (Q(sqrt 2) norm 46, Q(sqrt 17) 172,
Q(sqrt -7) 268, Q(sqrt 13) 828, Q(sqrt -2) 4338 — conjugate-distinct
conductors, so E never isogenous to E^sigma) were all tested:
sigma-2-congruent: NONE; sigma-3-congruent: NONE (trace-mod-3 mismatches
22-32 out of ~60 primes each; scripts `code/split_campaign_2026_08_26/glue_lmfdb11.m`,
`lmfdb11_c3.m`).  The s0 = 4/5 point over Q(sqrt 11) remains the only known
glueable 11-torsion configuration over any quadratic field.

## Complete LMFDB quadratic-torsion glue census (AVS directive, 2026-08-26)

lane_lmfdbglue.m swept ALL 23,358 LMFDB elliptic curves over quadratic
fields with torsion order >= 5 (18,000 quadratic-j + 5,358 rational-j;
6 minutes; data-integrity 0 warnings): sigma-2-congruent pairs 4674,
sigma-3-congruent 1734, chi_d-twisted-3 survivors 48 — but the torsion
distributions are entirely small: SURV2 tops out at [14] (x2), [10], [7];
SURV3 at [18], [15], [14]; SURVT3 all [6]/[2,6].  EVERY curve with odd
torsion part >= 11 in the LMFDB fails all congruence tests against its
conjugate (the five [11]-pairs and the single rational-j [13] curve
included).  Conclusion: no gluable configuration in the entire LMFDB
quadratic range can produce a torsion group outside the 76 realized — the
Q(sqrt 11) [11]-witness (conductor norm 6050, beyond the DB's norm-500
range over that field) is not an accident of DB coverage on the DB side:
the congruence condition is what is rare.  The abundant small-torsion
congruent pairs ([2,4], [6], ...) glue only to realized groups and were
recorded but not funneled.  Log product/logs/lmfdbglue_full.log.

## biellfam final numbers (17-shard fleet, 2026-08-26 close-out)

The x -> x^2 bielliptic family sweep (lane_biell_fam.m) completed over all
ten families {5,6,7,8,9,10,12,[2,4],[2,6],[2,8]} at TH = SH = 24:
5,105,564 (t,s)-pairs, 441 gate-1 candidates (reversed partner with torsion
order >= 5), 1,207,972 gain-gated funnel calls, **zero groups outside the
77 realized** (every funneled glue aborted against KNOWN or resolved
exact-known).  Logs product/logs/biellfam_F*_TH24*.log (per-family shards).
The campaign's box-search phase is closed; the structured frontiers that
remain are listed under "Honest odds" above.
