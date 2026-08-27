# prod #10: (2,6,6), order 72 — BFS reality check, bielliptic locus B(a,b), (b,v)-cover scan

Date: 2026-07-18. Production campaign step. All computation in
`/tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/prod10_266/`
(reusing the aborted-launch state found there: H=300 extra-root sieve + deep phase were
already complete; `bell.m`, `invol.m`, `a266_slice.m` were written but unrun).
Single-thread budget (briefly 2 during the sweep, dropped back when load > 24).

## Strategy recap (3 lines)

(A) Reality-check the Bruin–Flynn–Shnidman sqrt3-RM escape route against what (2,6,6)
actually needs; it fails structurally, so pivot to imposing extra-2 on the cubic-contact
(6,6)-core cover of the contact-6 chart, swept in the (b,v)-chart where the second
3-class is built in. (B) Derive and certify the bielliptic locus B(a,b) and convert the
empirical split-confinement into a precise interpolation statement. (C) Scan the pivot
chart at production height with B(a,b) evaluated in-loop and jackpot escalation off-B.

## (A) BFS sqrt3-RM reality check — route CLOSED as a (3,3)-substrate, pivot justified

What (2,6,6) = C2 x C6 x C6 needs: 2-rank 3 (quintic factor type [1,1,1,2] or finer;
sextic needs 4 rational Weierstrass points) AND two independent rational 3-torsion
classes (3-rank 2).

Facts from arXiv:2102.04319 (Bruin–Flynn–Shnidman, local copy scratchpad/full2.txt):

1. Their surface H_3 = P^2_{(a:b:c)} \ Delta parameterizes PP abelian surfaces with RM
   by Z[sqrt3] and full sqrt3-level structure eps: Z/3 x mu_3 = J[sqrt3]. The curve is
   y^2 = G1(x)^2 + lambda1*H1(x)^3 — itself a CUBIC-CONTACT model. D1 (from the H1
   roots on y = G1) is a rational point of order 3; D2 (from the H2 roots on
   y = G2/sqrt(-3)) has order 3 only over Q(sqrt-3) and generates the mu_3.
2. Hence J[sqrt3](Q) = <D1> = Z/3 EXACTLY, always: mu_3(Q) = 0. The level structure
   supplies precisely ONE rational 3-class, never two. The (-3)-quadratic twist swaps
   the roles (their Prop 3.5): Z/3 and mu_3 are never simultaneously rational.
3. There is NO second kernel "J[sqrt3']": 3 = (sqrt3)^2 ramifies in Z[sqrt3], (sqrt3)
   is the unique prime above 3 (sqrt3 * unit exhausts the associates), so the only
   intermediate subgroup scheme between 0 and J[3] from the RM is J[sqrt3] itself.
4. A second independent rational 3-class P2 must therefore lie in J[3] \ J[sqrt3] with
   sqrt3*P2 in J[sqrt3](Q) \ {0} = {+-D1} (if sqrt3*P2 = 0 then P2 in <D1>). So
   P2 in sqrt3^{-1}(D1), a torsor under J[sqrt3] = Z/3 x mu_3; its class lives in
   H^1(Q, Z/3) x H^1(Q, mu_3) = (cyclic cubic classes) x Q*/(Q*)^3, and BOTH components
   must vanish: "D1 is sqrt3-divisible in J(Q)". This is a 9-fold division-cover
   condition on H_3 — an arithmetic accident of exactly the same nature as the cubic
   contact accident on our contact-6 chart, not a free structure. (Weil-pairing check:
   e_3(D1,P2)^2 = 1 automatically by Rosati-self-adjointness of sqrt3, so isotropy is
   no obstruction — the torsor condition is the whole cost.)
5. On top, the BFS sextic generically has J[2](Q) = 0; 2-rank 3 costs a further
   4-rational-Weierstrass-points cover (degree up to 360 over the S6-generic surface;
   the built-in (1+sqrt3) Richelot kernel is rational only as a subgroup, not
   pointwise).

Empirical grounding (`bfscheck.m`): the paper's worked example C_{1,2,-1}:
y^2 = 8x^5-3x^4-2x^3-7x^2+4x+20 has factor type [1,4], torsion exactly [6]
(3-rank 1, 2-rank 1) and a STRICT simplicity certificate already at p = 11 —
simple as expected for RM, but with a single 3-class, exactly as the theory above
predicts.

Verdict: BFS substrate = one free 3-class + (division-torsor accident) + (2-rank-3
accident) on a 2-dim base. Our contact-6 chart = one free 3-class + free 2-rank 2 +
(cubic-contact accident) + (one split-quadratic half-accident). The BFS route is
strictly MORE expensive for (2,6,6); its only advantage is generic geometric
simplicity. PIVOT: sweep the cubic-contact cover of the contact-6 chart in the (b,v)
chart (second 3-class built in fiberwise) and let the split quadratic be the accident —
new height region, same locus, complementary to the (eps,r,b) extra-root sieve.

## (B) The bielliptic locus B(a,b) — derived, certified, interpolates all 10 points

Chart: f = h6^2 - (x-1)^6, h6 = 1 + a x + b x^2 + x^3, f = x*q1*q2.
Method (`bell.m`): sextic form F = Z^6 f(X/Z); Mobius involutions sigma(X:Z) =
(pX+qZ : X-pZ) (all involutions not fixing infinity; the affine chart sigma = t-x
separately); rank-1 minors of (F, F o sigma), saturate by det = -p^2-q, eliminate to
(a,b). Result: principal ideal of degree 20 which factors as

  (a-b) * (a+b+2) * C2 * Q2 * Q3 * C3 * C3s * (deg-6 factor),

with (certificates via GeometricAutomorphismGroup at generic rational points,
`valslice.m`, `bellcheck2.m`, `q2check.m`):

| component | equation | status | generic structure |
|---|---|---|---|
| diag | a = b | bielliptic C2xC2 (invol x -> 1/x) | torsion [2,6], ft [1,2,2] |
| C2 | (a+3)(b+3) = 4 | bielliptic C2xC2, invol x -> c/x, c = (a+3)/2 | torsion [2,6]; sporadic [6,6] |
| Q2 | a^2+ab-3a+b^2-3b-9 = 0 | GeomAut order 12 (D6) — split | torsion [2,6] |
| Q3 | a^2+(5/2)ab+9a+b^2+9b+27/2 = 0 | bielliptic C2xC2 | quintic FULLY SPLIT [1,1,1,1,1], torsion [2,2,2,6] |
| C3 (+swap C3s) | 3a^2b+9a^2+3ab^2+20ab+45a+b^3+6b^2+36b+53 = 0 | bielliptic C2xC2 | extra root BUILT IN: ft [1,1,1,2], torsion [2,2,6] |
| a+b+2 | — | SPURIOUS: disc f has factor (a+b+2)^6 | singular curves |
| deg-6 | see bell.log | unclassified, carries no known point | — |

disc(f) ~ (a+3)^2 (a+b+2)^6 * disc(q1) * disc(q2), disc(q1) = a^2-6a-8b-15,
disc(q2) = b^2-6b-8a-15. Rational parametrizations: C2: a = 2t-3, b = 2/t-3;
Q2: a = -3+6/(1+m+m^2), b = 3+6m/(1+m+m^2); C3 is a NODAL RATIONAL cubic, node
(-11/3, 3) (= the unique GeomAut-order-8 point found), parametrized by
s(m) = -(4m^2+16m+18)/(m^3+3m^2+3m), a = -11/3+s, b = 3+ms.

INTERPOLATION (the precise confinement statement): all 10 known [2,6,6] chart points
(H<=300 census) lie exactly on C2 u C3 u C3s — K1,K2,N1b,N2,N4 on C3; K3,N1a,N3,N5 on
C2 — and the SIMPLE (6,6) control point (133/39,-7/13) lies on NO component. Rational
involutions listed per class in invol.log. New H=300 classes certified split:
N4 = (-13673/16275,-251/217), model 244125000x^5-1015952475x^4+1585513050x^3
-1099734299x^2+286049400x, GeomAut C2xC2, new Qbar-class D; N5 = (-835/289, 265/8),
model 386201104x^5+5754193745x^4-905578032x^3+318571280x^2+1183744x, GeomAut C2xC2,
Qbar-class B. Data: `data/claude_prod_10_266_bielliptic_locus.txt`,
`data/claude_prod_10_266_split_models.txt`.

CONJECTURE (2,6,6)-confinement, now precise: every rational point of the [2,6,6]
locus of the contact-6 chart lies on B0 = C2 * C3 * C3s = 0. Verified for all points
to extra-root height 300 and (so far) all points of the (b,v)-cover scan. Since every
component of B0 is certified geometrically split, the conjecture implies: NO
geometrically simple [2,6,6] on the contact-6 chart at any height.

Structural bonus explaining the confinement mechanism: on C3 the extra-2 condition is
FREE (ft [1,1,1,2] identically, generic [2,2,6]), so on C3 only the 3-contact accident
is needed; on C2 the involution x -> c/x is rational with c = (a+3)/2 for ALL points
(no square condition), making C2 the cheapest bielliptic family — both are rational
curves, positive-rank elliptic fibrations over them produce the recurring accidents.

## (B') MAIN STRUCTURAL RESULT: the [2,6,6] census is THREE INFINITE RATIONAL FAMILIES

Discovered from the (b,v)-sweep data + fiber probes, then verified exactly; this
explains the entire census and the "locus keeps acquiring rational points" phenomenon.

1. C2-BASE FAMILY: a = 2v^2-3, b = 2/v^2-3 (C2 point with parameter t = v^2 = the
   involution constant c of sigma(x) = c/x). PROVEN over Q(v) (`c2sym.m`, `c2sym2.m`):
   the contact system has the unique rational branch
     U(v) = (2v^4-v^3-4v^2-v+2)/v,  M(v) = ((v-2)(v-1/2)(v+1))^2 (a square identically),
   so a rational sigma-stable cubic contact q = x^2+U(v)x+v^2 exists for EVERY v:
   generic torsion >= [6,6], all split (on C2). Extra-root discriminants factor as
   disc(q1) ~ (v^2-1)^2(v^2-4), disc(q2) ~ -(v^2-1)^2(4v^2-1), giving two RATIONAL
   [2,6,6] branches:
     A: v = k+1/k  and  B: v = (1-k^2)/(2(1+k^2))  ==> EXACT torsion [2,6,6].
   Known C2-side points: K3 (k=2 A), N1a (k=3 A / k=2 B), N3 (k=1/5 B), N5 (k=3/5 B).
   New members verified at k=5, 7/2 (A), 4/7 (B). Data:
   `data/claude_prod_10_266_c2_family.txt`.
2. C3-FAMILY: m(j) = -3j^2/(2j^2+1) on the nodal-cubic parametrization of C3
   (condition: -m/(2m+3) is a square), i.e. a(j) = (-j^6+3j^4+3j^2+2)/(j^6+j^4+j^2),
   b(j) = -(j^4+5j^2+3)/(j^4+j^2+1). PROVEN over Q(j) (`c3symj4.m`): TWO rational
   sigma-stable contacts, at v = (j^2+1)/(j(j-+1)), with U(j) = -(2j^4+3j^2+-j+2)/
   (j^4-+j^3) and multiplier M = L(j)^2 a square IDENTICALLY,
   L(j) = (j^4 +- j^3/2 + 3j^2/2 +- j/2 + 1)/(j^4 -+ j^3); extra-2 free on C3 ==>
   torsion containing [2,6,6] along the whole family (EXACT [2,6,6] verified at 11
   members). Known points: j = 1/2 (K1), 2 (K2), 3 (N1b), 4 (N2), 5 (N4), 1/3
   (class-A pt (1727/91,-289/91)); new members at j = 6, 7, 2/5, 8/3. Controls off
   the family have NO rational contact (`c3probe.m`, `c3family.m`). Data:
   `data/claude_prod_10_266_c3_family.txt`.
3. Swap images (a <-> b) of both.

Refined confinement conjecture: the [2,6,6] locus of the contact-6 chart is EXACTLY
these rational families (+ swaps) — an infinite, fully structured, entirely SPLIT
census. Any rational point off them (in particular off B0) would be the jackpot; none
found (sweeps below).

## (C) Production scan of the pivot chart (b,v), height 32

`slice_prod.m` (extends tier2 `a266_slice.m`): for each (b,v), solve the cubic-contact
system {F2 = F3 = 0 after linear elimination of a via F1} in (M,U) (0-dim saturated
ideal), keep rational solutions with M = L^2 square; reconstruct a; REQUIRE factor type
[1,1,1,2] (extra-2); classify against B(a,b) in-loop; [2,6,6] structure + 144-filter at
good p <= 67; strict simplicity pre-certificate p <= 97; inline exact phase
(TorsionSubgroup + contact-point orders); OFFBELL/JACKPOT escalation for any candidate
off B0. Validation: reproduces the simple (6,6) at (b,v) = (-7/13, 5/2) with
a = 133/39, L = 29/16, U = -9/4; v -> -v is NOT redundant (sign is pinned by
h3(0) = v^3/L, both signs swept). Smoke H=6: found C2-point (-19/9, 3/2) with contacts
at v = +-2/3; its curve is a NEW SPLIT [6,6] (banked, see data file). Prior coverage
(tier2): full (b,v) to height 16 (only the known nonsimple curve). This run: height 32
(params 1295 per axis, ~1.67M pairs), chunked by b-index for resume.

RESULTS — COMPLETE at height 32: 804,246 + 113,784 + 756,405 = 1,674,435 pairs
= 1295 x 1293 exactly (full coverage, chunks c1/c2/c2b, ~80 CPU-min total).
23 unique (a,b) points carry a rational cubic contact at (b,v)-height <= 32.
COMPLETE CLASSIFICATION (`cls.gp` + final two by hand):
  - 18 on the C2-base family t = v^2 (as the symbolic theorem predicts; the last two,
    (-115/81,-15/32) v=8/9 and (-235/121,25/32) v=8/11, have t = 64/81, 64/121),
  - 2 on C3 (= family j = 1/2, 2: the K1, K2 census points, found through their
    predicted contact data v = 5/3,-5 and 5/2,5/6),
  - 2 on Q2, the D6 split locus: NEW split [6,6] curve (11/7, -19/7) (GeomAut 12,
    torsion [6,6]; note the in-loop bell tagger omitted Q2 — post-hoc classified),
  - 1 OFF-EVERYTHING: exactly the known SIMPLE (6,6) point (133/39, -7/13), which
    turns out to carry THREE rational contacts (v = 5/2, -5/4, 5/6), matching its
    three non-built-in 3-subgroups.
[2,6,6] candidates (ft [1,1,1,2] + finite checks): only the family points K1, K2,
K3-swap; all HIT266 with torsion [2,6,6] exact, bell-classified, simple=false.
ZERO off-locus candidates, ZERO simplicity certificates, NO JACKPOT.

Interpretation: at (b,v)-height <= 32 the rational points of the cubic-contact cover
are EXACTLY (split-Aut locus points with forced contacts) + (the one known simple
(6,6)). The hunt for simple [2,6,6] on this chart = the hunt for a second
off-split-locus point of this cover that ALSO has an extra root; none exists in range.

## Files

Scratchpad prod10_266/: bell.m/.log (B(a,b) elimination), invol.m/.log (per-class
rational involutions), valslice.m (validation + interpolation + component certs),
bellcheck2.m (disc factorization, generic torsion, N4/N5 certs), q2check.m (Q2 = D6
locus), c3rat.m (C3 nodal), slice_prod.m + slice_prod_c1.log/c2.log (production sweep),
deep_h300.log + sieve_h300.txt (H=300 extra-root census, from aborted launch).
Repo data: data/claude_prod_10_266_bielliptic_locus.txt,
data/claude_prod_10_266_split_models.txt.

## Resume state

- Extra-root chart: rigorous negative now to parameter height 300 (26 parametrizations,
  5 Qbar-classes K,A,B,C,D, all [2,6,6], all split, zero strict certificates).
  Resume: `./sieve266 450 tables_exact_p67.txt out.txt` then deep.m (a266v scripts).
- (b,v)-cover: COMPLETE at height 32 (chunks: c1 = b-idx 1..622; c2 = 623..710;
  c2b = 711..1295; logs slice_prod_c1/c2/c2b.log). Extension to height H:
  `nice -n 10 magma -b height:=H ilo:=1 ihi:=<n> slice_prod.m > log 2>&1` (chunk by
  b-index; BDONE lines are the checkpoints). NOTE for any rerun: add Q2
  (a^2+ab-3a+b^2-3b-9) to BellTags in slice_prod.m — it was certified split (D6)
  only after the sweep launched, so off-B classification was done post-hoc (cls.gp).
- Next-step spec (highest value, cheap): restrict the 3-contact cover to the RATIONAL
  curves C2 (a = 2t-3, b = 2/t-3) and C3 (param above) — the [2,6,6]-on-B0 locus
  becomes rational points of an explicit COVER CURVE over t (resp. m); computing its
  genus/rank would either mass-produce split [2,6,6] or prove finiteness, upgrading the
  confinement conjecture to a theorem on B0. Off B0, local solvability of the
  3-contact cover restricted to the chart-minus-B0 is the remaining open door.
- BFS route: closed as a substrate for the SECOND 3-class (this file, part A). If ever
  revived, the object to construct is the sqrt3-division cover {x in J : sqrt3 x = D1}
  over H_3, NOT a second kernel.
