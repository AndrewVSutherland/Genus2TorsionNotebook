# Which torsion groups are realized by infinitely many simple Jacobians? The record.

**What this file is.**  The project's census (paper Table 1) lists the 72
finite abelian groups G for which we can document a genus-2 curve X/Q
with geometrically simple Jacobian and J(Q)_tors isomorphic to G, exactly.
A natural companion question, group by group, is: *do we know infinitely
many such curves?*  A column recording the answer (marked with infinity
symbols) lived in Table 1 for a while and was removed from the manuscript
on 2026-08-20 for reasons of space and scope.  This file is the permanent,
self-contained repository record of the answer.  The detailed campaign log
with day-by-day history is
[claude_census_infty_upgrades.md](claude_census_infty_upgrades.md); the
present file states the mathematics and the certificates.

**The bottom line.**  For **39** of the 72 groups we have a certified
infinite family of geometrically simple exact realizations: 21 recorded
before August 2026 (project charts plus literature families), and 18
certified in the campaign of 2026-08-14 to 2026-08-16 documented here.
For **12** groups a positive-dimensional construction is recorded but at
least one leg of the certification is open.  For **21** groups nothing
positive-dimensional is recorded at all.

---

## 1. What we require before writing "infinitely many"

Saying that a group G is realized infinitely often is a three-part claim,
and each part is certified separately.

**(L1) A marked family, generically exact.**  We exhibit a family of
genus-2 curves over a rational base (typically P^1 with parameter t) such
that the Jacobian of the generic fiber carries a *marked* subgroup
isomorphic to G over the function field Q(t) — meaning the torsion
divisor classes are written down algebraically in t, not found fiber by
fiber.  Markedness gives J(Q)_tors containing G for every rational fiber
outside a proper closed locus.  To rule out the generic torsion being
strictly *larger* than G we exhibit at least one fiber whose torsion is
exactly G; since excess torsion is a closed condition, exactness at one
fiber forces generic exactness.

**(L2) Generic geometric simplicity.**  One *strict* fiber suffices: a
fiber C and a prime p of good reduction such that the characteristic
polynomial of Frobenius on C mod p is irreducible of degree 4 and its
root pi satisfies [Q(pi^n) : Q] = 4 for all n up to 12.  The power
condition is what makes the certificate *geometric*: it rules out the
Jacobian becoming isogenous to a product after any base change (an
irreducible charpoly alone only proves simplicity over Q, a mistake this
project has had to retract once before).  Geometric simplicity of one
fiber propagates to the generic fiber.

**(L3) Infinitely many honest fibers.**  The base must be rational (so
rational parameters are infinite in number) and the modulus must be
nonconstant, witnessed by fibers with pairwise distinct G2-invariants.
This clause exists to exclude cheats: a family of quadratic twists of one
fixed curve has infinitely many members but constant G2-invariants, and
proves nothing about the moduli of realizations.  (Twists are structurally
useless here anyway: twisting preserves exactly the rational 2-torsion
and generically kills every element of order > 2.)

In the certificate tables below, "strict p" is the L2 prime for that
fiber; every displayed family also passed the distinct-G2-invariants
check across its listed fibers.  All computations are Magma
(`TorsionSubgroup` on integral models; L-polynomials for strictness);
scripts are `code/claude_census_infty_*.m` with logs under `results/`.

---

## 2. The eighteen groups certified in the August 2026 campaign

### 2.1 Factor-type families: [2], [2,2], [2,2,2], [2,2,2,2]

Rational 2-torsion is read off the factorization of f in y^2 = f(x): on
an odd-degree (degree-5) model the 2-rank equals (number of irreducible
factors) - 1, on an even-degree (degree-6) model it is (number of
factors) - 2.  So prescribing the factor type marks the full 2-elementary
group over the function field, and a free coefficient t moves the modulus.

| Group | Family | exact+strict fiber | strict p |
|---|---|---|---|
| [2]       | y^2 = x(x^5 + x + t)                | t = 2 | 17 |
| [2,2]     | y^2 = x(x-1)(x-2)(x^3 + x + t)      | t = 3 | 11 |
| [2,2,2]   | y^2 = x(x-1)(x-2)(x^2 + x + t)      | t = 3 | 29 |
| [2,2,2,2] | y^2 = x(x-1)(x-2)(x-3)(x-t)         | t = 5 | 11 |

(The [2,2] family needs the even-model rank rule: an earlier attempt
x(x-1)(x^4+x+t) has only three factors, hence 2-rank 1 — that slip and
its correction are recorded in the campaign log.)

### 2.2 Contact families: [3] and [10]

A rational 3-torsion class exists iff the sextic has the contact form
f = h^2 - lambda*Q^3 (the paper's order-3 normal form).  The family

    [3]:   y^2 = (x^3+x+1)^2 - t (x^2+1)^3

therefore has a marked 3-class over Q(t).  Certificate: exact fiber
t = 2, strict p = 13.

For [10], take the plain contact-5 chart at slice b = 1:

    [10]:  y^2 = (1 + a x + x^2)^2 - (a+2)^2 x^5.

The contact identity marks a class of order 5, and the slice has a free
bonus: x = 1 is a root of f *identically* (f(1) = h(1)^2 - (1+a+b)^2
vanishes when b = 1), so [(1,0) - infty] is a marked rational 2-torsion
class, and the marked group contains Z/10.  Certificate: all eight fibers
a = 1, 2, 3, -3, 5, 1/2, 7, -5 are exactly [10], strict primes
17, 17, 13, 11, 23, 11, 19, 11.  (The project's [20]-family lives on the
special slice b = (a^2-1)/2 of the same chart, where an extra halving
exists; the generic slice stays at [10].)

### 2.3 The A(8) chart: [8]

The named chart A(8) (see `named-charts-reference`; formulas as in
`code/agent_a2_24_composite8x3.m`) carries a symbolically verified marked
divisor class of exact order 8.  The slice p = 2, r = 3 with free
parameter t is a one-parameter rational family.  Certificate: exact fiber
t = 2, strict p = 13.

### 2.4 Automatic halving: [4], [2,4], [2,2,4]

The engine for order 4 is a quadratic whose root *is automatically minus
a square*.  Define

    q_{alpha,beta}(x) = x^2 + ((2 alpha beta + 1)/beta^2) x + alpha^2/beta^2.

Its root theta satisfies the identity  -theta = (alpha + beta*theta)^2,
and q_{alpha,beta}(0) = (alpha/beta)^2 is a square.  Consequence: on any
model f = x * (product of such q's), every component of the x-T descent
image of the 2-torsion class e = [(0,0) - infty] is a square (the
finite components are the -theta's; the self-component is the product of
the q(0)'s; cf. Stoll's descent formalism).  So e = 2D for a rational
class D, and D has exact order 4, *by construction, at every fiber*.
Empirical control: 23 of 23 fibers across the three families below show
exactly the predicted group.

**[4]:**  to keep the 2-rank at 1 the quadratic must stay irreducible, so
run the construction over K = Q(sqrt 2) with alpha = 1 + s*sqrt2,
beta = 1, and let Q4 in Q[x] be the norm (the product of q_alpha with its
conjugate).  Family y^2 = x * Q4(x).  All eight fibers
s = 1, 2, 3, -1, -2, 1/2, 5, -3 exact [4]; strict primes
19, 17, 19, 19, 17, 11, 11, 19.

**[2,4]:**  y^2 = x (x^2 + (2t+1)x + t^2)(x^2 + 3x + 1) — two rational
q-quadratics, (alpha,beta) = (t,1) and (1,1), with discriminants 4t+1 and
5 kept nonsquare so the 2-rank is exactly 2.  All eight fibers
t = 3, 5, 7, -2, 1/2, -4, 11, 13 exact [2,4]; strict primes
17, 11, 11, 13, 29, 11, 13, 37.  (The first slice tried had discriminant
25 — a square — and correctly came out as a [2,2,4]-family instead: the
disc condition is load-bearing.)

**[2,2,4]:**  y^2 = x (x+4)(x+9)(x^2 + (2t+1)x + t^2): 2-rank 3, and the
descent components of e at the rational factors are 4 and 9 — squares —
so the halving still goes through.  All seven fibers
t = 1, 3, -2, 5, 1/2, -4, 7 exact [2,2,4]; strict primes
13, 19, 17, 13, 29, 17, 17.

### 2.5 Square-branch and doubling charts: [2,2,2,4], [2,2,4,4], [2,2,2,8], [2,2,2,6]

These four come from the paper's own chart machinery; the families are
explicit rational curves inside the relevant moduli charts.

**[2,2,2,4]:**  on y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2) the class
[(0,0) - infty] is divisible by 2 with an explicit Mumford half (paper,
Lemma "2244"(b)).  Slice (a,b,c,d) = (1,2,3,t):
y^2 = x(x+1)(x+4)(x+9)(x+t^2).  Exact+strict fiber t = 5, p = 59.

**[2,2,4,4]:**  the paper's explicit rational curve inside A(2,2,4,4)
(the Proposition displaying a, b, c, d as polynomials of degree <= 9 in
s, e.g. a = (s-2)(s^2-s+1/2)(s^2-3/2)(s^4+s^2+9/4)).  Marked [2,2,4,4]
is proven there; certificate fiber s = 3, strict p = 37.  (The paper
separately proves rational points on A(2,2,4,4) are Zariski dense.)

**[2,2,2,8]:**  the explicit rational curve on the K3 surface
M(2,2,2,8) from the paper's Theorem:
a = -4t^2(t+1)/(t^2+t+1)^2, b = -t/(t+1), c = 1, d = t on the
square-branch chart.  Exact+strict fiber t = 2, p = 31.

**[2,2,2,6]:**  the M(2,2,2,6) chart (notes/m2226_order6_doubling.md):

    y^2 = x (x+2s^2-sn)(x+2s^2+sm-2sn-mn)(x+2s^2+sm-sn-mn)(2x-mn)(2x+4s^2-4sn-mn)

has six rational linear factors (full marked [2,2,2,2]) and a marked
order-6 class, so [2,2,2,6] is marked over Q(s,m,n); the chart is
rationally parametrized.  Five fibers (s,m,n) = (25,-26,-15), (7,2,1),
(3,-2,5), (11,4,-3), (5,7,2): all exactly [2,2,2,6], strict primes
23, 47, 41, 23, 67, five distinct G2-invariants.

### 2.6 Richelot images: [4,4] and [2,8]

Idea (Filip's): a certified family with full rational 2-torsion has 15
rational Richelot (2,2)-isogenies per fiber, computed by Magma's
`RichelotIsogenousSurfaces` in a deterministic order; the *image* torsion
is source-torsion divided by the kernel (plus possible gains), and since
the isogeny is algebraic in the family parameter, a partner index that
realizes group G at every sampled fiber gives a G-marked family with all
three legs inherited.  Tabulating all 15 partners across fibers:

**[4,4]:**  partner #5 of the [2,2,4,4]-family of 2.5.  Quotient
arithmetic: ((Z/2)^2 x (Z/4)^2) / <e1,e2> = [4,4] when the kernel is the
two free 2-generators.  Fibers s = 3, 5, 7, 9, 11: all exactly [4,4],
strict primes 37, 73, 61, 67, 79, five distinct invariants.  (8 of the
15 partners realize [4,4] at every fiber; #5 is the certified one.)

**[2,8]:**  partner #4 of the [2,2,2,8]-family of 2.5.  Quotient
arithmetic: ((Z/2)^3 x Z/8) / K = [2,8] for K missing the bottom of the
8-chain.  Fibers t = 2, 3, 5, 7: all exactly [2,8], strict primes
31, 73, 59, 83.  (12 of 15 partners realize [2,8] at every fiber.)

Reach of the method was mapped and is worth remembering: second-generation
partners of the [2,8]-family lead only back to [2,2,2,8] (the dual
isogeny), the [2,20]-family's fibers have 2-rank 2 with no usable rational
kernels, and order arithmetic excludes [2,4,4] (all 15 partners of
[2,2,4,4] give only [4,4], [2,2,4], [2,2,2,2]).  Richelot is prime to 3,
so it can never manufacture [3,9]-type structure.

### 2.7 The Elkies 5-torsion chart: [2,10] and [2,2,10]

Elkies (Contemp. Math. 796 (2024), Thm 2.3 / Cor 2.4) parametrizes
genus-2 curves with a typical 5-torsion point by a single quadratic
Q = q2 x^2 + q1 x + q0: the curve

    y^2 = g,   g = 4Q^3 + (x^2 - 6x + 1) Q^2 + 2x(x-1) Q + x^2

carries the 5-torsion class T = [D - K], D = {Q = 0, y = 0}.  A rational
Weierstrass point at x = lam imposes his condition (2.22), which after
his desingularizing substitution becomes the quadratic

    lam^2 q^2 + 2 lam q (2q-1)(q-1) + (q-1)^2 = 0,   q = Q(lam)/lam,

whose lam-discriminant is a square iff q(q-1) is.  That is a conic:
q = t^2/(t^2-1), m = t/(t^2-1), lam = (q-1)(-(2q-1)+2m)/q.  So each
rational parameter t buys one rational Weierstrass root.

**[2,10] (three Weierstrass points).**  Given t1, t2, t3, the three
interpolation conditions Q(lam_i) = q_i lam_i determine Q uniquely, and
g then has three rational Weierstrass roots *by construction*.  The
root-pair classes generate (Z/2)^2, and with T the marked group is
Z/5 x (Z/2)^2 = [2,10].  Validation: 6 of 6 random triples gave exact
[2,10].  Certified slice (t, 2, 3):

| t | torsion | strict p |
|---|---------|----------|
| 5 | [2,10] | 37 |
| 7 | [2,10] | 11 |
| 9 | [2,10] | 29 |
| 11 | [2,10] | 13 |
| 13 | [2,10] | 17 |

Five distinct G2-invariants.  Files: code/claude_census_infty_round6.m,
_round6b.m.

**[2,2,10] (four Weierstrass points).**  A fourth point makes the
interpolation overdetermined — one determinant relation among
t1, ..., t4.  On the generic slice (t1,t2) = (2,3) that relation is a
quintic in (t3,t4) of geometric genus 4, so no family lives there; but
its finitely many rational points produced the *first known examples* of
this group, at (t3,t4) = (-6,6) — exact torsion [2,2,10].  The shape of
the seed (t4 = -t3) is the key: taking the extra pair of Weierstrass
points to be the two roots over a *common* q = q(s) (parameters s, -s)
turns the interpolation into the congruence Qq = q*x + kappa*Cq(x) with
Cq(x) = x^2 q^2 + 2xq(2q-1)(q-1) + (q-1)^2, and eliminating kappa leaves
ONE equation whose honest component, for fixed t2, is a *nodal cubic of
genus 0* — for t2 = 3:

    a^3 + 5a^2 - (5/9) a b^2 + (1/3) b^2 = 0     (a = t1, b = s),

parametrized by lines through the node:  t1 = -(45 + 3 tau^2)/(9 - 5 tau^2),
s = tau * t1  (tau = 3 recovers the seed).  The four Weierstrass roots are
lam1, lam2 and the two roots of Cq, all rational; factor type
[1,1,1,1,2] gives 2-rank 3 and marked (Z/2)^3 x Z/5 = [2,2,10].

| tau | rational W-roots | torsion | strict p |
|-----|------------------|----------|----------|
| 2   | 4 | [2,2,10] | 37 |
| 5   | 4 | [2,2,10] | 13 |
| 7   | 4 | [2,2,10] | 73 |
| 1/2 | 4 | [2,2,10] | 59 |
| 4   | 4 | [2,2,10] | 47 |

Five distinct G2-invariants.  Files: code/claude_census_infty_round7*.m.
This row had no recorded realizations at all before 2026-08-16.

---

## 3. The twenty-one groups with previously recorded families

These were certified (or taken from certified literature) before the
August campaign; formulas and provenance are consolidated in
notes/infinite_families_inventory.md, with fiber certificates in
data/claude_census_family_fibers.txt.  In brief:

- **[7]** — contact-7 chart: h = 1 - (7/2)x + a x^2 + b x^3,
  f = (h^2 + (x-1)^7)/x^2; sampled fibers all exactly [7], strict.
- **[9]** — contact-9 chart: h = 1 - (9/2)x + (63/8)x^2 - (105/16)x^3
  + a x^4, f = (h^2 + (x-1)^9)/x^4.
- **[11]** — two published one-parameter families: Flynn's
  y^2 = x^6 + 2x^5 + (2t+3)x^4 + 2x^3 + (t^2+1)x^2 + 2t(1-t)x + t^2 and
  Daowsud-Schmidt's continued-fraction family; both marked by the class
  of infinity.
- **[6], [12], [2,12]** — the M(12) chart T12 = a x^2 - x + r,
  h = (x-r)(T12+1), f = a x^2 T12 (T12+1): the plain chart has generic
  exact [6]; the rational-Weierstrass subchart a = (1-z^2)/(4(r+1))
  upgrades to [12]; the line a = (1-r)/4 gives [2,12].
- **[2,6]** — contact-6 chart (h cubic, f = h^2 - (x-1)^6).
- **[14], [18]** — rational-root subfamilies of contact-7 and contact-9.
- **[5], [15]** — five-torsion families; [15] is inherited from the
  cyclic [30] family below (see the inventory for the exactness
  discussion).
- **[20], [2,20]** — the contact-5 + 4-torsion family and its extra-2
  loci.
- **[21]** — Leprevost's 1991 one-parameter degree-5 family
  f21 = A21^2 - k21 x^3 (x-1)^2 (polynomials in
  code/z21_leprevost_family_verify.m; marked class of order 21 verified
  symbolically over Q(t)).
- **[22]** — rational-branch subfamilies of the order-11 families.
- **[23]** — the Kuru-Sadek quadratic-order construction, genus-2
  specialization (code/order23_kuru_sadek_check.m).
- **[30]** — the simultaneous contact-5/contact-6 family.
- **[32]** — the reconstructed Elkies order-32 component
  (notes/elkies32_reconstruction.md).
- **[3,3]** — Bruin-Flynn-Testa 2014, Thm 6: rational parametrization of
  the pointwise A(3,3) chart.
- **[4,8]** — the tangent-cover family (best 2-primary family recorded).
- **[2,2,2,10]** — Elkies 2024: the Clebsch-Klein full-level-2 plus
  5-torsion two-parameter family.

---

## 4. The twelve groups with a recorded construction but an open leg

- **[13], [17], [19]** — Leprevost's 1991 C.R. notes construct families
  with 13-, 17-, 19-torsion; the formulas are not in the repository and
  the papers were not freely obtainable during the campaign.  Once
  transcribed, the standard fiber certificates are expected to go
  through without obstruction.
- **[26], [28]** — same situation for the Platonov-Petrunin (Doklady
  2012) constructions.
- **[16]** — the recorded positive-dimensional annotation could not be
  traced to a repository note; provenance must be identified before
  anything can be certified.
- **[3,9]** — only isolated examples.  Richelot-type (2,2)-isogenies can
  never create 3-power structure, so the one systematic tool of the
  campaign provably does not apply; a (3,3)-isogeny analogue
  (Bruin-Flynn-Testa machinery) is the natural next tool.
- **[6,6]** — the contact-6 x contact-6 locus is positive-dimensional
  with one certified good point, but no rational curve on the locus is
  known (leg L3 open).
- **[2,2,8]** — needs 2-rank 3 plus an order-8 class.  The A(8)-chart
  with two extra splittings has codimension >= 2 and produced no seeds;
  the honest route is a second-level halving criterion (halve the
  order-4 class again), which was designed but never derived.  Richelot
  gains from neighboring families were checked and do not occur.
- **[2,4,4]** — needs the pair-class halving criterion (the
  half-divisor Proposition) transplanted to a 2-rank-3 chart; also
  provably NOT reachable by Richelot images from [2,2,4,4] (all 15
  partners tabulated).
- **[2,4,8]** — the one-split tangent subcover has a verified exact hit,
  but the genus/rank of the one-split base curve was never computed, so
  infinitude is unproved (leg L3).
- **[2,2,2,12]** — the record group (order 96).  Its realizations come
  from the symmetric surface S: x1^2+x2^2+x3^2 = x4^2+x5^2,
  x1^4+...+x3^4 = x4^4+x5^4 in P^4, of general type, with three known
  good orbits.  The campaign closed every low-degree channel to
  infinitude: S contains NO nondegenerate line (complete Grassmannian
  enumeration; the only lines are 144 node-pair lines, all inside
  degeneracy hyperplanes) and NO nondegenerate conic (128 distinct
  conic-planes, all degenerate).  Among all 248 hyperplane sections
  through >= 8 of the 36 nodes, exactly two elliptic-curve orbits meet
  the nondegenerate locus, with Jacobians y^2+xy+y = x^3-19x+26
  (torsion [2,6]) and y^2+xy+y = x^3+x^2-10x-10 (torsion [2,4]) — and
  BOTH have Mordell-Weil rank 0, their rational points being exactly
  the nodes they pass through.  Infinitude for this group is genuinely
  open and, on this evidence, may well be false (Bombieri-Lang
  pessimism: the three known orbits could be close to all there is).

---

## 5. The twenty-one groups with nothing positive-dimensional recorded

[24], [25], [27], [29], [31] (all known examples have RM), [33], [34],
[36], [39], [40], [2,14], [2,16], [2,18], [2,22] (RM), [2,26], [2,28],
[3,6], [2,2,6], [2,2,12], [2,2,14], [2,2,20].

These have one or a handful of isolated examples (databases, Leprevost,
Elkies 2002, Platonov-Petrunin, Howe, this project) and no family.

---

## 6. Methods notebook (what transfers to future rows)

1. **Factor types** buy the 2-elementary part for free; the odd/even
   rank rules are load-bearing and were mis-stated once (see 2.1).
2. **The automatic-halving identity** (2.4) manufactures order-4 classes
   with no per-fiber conditions; iterating it one more level is the
   designed-but-unfinished route to [2,2,8].
3. **Richelot partner tabulation** (2.6) converts one certified
   full-2-torsion family into up to 15 quotient families at zero search
   cost; its reach is fully mapped and its limits are provable (prime to
   3; order arithmetic).
4. **The Elkies chart** (2.7) is a machine for 5-part groups: each
   rational Weierstrass point costs one conic parameter, and the
   antidiagonal congruence trick collapsed an overdetermined locus from
   genus 4 to genus 0.  Any future group with a 5-part should try it
   first.
5. **Twists never help** with L3 (constant G2-invariants), and they
   destroy all torsion of order > 2 generically.
6. A negative worth keeping: on the [2,2,2,12] surface, *every* rational
   curve of degree <= 2 and every nodal-section elliptic curve is now
   excluded (rank 0 or degenerate).  Any future attack must find
   higher-degree or non-nodal curves on S, or a different construction
   entirely.

## 7. Provenance

- The paper's last pre-deletion column state (21 infinity / 30
  parenthesized / 21 dash — predating the campaign's 18 upgrades) is in
  an earlier revision of the paper source (working-repo git history).
- Campaign proofs and history: notes/claude_census_infty_upgrades.md.
- Scripts: code/claude_census_infty_certify*.m (round 1),
  _round2b/c.m (automatic halving), _round3*.m ([10] and loci),
  _round4*.m / _round5.m (Richelot), _round6*.m ([2,10]),
  _round7*.m ([2,2,10]); surface hunt: code/claude_22212_lines_on_S.m,
  _curves_on_S.m, _conic_screen.m, _nodal_sections*.m,
  _elliptic_section.m.  Logs under results/ with matching names.
- Family inventory for Section 3: notes/infinite_families_inventory.md;
  fiber certificates: data/claude_census_family_fibers.txt.
