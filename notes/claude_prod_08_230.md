# prod #8: (2,30) production — infinite simple [15] family + five simple [30] curves on Elkies' A_1(5)

Production session on target (2,30) (order 60), continuing
`notes/claude_top10_09_230.md` (6x5 combination GO). Strategy recap (3 lines):
(A) settle the flagged [5,10] member at q=(3/5,1/2,3/5); (B) mod-p prescan of
the {T5,T3} cubic-contact locus on the q2=-1/4 surface; (C) impose the 3-part
algebraically and cut towards (2,30).

Workdir (all scripts/logs):
`/tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/prod230/`

## Headline results

1. **Task A settled: the [5,10] member is SPLIT** (not a (5,10) realization).
   Integral model `y^2 = 3060x^6 + 1500x^5 + 8305x^4 + 8305x^2 - 1500x + 3060`,
   exact torsion [5,10]; no simplicity cert for p <= 500 (all 91 L-polys
   biquadratic type [2]); GeometricAutomorphismGroup <8,3>; Q-rational
   involution `(x,y) -> (-1/x, y/x^3)`; genus-1 quotients over Q are
   **LMFDB 66.c3 and 66.c4** (same isogeny class! conductor 66, torsion Z/10,
   rank 0). So J ~ E x E' with E,E' isogenous rank-0 curves — proven split via
   Q-rational quotient maps. (Script `taskA_510.m`, output `taskA_510.out`.)

2. **Task B done (previous launch): the {T5,T3} locus is live at every prime.**
   `contact_scan.c`: on the q2=-1/4 surface at p=7..23, per-fiber solvability of
   `h3^2 - F = a3^2 k^3` over F_p EXACTLY coincides with 3 | #J(F_p)
   (mismatch=0 at all primes), and per-fiber solution counts are always in
   {0,1,4,13,40} = (#J[3](F_p)-1)/2 — the contact system captures precisely the
   +-pairs of nontrivial 3-torsion classes. Live (contact AND 2-rank>=2) fibers
   exist at every p (e.g. (6,2) at p=7). Outputs `contact_p7.out`, `contact_p11_23.out`.

3. **NEW: pure-3 sieve finds the missed stratum.** The old level-B sieve
   demanded 2-rank>=2 at every prime and found nothing; a new table-based
   3-only sieve (`scan230b.c`, mode `sieve5t`: precomputed kill tables over
   F_p^2 for p in {5,...,59}, ~ns/pair) instantly found 76 survivors at
   height <= 60 (den <= 8), 1143 at height <= 600.

4. **JACKPOT-ADJACENT: an infinite, provably-[>=15], generically simple family**
   on the q2=-1/4 boundary of Elkies' A_1(5): the line `q1 = 1, q0 = c^2`.
   With `Q = -x^2/4 + x + c^2`, `Q' = -(x-2c)(x+2c)/4`, the curve
   `y^2 = F(x) = (Q'-xQ)^2 + 4Q^2Q'` has the four rational points
   `(±2c, ±4c^2)` and TWO torsion classes supported on them:
   - 5-torsion: Mumford `(x^2-4c^2, 4c^2)` (order 5 verified for many c),
   - 3-torsion: Mumford `(x^2-4c^2, 2cx)`, proven for ALL c by the exact
     contact identity over Q(c) (Groebner, `family15.m`):
     `h3^2 - F = (1/64c^2)*(x^2-4c^2)^3`, `h3 = x^3/(8c) - (c/2)x^2 + (3c/2)x + 2c^3`.
   So torsion contains Z/15 for every c with F squarefree; exact torsion [15]
   verified at c = 1, 2, 3, 7, 1/2 and all square q0 <= 60 in the sieve box.
   Simplicity certified (3 primes each, deg-4-irreducible L-poly AND
   12th-power transform irreducible) at q0 = 1, 4, 9, 2, 49/4.
   For non-square q0 the same identity holds over Q(sqrt(q0)) with the class
   anti-rational (sigma D = -D): explains the whole line q1=1 passing 3|#J at
   every prime while having exact torsion [5] (systematic sieve false-positive
   analogous to the known q1=-(3/4)q0 rank-2 line).
   Bonus sporadic simple [15]s off the line: (q0,q1) = (-36, 11/2), (-2, 16/7),
   and (2, 1) (q0=2 non-square: a DIFFERENT 3-torsion class than the contact-k
   above; all three certified simple).

5. **Five certified-simple exact [30] curves** — the odd-degree extra-2 cut.
   The condition "F(x; c^2, 1) has a rational root r" defines a plane curve
   that is irreducible of **GENUS 0** (`family15.m` part iv), with rational
   points found at height <= 1000:
   (r, c) = (45/4, 27/8), (245/18, 49/9), (845/72, 169/36), (637/36, 245/72),
   (245/36, 49/24) [and sign mirrors]. ALL FIVE give factorization type (1,4),
   exact torsion **[30]**, and 3-prime simplicity certificates (`family30.m`):
   ```
   c=27/8   r=45/4    torsion [30] certs 19,23,31
   c=49/9   r=245/18  torsion [30] certs 23,37,41
   c=169/36 r=845/72  torsion [30] certs 11,17,19
   c=245/72 r=637/36  torsion [30] certs 11,19,59
   c=49/24  r=245/36  torsion [30] certs 17,19,29
   ```
   Genus 0 + infinitely many points => an infinite [30]-bearing family
   (second known infinite [30] construction, first inside A_1(5)).
   **(2,30) is now exactly one factorization condition away**: on the
   parametrized genus-0 curve, the quartic cofactor G4 = F/(x-r) must factor
   (rational root => type (1,1,3) => 2-rank 2 => [2,30]).

## Structural facts for the record

* The fiber of the contact system over a generic rational (k1,k0) has its
  l-eliminant factoring as degrees [1, 56, 112] (5 random fibers, `probe_fibers.gp`),
  the linear factor being the constant parasite l = -1/108 (never a square).
  So no low-degree branch over the (k1,k0)-plane; the rational points come from
  the special line s=0 (k = x^2 - 4c^2) found by the sieve.
* F(x; c^2, 1) = -(1/8)(x^5 - (2c^2 + 9/2)x^4 + 8c^2x^3 + (16c^4-12c^2)x^2
  - 48c^4 x - 32c^6 - 8c^4)... (see `family15.gp` output for exact form);
  irreducible over Q(c), lc = -1/8 < 0.

6. **Clean parametrization of the [30] family.** Projecting the genus-0 root
   curve from its multiplicity-4 singular point (0,0) via c = t*r gives residual
   conic `A r^2 + B r + C` with `A = -4t^2(4t^2-1)^2`, `B = -96t^4+16t^2+2`,
   `C = -(4t^2+3)^2`, discriminant `Delta = -4(4t^2-1)^5`; rationality <=> the
   conic `w^2 = 1-4t^2`, parametrized by `t = g/(1+g^2)`. Result
   (`family2230.m`):
   ```text
   c(g) = g(g^2+3)^2 / (2(g^2-1)^2),   r(g) = (g^6+7g^4+15g^2+9)/(2(g^2-1)^2)
   ```
   (g=3 -> (c,r) = (27/8, 45/4)). EVERY rational g (12 tested, g != 0,+-1)
   gives exact torsion [30] with a 3-prime simplicity certificate
   (`data/claude_prod_08_230_curves.txt`: 12 certified [30] models).
   This is an infinite [30]-bearing family with certified-simple members at
   will — the second known [30] construction, first inside A_1(5).

7. **The (2,30) gate = two explicit higher-genus curves (both computed):**
   - (1,1,3) route (second rational root of the quintic): the curve
     {G4(rho; m) = 0} over the parametrized family is irreducible of
     **genus 5** (`family30b.m`).
   - (1,2,2) route (rational quadratic factor of the quartic cofactor G4):
     S122 = {(a,b,g) : x^2+ax+b | G4(x;g)} is 1-dimensional with a single
     irreducible component (`family2230.m`); genus run: `hunt2230b.m`.
   - Point hunts so far: 138 rational points on the root curve to height
     2*10^6 — all 69 distinct-c members give factype (1,4) (no extra factor,
     `hunt2230.m`); duplicate-c^2 scan + same-c correspondence curves
     (`hunt2230b.m`, `same_c.m`) pending below.

8. **Generators confirmed numerically** (`confirm_gens.m`): on the integral
   model `y^2 = den^2*F(x; c^2, 1)`, `D3 = (x^2-4c^2, 2c*den*x)` has order 3,
   `D5 = (x^2-4c^2, 4c^2*den)` has order 5, `D3+D5` order 15 — verified at
   c = 1, 2, 3, 1/2, 7, 27/8 (the last being a [30] member).

9. **Sporadic [15] 3-torsion data** (`sporadic_gens.m`): the k-quadratics of
   the sporadics' 3-parts are x^2-8x+24 (q0=2,q1=1), x^2-(39/2)x+108
   (-36,11/2), x^2-(60/7)x+8 (-2,16/7) — unrelated to the family's
   k = x^2-4q0 and to each other: genuinely isolated points of the {T3,T5}
   locus (no second line/family visible).

10. **Same-c correspondence closed at small height**: both branches
   {c(g')=c(g)}/(g'-g) and {c(g')=-c(g)}/(g'+g) are irreducible plane curves
   of degree 8 and **genus 5**; PointSearch to 10^5 finds only trivial points
   (`same_c.m`). No duplicate c^2 among all 138 root-curve points to height
   2*10^6 (`hunt2230b.m`). The swap involution (u,v)->(v,u) preserves both
   branches, so each has a symmetric quotient of genus <= 3 in (e,p)=(u+v,uv)
   coordinates — the designated Chabauty target (`sym_quot.m`).

11. **Symmetric quotients** (`sym_quot.m`): minus branch {c(g')=c(g)} quotient
    is the plane quartic `e^4 - 2e^2p^2 - 12e^2p + 6e^2 + p^4 - 4p^3 - 2p^2 + 12p + 9 = 0`
    ((e,p)=(g+g', gg')) of **genus 1** (two nodes at infinity [1:+-1:0]);
    plus branch quotient has genus 3. Rank analysis of the genus-1 quotient in
    `equot.m` — if it has positive rank, its point orbit is an infinite source
    of [2,30] candidates (lift condition: e^2-4p a nonzero square).

12. **Deep sieve H=2000, den<=10** (`sieve5t_h2000.out`): 10794 survivors,
    8170 off the q1=1 line — mostly mod-p accidents at this box size
    (density scaling); needs an extended-prime refilter (p up to 127 fits the
    chitab) before exact checks. Not yet processed.

13. **THEOREM (unconditional): the (1,1,3) route is CLOSED.** The genus-1
    minus-branch quotient is birational to the elliptic curve
    `y^2 = x^3 - x^2 + x` (**conductor 24, rank 0** by RankBounds, torsion
    Z/4, `equot.m`), whose 4 rational points correspond exactly to the 4
    trivial affine points (e,p) = (-2,1), (0,-1), (2,1), (0,3). Since the
    g-parametrization of the root curve is birational (g <-> 1/g resolves the
    two conic roots r+-, exceptional set only c=0/infinity), a second rational
    root of any member's quintic would force a nontrivial rational point on
    the minus branch {c(g')=c(g)} — which does not exist. Hence NO member of
    the g-family has quintic factorization type (1,1,3): within this family,
    (2,30) can only arise via the (1,2,2) route (C122). (The plus branch
    {c(g')=-c(g)} is subsumed by the same argument: a plus-branch pair would
    also produce a minus-branch pair through (rho, -c) -> (rho, c).)

14. **(1,2,2) route: genus 5 as well, empirically closed to large height.**
    Eliminating b by hand (b = (p1-aS)/(p3-2a), S = p2-a(p3-a), then bd = p0)
    gives a plane model PC(a,g) of the quadratic-divisor locus C122, bidegree
    (6,30), EVEN in g, descending to PC2(a,h), h = g^2, bidegree (6,15) with
    a^6-coefficient -(h-1)^14 (`c122plane.gp`, `pc2_for_magma.txt`). PC2 is
    irreducible of **genus 5** (`pc2genus.m`). Direct hunt: the quartic
    cofactor G4(x;g) was factored over Q for **175,280 rational g**
    (|num|<=1200, den<=120): **zero splits** (`g4hunt.gp`, `g4hunt2.gp`).
    So within the g-family, (2,30) requires a rational point on one of two
    genus-5 curves; none exists at small height, and the (1,1,3) one is fully
    closed (item 13).

## Verdict

**(2,30) does not fall this session, but the target moved decisively closer
and produced two new infinite families with certified-simple members:**

* Exact **[15]** (geometrically simple, certified): infinite family
  q0 = c^2, q1 = 1, q2 = -1/4 on Elkies' A_1(5); explicit torsion generators
  (x^2-4c^2, 2cx) [order 3] and (x^2-4c^2, 4c^2) [order 5]; 3 isolated extra
  members. Possibly the first recorded exact-[15] simple examples (literature
  check pending: Nicholls' Z/30 contains Z/15 as a subgroup but exact [15]
  simple examples are not in the repo's records).
* Exact **[30]** (geometrically simple, certified, 12 members in data):
  the g-subfamily c = g(g^2+3)^2/(2(g^2-1)^2) — every tested g gives [30].
* The only two ways this family reaches (2,30) are rational points on two
  explicit genus-5 curves; the (1,1,3) one is unconditionally CLOSED
  (rank-0 conductor-24 elliptic quotient), the (1,2,2) one has no points to
  height ~10^3 in g. A first-ever (2,30) therefore most likely needs either
  a rational point on PC2 (Chabauty target) or a different order-30
  construction (interior deg-6 chart of A_1(5), Route 2 A(6)+5, or new
  sporadics from the H=2000 sieve refilter).

## Resume state / next steps

1. **PC2 Chabauty (the remaining gate):** PC2(a,h) in
   `prod230/pc2_for_magma.txt` (genus 5, irreducible). Look for extra
   involutions/quotient maps (the a <-> (complementary quadratic) swap acts:
   a' = p3(h) - a — check if it preserves PC2 and quotient to lower genus);
   then RankBounds/Chabauty on a quotient. Any rational point with h = g^2,
   g in Q, h != 0,1 gives a (2,30) CANDIDATE — jackpot protocol.
2. **Sieve refilter:** `prod230/sieve5t_h2000.out` (10794 survivors, 8170
   off-line, unprocessed). Extend `scan230b.c` tprimes to p <= 127 and
   refilter before exact checks (expect O(10-100) true candidates; each new
   off-line [15] is a new isolated point of the {T3,T5} locus and a
   potential seed).
3. **Interior chart:** repeat the whole construction on the deg-6 chart
   (q2 free): the analogous family ansatz is Q' with rational root pair
   symmetric about the chart — the q1=1 trick came from Q' = Q - x becoming
   even; look for the analogous locus with TWO stable quadratic factors of
   the sextic F for the extra (2,2).
4. **Literature check:** is exact [15] simple new? (LMFDB g2c torsion query +
   Nicholls thesis.)
5. Unfinished/killed jobs: taskC_setup.m and taskC_elim.m (full-ideal
   saturations, superseded by the direct construction); hunt2230b.m part (ii)
   3-space C122 genus (superseded by pc2genus.m).

## Files

```text
notes/claude_prod_08_230.md                     this file
data/claude_prod_08_230_curves.txt              22 certified curves: 12x[30], 3 sporadic [15], 5 family [15] samples (+ generators/parametrization in header)
code/claude_prod08_scan230b.c                   table-based pure-3 sieve (sieve5t mode)
code/claude_prod08_family15.m                   family verification: orders, symbolic contact Groebner, simplicity certs, root-curve genus
code/claude_prod08_family2230.m                 g-parametrization + C122 construction
code/claude_prod08_same_c.m                     same-c correspondence curves (genus 5 + 5)
code/claude_prod08_confirm_gens.m               Mumford generator confirmation
scratchpad prod230/: taskA_510.{m,out} contact_scan.c contact_p*.out scan230b.c
  sieve5t_h60/600/2000.out check3*.out family15.{m,gp,out} family30{,b}.{m,out}
  hunt2230{,b}.{m,out} same_c.{m,out} sym_quot.{m,out} equot.{m,out}
  c122plane.gp pc2_for_magma.txt pc2genus.{m,out} g4hunt{,2}.{gp,out}
  param_clean.gp probe_fibers.gp elim_res.gp datawrite.m sporadic_gens.m
```

## Caveats

* Simplicity certs are per-member (3 good primes each, deg-4-irreducible
  L-poly with irreducible 12th-power transform); "infinite family of simple
  members" means: certified members at will, not a uniform proof of generic
  simplicity (End jumps are on thin subsets anyway).
* Exact torsion [15]/[30] verified member-by-member (Magma TorsionSubgroup on
  integral models); the family-level statements proved symbolically are the
  CONTAINMENTS Z/15 (resp. Z/30) <= J(Q)_tors.
* The (1,1,3) closure argument uses birationality of the g-parametrization
  with exceptional set {c=0, r=+-2c, infinity}; r=+-2c never occurs
  (F(+-2c)=16c^4 != 0), c=0 is degenerate. Double-checked numerically:
  no duplicate c^2 among 138 root-curve points to height 2e6.
* The old level-B sieve requirement (2-rank>=2 at all p) was the reason the
  entire [15] stratum was invisible in the test session — lesson recorded.

