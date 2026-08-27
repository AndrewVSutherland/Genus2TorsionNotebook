# Split-Jacobian torsion census campaign (started 2026-08-13)

Goal (Filip): companion table to the geometrically-simple census in the
paper, for GEOMETRICALLY SPLIT genus-2 Jacobians over Q.  Sources: LMFDB,
Howe--Leprevost--Poonen (HLP 2000, primary literature), Frengley's
N-congruences repo (N > 3 gluing), and new systematic searches.

## Source 1: LMFDB (exact torsion, conductor <= 10^6)
Query: db.g2c_curves, is_simple_geom = False: 2926 curves, 31 distinct
torsion groups, smallest-conductor witness recorded per group (script
code/split_census_aggregate.py will regenerate).  Notable: [19] at
conductor 169 (torsion order impossible for the elliptic factors --
purely a gluing phenomenon).  NO cyclic [16] in LMFDB split data.

## Source 2: HLP 2000 (arXiv math/9809210; Forum Math 12 (2000))
Table 1 (genus 2), 20 groups, mostly CONTAINMENT with families:
  [20],[21],[3,9],[30] over P^2; [35],[45],[60],[2,2,24] over pos-rank
  elliptic curve; [40],[5,10],[6,12],[2,6,6],[2,2,4,8] over pos-rank
  elliptic surface; [6,6],[3,12],[2,24],[8,8],[2,4,8] over P^2;
  [7,7],[63] single curves (P^0) with EXACT torsion proven in-paper.
Both exact singletons INDEPENDENTLY VERIFIED here
(code/split_hlp_verify.m, results/split_hlp_verify.log):
  y^2 = x^6+3025x^4+3232987x^2+869675859           -> [7,7]
  y^2 = 897x^6-197570x^4+79136353x^2-146398496     -> [63]
KEY CROSS-CHECK for the simple-census paper: HLP's order-96 group is
[2,2,24] and order-128 is [2,2,4,8] -- NOT [2,2,2,12]; the paper's
"first [2,2,2,12] simple or split" claim survives.

## Source 3: Frengley N-congruences repo (github/SamFrengley, v0.1.0)
Z(N,r) Hilbert-modular-surface models for N <= 19; genus-2-curves/N.m
holds example glued curves (optimal degree-N covers of (N,-1)-congruent
pairs) for N = 5,7,15,16,17.  Their torsion computed here
(code/split_frengley_torsion.m): [2],[3],[],[],[]  -- small, as expected
(examples were not torsion-optimized).  A torsion-aware use of the
Z(N,-1) models (rational points whose BOTH projections have rational
torsion) is the natural follow-up campaign.

## Source 4: new bielliptic sweep (2-gluing leg)
code/split_bielliptic_sweep.m: y^2 = a x^6 + b x^4 + c x^2 + d over
|coeffs| <= 12 (3 shards, gcd(#J(F_p)) >= 12 prefilter, exact torsion on
survivors).  Split-by-construction (x -> -x involution).  Hits so far
include cyclic [16] (NEW: absent from LMFDB split data and from HLP),
[2,2,8] (order 32), [24], [21], [15], [4,4], [2,8], [2,6], [12] with
tiny witnesses, e.g. [16]: y^2 = -8x^6-7x^4+6x^2+9.

## TODO
- finish shards; aggregate into paper-style longtable (draft generator
  code/split_census_aggregate.py)
- extend sweep height; consider odd-degree double covers (y^2 = x*g(x^2))
- 3-gluing systematic leg (HLP Lemma 11 parametrization)
- torsion-aware Z(N,-1) search for N = 5,7 (Frengley models)
- merge literature beyond HLP if found (user said HLP primary)

## 2026-08-14: HLP-Prop-4 gluing sweep — NEW GROUP [5,5]; [35] upgraded to exact

code/split_hlp_glue_sweep.m implements HLP Proposition 4 (explicit glued
sextic for E,F with trivial rational 2-torsion and isomorphic 2-division
fields) over Kubert-parametrized X_1(5)/X_1(7)/X_1(9) Tate curves,
disc-mod-squares + cubic-field-isomorphism funnel, all 6 root matchings,
h-rationality test, exact torsion.  VALIDATION: rediscovers HLP's [7,7]
at (t,u)=(7,-14/13) and [63] at (-16/3,4), both exact. ✓✓

NEW: exact [5,5] (in no source: not LMFDB-split, not HLP, not bielliptic
sweep) -- FOUR curves (code/split_glue_extract.m):
  (t,u)=(-1,11):  y^2+(x^3+x^2+x+1)y = x^5+2x^4+5x^3+2x^2+4x-2
                  Jac ~ glue of 11a1 = X_0(11) and 11a3 = X_1(11) (!)
  (t,u)=(-2,8):   y^2+(x^2+x)y = x^6-3x^5+7x^4-10x^3+7x^2-3x+1   (50b1,50b2)
  (t,u)=(-7,7/5): 175a2 x 665d1;  (10,11/2): 110a1 x 550k2.
UPGRADED: [35] from HLP-containment-family to FIVE exact witnesses, e.g.
  (13/2,-1): y^2+(x^2+x)y = 10x^6+30x^5+116x^4+182x^3+166x^2+80x+160
  (26b1 x 286d1).
Second-63 hunt: in the validation box exactly ONE [63] curve (HLP's) --
their uniqueness suspicion survives; wider 9x9/7x9 box running
(results/split_glue_99wide.log), bielliptic H=25 sweep running.

## 2026-08-14 (later): corrections + 3-gluing sweep launched

CORRECTION to the [5,5] claim: our 11a1x11a3 glued curve IS X_0(22)
(G2Invariants match HLP Cor 5's y^2 = -2x^6-10x^4+26x^2+242, which HLP
identify as X_0(22)); [5,5] = J_0(22)(Q)_tors is CLASSICAL.  Census row
should read: group known via X_0(22) (Ogg/Ligozat era; HLP Cor 5),
absent from LMFDB-g2c and HLP Table 1; three ADDITIONAL witnesses new
here (175a2x665d1, 110a1x550k2 non-isogenous pairs; 50b1x50b2).
NEGATIVE (clean): Cor-7 self-gluing route to [9,9]: NO square-disc
2-division cubics on X_1(5)/X_1(7)/X_1(9) families in |num|<=24,den<=16
(results/split_selfglue.log) -- likely modular obstruction (mod-2 image
in A_3 incompatible with the torsion); [9,9] via this route is cold.

3-GLUING (user request) via BHLS arXiv:1403.6911 Appendix: the complete
rational family C_{a,b,c,d,t}: t y^2 = (x^3+3ax+2b)(2dx^3+3cx^2+1),
12ac+16bd=1, with explicit elliptic quotients f1, f2.  Sweep
code/split_glue3_sweep.m (3 shards): torsion of both quotients as fast
prefilter (product >= 25), exact genus-2 torsion on survivors.  Targets:
[56](7x8), [70](7x10), [2,40](8x10), [10,10](10x10), plus whatever the
family finds that pair-thinking missed.

## 3-gluing round-2 verdict (2026-08-14/15)

Blind BHLS-family grid (Hn=12,Hd=6, 10 twists, ~2.2M members): glued
torsion reached [20], [2,10], [3,6] (order-18 with a genuine coset gain
over the naive quotient), many [6]/[2,6] -- NO new groups.  The 7-, 8-,
10-torsion quotient cells are codim-1 loci the grid cannot land on.
ROUND-3 DESIGN (recorded, not yet run): solve the slice
{(a,b,c,d) : E_1 = (t y^2 = f_1) has 7-torsion} by matching j(f_1-curve)
to the X_1(7) j-map -- a surface in (a,b,c;s); search THAT with the
7-divisibility mod-p sieve; likewise 8/10-slices.  Targets remain [56],
[70], [2,40], [10,10].
H=25 bielliptic shards died at 3GB (G2-invariant dedup set); relaunched
with Hash dedup + 8GB (results/split_bielliptic_H25b_s*.log); partial
H25 logs already show [21],[2,8],[3,6] at height 25.

## 2026-08-15: DB-driven reducible-2-gluing sweep — THREE new groups + five exactness upgrades

code/split_glue2_db_sweep.m over the LMFDB pool (739 curves, torsion
[8],[10],[12],[2,4],[2,6],[2,8]), 2-division-bucket matching, cached
point-count gcd prefilter, Prop-4 gluing, exact torsion.  33k+ glued
curves, 15 distinct exact groups (mid-flight aggregate).

NEW SPLIT GROUPS (census 49 -> 52):
  [2,2,12]  (5 pairs; best 15.a6 x 90.c3):
    y^2+(x^2+x)y = x^6+3x^5-11x^4-27x^3+44x^2+58x-69
  [2,2,2,8] (2 pairs; 21.a5 x 210.e6):
    y^2+(x^2+x)y = x^6-3x^5-9x^4+22x^3+30x^2-42x-45
  [2,2,4,4] (210.c5 x 2310.o4):
    y^2+(x^2+x)y = 60x^5+1000x^4-671x^3-5657x^2+867x+4913
  (NB: [2,2,2,8] and [2,2,4,4] are also new SIMPLE realizations in the
  paper -- both groups now exist on both sides of the divide.)

EXACTNESS UPGRADES of HLP containment rows (explicit tiny witnesses):
  [40]:   y^2+(x^3+x)y = 2x^4+12x^2+66            (48.a6 x 66.c4)
  [60]:   y^2+(x^2+x)y = x^6-7x^5+25x^4-24x^3+25x^2-7x+1  (150.c3 x 90.c7)
  [5,10]: y^2+(x^2+x)y = 9x^6+3x^5+15x^4+32x^3+11x^2-135x-135
  [2,24]: y^2+(x^3+x)y = -x^4-6x^2+15             (48.a3 x 30.a6)
  [6,12]: sweep-log witness 8190.bx1 x 155610.fa1 (extraction retry TODO)

Structure note: different matchings of the SAME pair give different
exact groups (e.g. the [2,2,12]-pair also produces [2,12]-only curves)
-- the choice of 2-torsion matching genuinely matters.
No order > 72 so far; no split [2,2,2,12] (the paper's claim stays safe);
(2,6)x(2,6), (2,6)x(2,8), (2,8)x(2,8) blocks still streaming.

## 2026-08-15 (CORRECTION): paper/split_torsion_table.tex supersedes the sweep "news"

The paper repo contains split_torsion_table.tex (75 groups,
generated from the EXTENDED Booker-Sutherland snapshot 2026-08-12 plus an
earlier project stage's gluing constructions; companion of
torsion_realizations.tex).  Verification here
(code/split_table_verify.m): ALL NINE displayed non-database equations
([45],[63],[70],[6,12],[7,7],[8,8],[2,2,24]=96,[2,2,4,4],[2,2,4,8]=128)
have EXACTLY the stated torsion. ✓

Consequences for this campaign's claims (ERRATUM discipline):
- Every group found by the sweeps (31 distinct) is ALREADY in that
  table; my "new split groups" announcements ([16],[2,2,12],[2,2,2,8],
  [2,2,4,4],[2,2,16], the [35]/[40]/[60]/[5,10]/[2,24]/[6,12] exactness
  upgrades) were new only vs the PRODUCTION-LMFDB+HLP baseline, not vs
  the extended DB + prior-stage table.  Census bookkeeping: use the
  table's 75 as the base count.
- Independent-confirmation value: several glued curves REPRODUCE the
  extended-DB minimal witnesses exactly (my [40]-glue = 3168.i.1, [60] =
  13500.r.4, [2,24] = 1440.c.2, [35] = 7436.a.2) -- the Prop-4 pipeline
  and the DB validate each other.
- SURVIVING new contributions: (a) a much smaller [2,2,4,4] witness
  (210.c5 x 2310.o4: y^2+(x^2+x)y = 60x^5+1000x^4-671x^3-5657x^2+867x
  +4913) vs the table's 10-digit-coefficient model -- suggest swapping
  into the table; (b) additional small witnesses for several rows;
  (c) the negatives: [9,9] resists three routes; [63]-uniqueness in
  large boxes; A3-self-gluing empty on X_1(5)/(7)/(9); (d) the BHLS
  (3,3)-family sweep machinery + round-3 targeted design; (e) [5,5] =
  X_0(22) contextual identification.

## 2026-08-15 (final): DB-gluing sweep curtailed at 464k gluings

464,035 glued curves, 17 distinct exact groups, ALL inside the 75-group
table -- the sweep became a witness mill; shards stopped mid-(2,4)-block.
Final standing of the campaign: strong mutual validation of the table
(several table witnesses reproduced exactly by independent gluing),
better [2,2,4,4] witness (adopted into the paper's split table),
citable negatives ([9,9] three-route wall, [63]-uniqueness boxes,
A3-selfglue empty, bielliptic-H25 completeness), and the reusable
machinery (Prop-4 DB gluing, BHLS (3,3)-family sweep).
