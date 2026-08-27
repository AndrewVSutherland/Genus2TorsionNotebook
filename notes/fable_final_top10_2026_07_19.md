# Fable final ranked top-10 and strategies — 2026-07-19

Baseline: the 2026-07-18 frontier note, PR #2 (record + top-10 campaign),
PR #3/#4 (second [2,2,2,12] curve + T5 structure theorems; [4,16] tier-2
cold to w-height 1e6), the [21]/[49]/[63] odd-torsion investigations (no
new realizations; Z/63 exists only split via HLP; [3,21] locally killed to
h=5000), the overnight [2,2,2,24] component closures (S1/W1/V/U2), and
this session's curve-2 audit (below).

## New today: curve-2 Richelot/divisibility audit (fable_curve2_audit.m)

Curve #2 ((s,m,n)=(2208,-8303,-7200), A=[1,1,1,2,2],
B=[25648128,-36568896,-52466496,-59781600,23309856]): exact [2,2,2,12]
reconfirmed (order 96).  Depth-3 rational Richelot component: 18 vertices,
census exactly 1 x [2,2,2,12] + 12 x [2,12] + 5 x [2,2,6] — IDENTICAL to
the record's component.  Divisibility: exactly one of the 15 nonzero
2-classes halves (the existing 4-direction), none of the 16 order-4
classes halves (no 8-chain), none of the 32 order-12 generators halves.
Conjecture (component rigidity): every T5-pencil [2,2,2,12] fiber has
this same 1-12-5 census and divisibility pattern.  Consequence: no
fiberwise route to [2,2,2,24]/[2,2,4,12]/[4,12]; those targets are global.

## Ranked top-10 open targets with strategies

1. [8,8] (64).  Genus gate PASSED (stage-1 cover fibers genus 0/1); local
   picture healthy everywhere probed.  Strategy: write the remaining
   conic-type lift layer symbolically over the genus-0/1 stage-1 fibers;
   compute the conic's obstruction primes (Hasse); parametrize rational
   points; exact-test second-halving lifts with an [8,8]-containment gate;
   certify.  The one remaining layer is explicit — most actionable target.
2. [2,24] (48).  Both loci dense, intersection empty at 6,426 exact
   checks, no local obstruction: either thin or globally incompatible.
   Strategy: prove-or-break the suspected 2-descent incompatibility on the
   2-parameter certified-simple [2,12] family: write the +3-contact cover,
   decompose (saturate split/singular components), compute genus/rank of
   the residual; a rational point realizes [2,24], an ELS/Sha certificate
   explains the emptiness.  Richelot moves cannot help (fixed seeds are
   kernel-poor; codomain redistributes down).
3. [3,12] (36).  Smallest open order.  The M-square-class dichotomy {1,-3}
   is conjecturally the localized obstruction; nonsplit carrier points
   exist with (3,12) rational over Q(zeta_3).  Strategy: prove the
   dichotomy symbolically; construct the M=1 stratum cover; boundary-first
   CRT with prescribed bad reduction at {5,7} (both order-12 charts are
   chart-level dead at 5 and 7, alive at 11).
4. [2,2,24] (96).  Ten split HLP anchors.  Strategy: two-parameter section
   plane through two anchors, transverse to the visible Humbert branches;
   eliminate/saturate; keep only low-degree residual components; funnel +
   certify.  Avoid the M_1(8,2,2)+3 shortcut (chart-level dead at 7,11,13).
5. [6,12] (72).  Strategy: finish the certification chain (exact quotient
   map to E8, bigonal transport), then the precision-safe p=37 Coleman
   driver plus the Abel-Prym Mordell-Weil sieve.  Prym rank 1 < 2: this
   lane will either produce the curve or close, and both are progress.
6. [2,2,16] (64).  The p=7 wall is chart-fundamental, p=11 is
   cover-specific, 13/17 alive.  Strategy: build the 7-adic boundary chart
   of the halving cover, CRT-lift {7-boundary, 11-cover} conditions,
   bounded exact tests on the lifted strata.
7. [4,12] (48).  Strategy: (a) mod-7 boundary blowup of the degree-8
   reduced s=m^2 equation on full M(2,12); (b) new: write the "[4,12] on a
   Richelot codomain" condition over the T5 pencil — a [2,12]-type leaf
   whose dual-kernel 2 lands inside the 12-direction — and search pencil
   members, starting at the near-miss classes 3/6/10/-1 (14,903-point
   table in data/claude_sib_t5_nearmisses.txt).
8. [2,6,6] (72).  Demoted: split confinement is mechanistic (two infinite
   split families explain the census; the mu_3 part is never rational in
   the BFS escape).  Strategy: only a genuinely new chart can work —
   simultaneous double-3-contact compatibility variety plus one rational
   2-torsion; decompose first, search only low-genus pieces; CRT at {5,7}.
9. [2,2,2,24] (192).  Both known [2,2,2,12] fibers have non-divisible
   12-generators (record audit + today's curve-2 audit); S1/W1/V/U2
   components of the global cover are closed.  Strategy: pencil-level
   divisibility: impose the halving square-condition along the T5 pencil
   and enumerate MW lattices member-by-member (the member-closure
   machinery of data/claude_sib_member_closures.txt is exactly this).
10. [2,2,4,12] (192).  Both known fibers have exactly one divisible
    2-class.  Strategy: write the 15 Stoll square-class halving conditions
    symbolically over the A(2,2,2,12) chart (the prod-07 delta machinery
    is reusable verbatim), quotient by level-2 symmetry, and characterize
    the "second divisible class" cover; kernel-first Richelot lowers the
    degree.

Alternates (in order): [2,2,4,8] (boxed on three sides — only the HLP
section-surface route remains), [4,16] (tier-2 cold to w-height 1e6),
Z/48, [60], Z/35, Z/5 x Z/5, Z/63 (split exists), [3,21], (7,7).

## Standing rules worth keeping

- Fan every new simple curve with torsion order >= 48 through a depth-1
  Richelot audit the moment it is found (seconds; jumps are real but
  directional — you must land on the special fiber).
- 3-part-mixed targets: CRT/bad-reduction-first at the walled small
  primes; never plain affine scans.
- The T5-component-rigidity conjecture above is cheaply falsifiable on
  every future [2,2,2,12] hit — test it each time.
