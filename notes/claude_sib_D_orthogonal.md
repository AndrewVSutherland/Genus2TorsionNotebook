# Lane D: orthogonal charts + the (2,2,2,12) hit curve's neighborhood

Date: 2026-07-18.  Sibling-hunt session following the first (2,2,2,12)
realization (notes/claude_prod_02_22212.md).  Scratchpad: `sib_D/`
(surfS.gp, batchscan.sh, memberjobs.sh, neighborhood.m, twists.m,
g2dedupe.m, runchunks_ext.sh); reused from `prod_02_22212/`: tor22212 +
t127.bin + postfilter, t5sweep, hitverify.m.  Data banked:
data/claude_sib_D_neighborhood.txt, data/claude_sib_D_t5_surface.txt.

## 0. Headlines

1. **No new (2,2,2,12) curve found** — but TWO NEW REPRESENTATIONS of the
   known hit curve were found and fully verified at u-heights (7105,
   13872) far beyond all previous sweeps, via a NEW search method (the
   C_rho' surface scan, below) that is ~1000x cheaper per unit of
   u-height than box sweeps.  The known curve now has FOUR (u,rho')
   representations on the T5 pencil:
   (-97/48,-49/240), (133/145,289/240)  [discovery session]
   (3637/7105,-49/240) -> (s,m,n)=(-527365,174576,-1030225)  [NEW]
   (6767/13872,289/240) -> (s,m,n)=(81204,-981215,166464)    [NEW]
   Both new ones: TorsionSubgroup exactly [2,2,2,12], simplicity certs
   p=37,73,113, G2Invariants identical to the hit (g2dedupe.log) —
   i.e. each hit member carries a second u-point (a second involution).
2. **The C_rho' criterion** (exact, one line from the proven product
   identity): a T5 hit on member rho'=rn/rd requires a rational point on
   the genus-1 quartic C_rho': y^2 = (q rn - rd)(q rn - (2u-1) rd),
   q=4u^2-6u+3, and ON C_rho' the identity makes X4 free: X1,X2,X3
   squares => full hit.  PARI hyperellratpoints enumerates C-points to
   u-height 3000 in ~1ms/member => member-boxes of 10^5-10^6 members are
   feasible at u-heights 20-1000x beyond the discovery sweeps.
3. **Isogeny neighborhood of the hit**: all 15 Richelot codomains and
   all 17 TwoPowerIsogenies codomains are Jacobians (no products — as
   forced by simplicity) with torsion order EXACTLY 24: [2,12] x12,
   [2,2,6] x5.  Torsion degrades 96 -> 24 uniformly (the extra Z/2 x Z/2
   and the 4|12 both drop); no new torsion group appears in the 2-power
   neighborhood.  Bad primes of the minimal model:
   {2,3,5,7,11,13,17,19,23,29,53,109,193,607,1013,1423} — confirms the
   forced-bad-reduction theorem (bad at 3,5,7,11,13).
4. **Quadratic twists**: all 37 squarefree |D|<=30 twists have torsion
   EXACTLY [2,2,2,2] (full 2-torsion is twist-invariant; 4- and 3-parts
   die in every twist tested).
5. **Failing-class constancy law** (new empirical law, 4 members x 2
   points): the failing square class at a 3-of-4 near-miss equals
   class(F(u)) (F = the C_rho' quartic) and is CONSTANT along each
   member's near-miss points: -1/143 -> 3; -25/551 -> 6; -169/1431 -> 10;
   841/697 -> -1; hit members -> 1.  Two NEW near-miss points found (the
   "pair law": every known member has exactly 2 points): u=89/169 on
   -169/1431 and u=633/841 on 841/697 (both X4-fail, same class as
   partner; note 169=13^2, 841=29^2 — near-miss u tend to have square
   denominators).
6. **MW ranks of C_rho'**: all six known members have ellrank = 3
   (841/697: bounds [2,3]), torsion Z/2 — NO rank-0 member kill exists;
   the true obstruction lives on the X1X2X3-cover of C, not on C.
6b. **Hit-representation structure** (exact, gp-checked; details in
   data/claude_sib_D_t5_surface.txt): the 4 representations pair under
   the cross-member involution with q(u)q(u') an exact square per pair;
   at every representation class(A)=class(B)=+-15 mod squares
   (A=q rho'-1, B=q rho'-2u+1), with A,B numerators squares of BAD
   primes of the curve (1013^2, (19*23)^2, (17*1013)^2, (17*19*23)^2)
   and 2u-1 numerator +-11^2/+-13^2 (the forced-bad primes).  The
   second representation on each member sits at canonical height
   31.4433... (identical on both members) relative to the first.
6c. **LMFDB context**: the largest Jacobian torsion order among ALL
   g2c_curves in the LMFDB is 39 (one curve, cond 1116); order 96 is
   far beyond anything catalogued.
7. **Enumerator extension (quintic chart y^2=x(x+a^2)...(x+d^2))**: the
   unconditional-kill sieve (-P127 tables + 96|#J postfilter p<=397) was
   extended beyond d=2000 in chunks of 100 (see Section 3 tallies and
   resume state).

## 1. The C_rho' surface scan (new tool)

`sib_D/surfS.gp`: for member rho'=rn/rd let
P(u) = (q(u) rn - rd)(q(u) rn - (2u-1) rd)  (integer quartic, q=4u^2-6u+3).
V'1V'2V'3V'4 == P/rd^2 mod squares (proven identity), so a hit needs
y^2=P(u) solvable (z != 0) — and then X4 is automatically square iff
X1X2X3 is.  qscan(rn,rd,H): hyperellratpoints(P,H), drop u in
{0,1,1/2} and the X3=0 or X4=0 degenerate sections, exact-test X1,X2,X3
(integer-cleared forms from t5sweep.c).  Validated: re-finds all 4
representations of the hit and nothing else on the six known members.

Scans completed (u-height H, all counts EXCLUDE degenerate sections):
* six known members at H=1e5 and the two hit members at H=1e6:
  only the hit's 4 representations; near-miss members carry NO hit
  to u-height 1e6/1e5 despite 44-62 C-points each.
* Pass A: ALL coprime members |rn|<=300, rd<=300 at H=3000:
  109587 members, full hits = only the hit's 2 low-height reps.
* Pass B: rn=+-a^2 (a odd <=199), rd<=2400, H=10000: see batch_passB.log.
* Pass C: rn=+-a^2 (a odd <=999, coprime), rd=240, H=200000 (the "magic
  denominator" law): see batch_passC.log.
* Pass D: ALL coprime members |rn|,rd <= 2000 at H=1000: see
  batch_passD.log.
(Final tallies of B/C/D are appended to data/claude_sib_D_t5_surface.txt.)

## 2. Isogeny/twist neighborhood details

data/claude_sib_D_neighborhood.txt has the 17 codomain minimal models
with exact torsion ([2,12] or [2,2,6], order 24 in all cases) and the
37 twist results.  Notes: (i) the Richelot set (15) is contained in the
TwoPowerIsogenies set (17, two extra at depth 2); (ii) codomain
G2Invariants all differ from the hit (proper isogenies); (iii) since
torsion order 24 groups on simple Jacobians are already realized in the
repo's tables, the neighborhood yields NO new group; (iv) an
interesting negative: no codomain keeps [2,2,2,*] — the full 2-torsion
never survives a (2,2)-isogeny here.

## 3. Enumerator extension (chunks past d=2000)

Driver `sib_D/runchunks_ext.sh`: chunks [lo,lo+99] from 2001, 3 threads,
cached t127.bin, postfilter after each chunk, markers done_dLO_HI.
Tallies (survivors = -P127 stage; pass = after 96|#J postfilter p in
[131,397]; ALL PASS=0 so far => unconditional kill extended):
see sib_D/driver.log and pf_d*.log; summary at return time in the
structured output.  Resume:
  bash /tmp/claude-1000/-home-claude-torsion-jac/3a79ef58-88cb-43d0-927d-10d9ae3e49a0/scratchpad/sib_D/runchunks_ext.sh
(idempotent; skips done chunks; edit the seq bound to extend past 4000).

## 4. Next steps (ranked)

1. **Scale the C_rho' scan** — it is the cheapest known probe of the T5
   hit locus.  Concretely: (a) port qscan to C with threaded member
   loops (current gp: ~1ms/member at H=3000, single thread); (b) push
   pass-A-style boxes to |rn|,rd <= 2000-5000 at H >= 10^4; (c) the
   (u,g) rational chart of the surface S covers members of unbounded
   height at bounded (u,g)-height — coordinate with the lane that owns
   the (u,g) sweep to avoid overlap.
2. **Explain the second u-involution** (u=-97/48 <-> 3637/7105 on the
   same member; both hit members have it): on C_rho' (rank 3) the hit
   points differ by an explicit MW element; identifying it would turn
   the 4-representation fiber into a group-law orbit and possibly
   generate representations on OTHER members (new curves!).
3. **Prove the failing-class constancy law** (class(F) constant on the
   X1X2X3-cover of each member) — it would let one compute kappa(rho')
   a priori and enumerate ONLY kappa=1 members (the hit-capable ones),
   collapsing the member search dimension.
4. Enumerator: resume chunks (see Section 3) toward d=4000-6000.
5. Isogeny neighborhood: (3,3)-isogenies (not available as a Magma
   intrinsic; would need the Bruin-Flynn-Testa machinery) — the 3-part
   is untouched by the 2-power graph, so a (3,3)-neighbor could a
   priori keep order-96 torsion or trade 12 -> 4.

## Addendum (Codex review, PR #4): checked-in scripts
The exact-commands sections above reference the session scratchpad. The key scripts are
now checked in under code/claude_sib_lanes/D/ (same filenames); sweep binaries rebuild
with the gcc lines given in the commands. Scratchpad paths remain valid only on the
discovery machine.
