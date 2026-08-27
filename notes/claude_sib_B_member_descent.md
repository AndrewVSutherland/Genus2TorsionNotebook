# Lane B: same-member Mordell-Weil descent on the T5 pencil — (2,2,2,12) siblings

Date: 2026-07-18.  Sibling-hunt agent, Lane B (member MW descent).
Scripts/logs: scratchpad `sib_B/` (recon.m, sibenum.m, g2cmp2.m, sanity.gp,
j*.log).  Inputs: notes/claude_prod_02_22212.md,
data/claude_prod_02_22212_structures.txt, prod_02_22212/ scripts.

## 0. Headline

1. **The lift5.m blocker is solved.**  Descent-based `Generators` hangs even
   on the minimal models of these covers; the fix is to AVOID descent
   entirely: harvest small rational points on the quartic cover itself
   (`Points(C : Bound := 1e5)` — every live cover has 20-80 of them), push
   them through cover -> E -> MinimalModel(E), then `ReducedBasis` +
   `Saturation(...,11)`.  This attains the full rank (verified against
   RankBounds on every cover run) in < 2 s per cover.
2. **A THIRD representation of the known hit curve was found and fully
   verified**: u = 3637/7105 on the hit member rho' = -49/240,
   (s,m,n) = (-527365, 174576, -1030225),
   B = [12923079325, -442594240440, 100710366685, 179852559600, -880910496000],
   TorsionSubgroup = [2,2,2,12] exact, simplicity certificates p=37,73,113,
   **G2Invariants identical to the known hit** (g2cmp2.m) — same curve, new
   model.  u-height 7105 explains why the old member scan (u-height <= 4000)
   missed it.  NOT a fiberwise 2:1 identification: a collision scan
   (pairfind.m) over 86 small u on the member found NO other G2-equal pair,
   so the double intersection of the (2,2,2,12) locus with this member is a
   special coincidence of the hit curve (analogous to its cross-member second
   representation), the hit curve now having >= 3 pencil representations.
   Every future FULL candidate must be G2-deduped before being called new.
3. **All live-member cover data** (recon.m): on every hit/near member BOTH
   quartic covers W3, W4 have square leading coefficient (lc(W4) =
   4(2rho'-1)^2 identically; lc(W3) = 16rho'^3(rho'-1), square exactly on the
   +-odd^2-numerator members), hence 2 rational points at infinity = free
   origin.  Minimal models and ranks:

   | member | cover W3 | cover W4 |
   |---|---|---|
   | -49/240  | cond 512295, rank 2, Z/2 | cond 940695, rank 2, Z/2 |
   | 289/240  | cond 940695 (= partner W4!) | cond 512295 (= partner W3!) |
   | -1/143   | cond 98670, rank 1, Z/4 | cond 2882880, rank 2, Z/2 |
   | -25/551  | cond 188359350, rank 2, Z/2 | cond 94419360, rank 2, Z/2 |
   | -169/1431| cond 558110670, rank 2, Z/2 | cond 127327200, rank 1, Z/2 |
   | 841/697  | cond 932021430, rank 3, Z/2 | cond 9372363840, rank 1, Z/2 |

   The pencil involution swaps (member, W3) <-> (partner member, W4) at the
   level of elliptic minimal models — enumerating both covers of ONE member
   of an involution pair covers both members' moduli.
4. **G-cover trick**: with G(u) := (q rho' - 2u + 1)(q rho' - 1) (the proven
   square-class of W1 W2 W3 W4), lc(G) = 16 rho'^2 is always square, so
   z^2 = G(u) is a third elliptic probe of the same member; on the hit member
   it has rank 3 (gens 4.41/13.66/4.79, cond 1619876790).
5. **NEW PROVEN IDENTITY (conid2.gp, symbolic in generic rho'):**
   `W1*W2 = rho'(rho'-1) * (4 rho' u^2 - (6 rho'+2) u + (3 rho'+1))^2`.
   Consequences:
   * **DEAD-member certificate: a T5 member can carry a hit ONLY IF
     rho'(rho'-1) is a rational square** (else the two conics are never
     simultaneously square).  This PROVES the empirical "+-odd^2 numerator"
     law of the discovery note.
   * **The live locus is the rational curve rho' = q^2/(q^2-p^2)**
     ((p,q) = (17,7) -> -49/240, (7,17) -> 289/240, (12,1) -> -1/143,
     (24,5) -> -25/551, (40,13) -> -169/1431, (12,29) -> 841/697) — the
     exact T5 analog of the T3 w-locus 1-2rho = w^2.  Fresh live members at
     small (p,q) ((2,1) -> -1/3, (3,1) -> -1/8, (3,2) -> -4/5, ...) are
     virgin sibling territory (previously only box-swept to u-height 150);
     smallmem.sh probes them with the lattice machinery.
   * On live members W1 <=> W2 (density-1/2 "conpass" on the W4/G covers),
     and since W1 W2 == square: **G == W3*W4 mod squares — any two of
     {W3 sq, W4 sq, G sq} imply the third**.  A member hit is exactly:
     u on cover W3 AND on cover W4 AND one conic (W1) square.  Full hits =
     rational points of the fiber product C3 x_u C4 (higher genus, Faltings-
     finite) with one extra conic condition — explains hit scarcity per
     member and makes Chabauty on the fiber product the closure route.
   * Also W1*W4 == W2*W4 == rho'(2u-1)(2 rho' u - rho' - 1) mod squares
     (a per-member quadratic class K(u); on -1/143's W4-lattice K passes on
     ~every point — j4 log, unexplained correlation worth a proof attempt).
6. No genuinely NEW (2,2,2,12) curve so far (see tallies below); the deep
   lattice tails are exponentially thin, as expected from the fiber-product
   picture.

## 1. Machinery (sibenum.m)

`magma -b rn:=RN rd:=RD cv:=CV nmax:=N nm2:=N2 tmax:=T sibenum.m` with
cv in {3,4,5} (5 = G-cover; tests all four W's).  Pipeline: integral quartic
model -> origin at infinity -> EllipticCurve(C, org) -> MinimalModel FIRST ->
RankBounds -> small-point harvest -> ReducedBasis -> Saturation -> drop
height-0 torsion strays -> lattice enumeration (rank 1: |n| <= nmax
incremental; rank 2: (2*nm2+1)^2 box, center-out rows so a time cap only
loses far shells; rank 3: |n_i| <= 10 box), each point mapped back to u via
Inverse(MinimalModel map) then Inverse(quartic->E map) with try/catch, then
exact IsSquare of the other conditions.  FULL candidates print exact u; the
degenerate u=1 (m=0, B4=0, degree drop) always appears and is discarded.
Validations: known hit u=-97/48 and partner rep u=133/145 recovered on the
hit members; known near-misses u=17, 13/4 (-1/143) and u=10 (-25/551), u=-4
(-169/1431), u=43/52 (841/697) recovered on theirs.

## 2. Enumeration tallies

(see final section of this note and j*.log for exact counts; filled at end
of session)

## 3. Files

* scratchpad sib_B/: recon.m (cover survey), sibenum.m (enumerator),
  g2cmp2.m (dedupe + reduced model), sanity.gp (W-formula validation),
  verify_3637_7105.log (jackpot protocol of the third representation),
  j*.log (per-cover enumeration logs), smoke.log.
* data/claude_sib_member_closures.txt — per-member/cover enumeration
  certificates (depth statements; rank-0 closures if any).

## 4. Next steps

1. Finish/extend the rank-2 boxes (nm2 > 100) and rank-1 chains (nmax > 400)
   if more budget; the G-cover (cv=5) rank-3 box deserves a deeper
   enumeration (its own N3 is hardcoded 10).
2. The 2:1 moduli identification on a member (u=-97/48 <-> u=3637/7105)
   should be made explicit (Richelot/model isomorphism as a map u -> u');
   quotienting by it halves all future member searches and may expose the
   member's true "new-curve" locus.
3. The fiber-product curves (W4-cover x_u G-cover) per member are explicit
   genus <= 5 curves; TwoCoverDescent/Chabauty on them could PROVE a member
   has no second curve, upgrading depth statements to closures.
4. Other T5 members with rd = 240-type structure (rho'-numerator +-odd^2):
   generate fresh members and run recon.m + sibenum.m on them (the whole
   pipeline is ~1 min per member/cover at nm2=40).

## Addendum (Codex review, PR #4): checked-in scripts
The exact-commands sections above reference the session scratchpad. The key scripts are
now checked in under code/claude_sib_lanes/B/ (same filenames); sweep binaries rebuild
with the gcc lines given in the commands. Scratchpad paths remain valid only on the
discovery machine.
