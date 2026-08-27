# Review of notes/ranked_unrealized_torsion_plan_2026_07_17.md (GPT 5.6 codex)

Date: 2026-07-17. Method: six verification agents read every underlying note the plan cites
(agent_z5x5_*, agent_Z35_*, agent_A2_24_*/agent_a2_24_*, agent_Z48_*, agent_m18_416_*,
contact30/order-40/2220/contact9 notes, contact6_m612_* campaign) and audited each claim.
Companion: claude_top10_ranking.md + claude_top10_00_summary.md (same date).

## Overall verdict
A competent, well-structured survey of the committed 2026-07-02 "main_four_target" program —
its top four ([5,5], [35], [2,24], [48]) are exactly that program's four lanes — with a genuinely
useful certification checklist. But its confidence calibration is systematically inflated (it
reads "no local obstruction" as positive evidence, the repo's classic locally-alive/globally-dry
signature notwithstanding), it contains one outright misreading, one plan item directed at a
proven-impossible lane, several already-built or already-refuted plan steps, and one scoping
error that wastes a top-10 slot. It also omits, without discussion, most of the strongest
open targets by the project's own evidence standards.

## Verified findings, per target
- **[5,5] (rank 1, "high")** — VERIFIED FACTS: the b2=0 full-norm slice really does have smooth
  F_7/F_11 points with unique Hensel lifts to 7^7/11^6. OMITTED by the plan: the CRT
  height-frontier theorem (any rational point through the studied charts needs coordinate
  height > 10^6; all 122,544 CRT-consistent tuples died at the first sieve prime); the F_3 open
  locus is empty; and TWO sibling (5,5) routes (contact5×contact5, y−H degree-2) were proven
  empty over Q by exact algebra despite abundant F_p points. The saturated slice is
  ZERO-dimensional, so its Qbar points are finitely many fixed algebraic numbers — the smooth
  mod-p points are probably simply irrational, and the plan's Coppersmith step attacks them
  anyway. Fair confidence: MEDIUM at best. The decisive (cheap) step is the over-Q Groebner
  decomposition of the slice, which the plan treats as a stepping stone but is a decision
  procedure. HOWEVER: the target promotion itself is partly right — (5,5) has order 25, the
  smallest open group (a fact my own ranking missed; see correction below).
- **[35] (rank 2, "medium-high")** — branch data accurate (18 → 2 central branches certified to
  scaled 3^8). But: plan step 1 would REBUILD the compressed automaton that already exists
  (code/agent_Z35_compressed_automaton.m) using a state the notes PROVED non-Markov; the lift
  tree is being pruned (multipliers 27,…,27,3), which is accumulating obstruction, not
  candidates; no certified Z_3 point exists; and the plan ignores the complementary A_1(5)
  result (no local obstruction on the full threefold; F_3 supports #J=35 only via
  (x^2+x+3)(x^2+3x+3)). Right rank, partly wrong reasons; the two routes are consistent and
  the A_1(5) production sieve is the better next move, with the symbolic obstruction
  derivation (plan step 2) the one genuinely good contact-lane item.
- **[2,24] (rank 3, "medium")** — the four-fiber closure and height-5 negative are accurately
  reported, but "mechanically clear" contradicts the notes' own verdicts ("high-height /
  construction-only"; "record simple [2,24] as open/hard"); the proposed mod-2 prefilter
  already exists twice; and the plan omits the composite lane's results (two certified simple
  Z/24 at H=12; 991,275 W-split curves with zero 24s — 3-torsion anti-correlated with 2-rank 2).
  The only route any note calls live is the M(12)-surface third component (25 unfitted points).
  Fair: medium-low, rank ~5-6.
- **[48] (rank 4, "medium")** — RTHeight-4 coldness and the gate-fooling incident are accurate,
  but 100% of in-lane evidence is negative (the 16-part and 3-part have NEVER co-occurred in
  ~4.4M tests plus reverse-order scouts), the repo's own fifth-pass had demoted Z/48 to
  "bounded background", and the plan omits Howe's split C48 family — the only real existence
  evidence, and (ironically) the strongest argument [48] has. Fair: low-medium, rank 8-10.
  Cheap decisive diagnostic: sweep the A16 equations over F_5/F_7 for 3|#J co-occurrence.
- **[60] (rank 5, "medium")** — MATERIAL ERROR: plan steps 1-3 and work-order item 4 direct
  effort into the order-30→60 halving lane that contact30_to_60_halving_obstruction.md closes
  with a COMPLETE local impossibility theorem over Q_11 ("not a bounded-search failure").
  The extra-2 lane is closed by unconditional theorems. Only M(12)+5 and order-20+3 are live,
  both with deep bounded negatives (2.39B pairs; h=100000 sieve). Fair: low-medium.
- **[8,8] (rank 6, "medium-low")** — CATEGORY ERROR in the stated rationale: the "[8,8]
  degree-pattern rows" in the A(2,24) scans are multisets of irreducible-factor degrees of a
  degree-16 halving resultant (16 geometric order-24 halves in two Galois orbits of size 8) —
  they carry zero information about Z/8×Z/8 torsion; if anything they are evidence of
  difficulty for [2,24]. Its suggested M_1(8,2,2) substrate has only ONE 8-chain and
  known strongly-locally-constrained halving (0/6,048,375 local survivors at h=10). The level
  is defensible on the real evidence (split P^2 family; no local obstruction to 43), but the
  route should be the Nicholls Prop 5.9.6 (4,4)-family (validated 2026-07-17).
- **[4,16] (rank 7, "medium-low")** — the plan is two generations behind: its step 1 is the
  route run since 07-02 (rigorous negative to height 800, corrected all-twist kill tables,
  suspected fiberwise Sha-type obstruction, R=-25/4 and R=-29/8 fiber-closure certificates),
  and its step 2 (Elkies-[32] subloci) is disfavored by committed negatives ([2,32]: 0/30.4M
  with a structural p=13 obstruction; the printed member's order-32 classes provably don't
  halve). The same-day per-R elliptic solve makes production search ~4000x cheaper and already
  produced a new certified-simple [4,8] at height 5093. [4,16] deserves a HIGHER rank than 7.
- **[2,2,20] (rank 8, "medium-low")** — SCOPING ERROR: realized. Exact torsion [2,2,20]
  (order 80) at the contact-5 point t=-8233/7225, Lombardo + p=71 12th-power certificates
  checked in (how_we_found_2220_examples.md); listed in the NotesAndTodo.tex known-groups
  table (Filip's email 6/28/26); Chabauty even proves the seed unique in its locus. At most a
  certificate-bundle packaging task; the slot should go to a genuinely open group.
- **[72] (rank 9)** — claims roughly accurate; omits that the contact-9 root subfamily's own
  note concludes it "appears to stop at [18]" (0 halvings in 15,718 fibers). Low-medium is
  defensible; the halving-cover formulation is genuinely new work.
- **[80] (rank 10)** — "order-40 family" does not exist: the order-40 locus is ONE isolated
  curve (t=-1/3) on a high-genus 41-node cover, and that curve's order-8 class provably does
  not halve; [2,40] and [40]+3 are boundary-forced and empty to h=3000. [80] is known for NO
  abelian surface construction at all (not even products). Fair: low; bottom rank correct,
  framing not.

## Omissions
The plan silently omits (3,12), (2,2,2,12), (2,2,4,8), (2,6,6), (2,2,16), (2,30), and the
[6,12] lane it sits two commits above. For (3,12) and (2,2,2,12) no demoting negative exists
anywhere in the tree — omitting them while ranking [72]/[80] is the plan's largest
prioritization error. (Its top-4 = the committed 07-02 program; the untracked same-day
claude_top10_* dossiers were invisible to it, which excuses missing the newest upgrades but
not the committed-evidence inconsistencies.)

## Correction to my own ranking (credit where due)
claude_top10_ranking.md called (3,12) "the smallest open group anywhere" (order 36). Wrong:
(5,5) (order 25), Z/31, and Z/35 are smaller open groups. The small-order heuristic is a
genuine argument for (5,5) that my list missed entirely — GPT's #1 pick has real merit as a
TARGET even though its "high" confidence and its specific plan do not survive scrutiny.

## Revised merged top-10 (2026-07-17, post-review)
1. **(3,12)** — unchanged: split-abundant, no negative anywhere, carrier surface now known
   genus-4-fibered, new p=2/p=5 conductor theorem.
2. **(2,2,2,12)** ↑ from 4 — the T3-pencil correlation mechanism (three of four conditions
   collapsing to one positive-rank elliptic curve per pencil member, ranks exactly 2 and 3)
   is the strongest new positive signal in the campaign; hit-existence is now a finite
   symbolic question (rho-scan).
3. **(4,16)** — per-R solve built, wall down, new in-family point at height 5093; twisted-class
   conic and M_1(8,2^w) deformation untried.
4. **Z/35** ↑ from 10 — two independent routes now GO (A_1(5) threefold local landscape +
   the 3-adic central branches), convergent and consistent.
5. **(5,5)** NEW (from GPT 5.6's list) — smallest open group (order 25); one live route
   (b2=0 full-norm branch) with smooth p-adic charts; tempered by the >10^6 height wall,
   two exact sibling-route kills, and zero rational structure to date. First step is the
   over-Q slice decomposition (a decision procedure, hours).
6. **(2,24)** ↓ from 2 — all three programs (A(12) descent, composite/W-split, M(12) main+G1)
   end in demotion; hope rests on the unfitted third component.
7. **(2,2,4,8)** — flagship unchanged: twisted second component now enumerated to d<=3000
   (0 hits, statistics justify 2-3 more orders), HPL threefold and genus-3 split-locus curve
   (elliptic Chabauty over Q(sqrt 2)) still live.
8. **(2,30)** ↑ — 6x5-combination route GO on A_1(5) (room at every prime; five new simple
   [2,10] curves banked); 3-part must be imposed algebraically (cubic contact on q2=-1/4).
9. **(8,8)** — Nicholls (4,4)-substrate validated; decide by the stage-1 cover genus gate.
10. **(2,6,6)** ↓ — split-locus confinement certified to H=150 on the contact-6 chart;
    stays only on the strength of the untried sqrt3-RM escape route.
Near-list: (2,2,16) (suspected Brauer-Manin; parity audit pending), [48] (Howe split family
= only real evidence), [60] (two live lanes only), (6,12) (Prym(E8/E4) rank pending; also
subsumed by #1), (2,2,24), Z/31/37/38, [72], [80].
