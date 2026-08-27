# Review of ChatGPT_summary.tex (GPT 5.6 Sol) — 2026-07-20

(The reviewed draft — an early survey that the paper was ultimately
written independently of — is not included in this snapshot; the claims
checked are restated below.)

Companion to the earlier review of Sol's ranked plan (notes/claude_review_gpt56_plan.md).
Verification script: code/claude_sol_review_verify.m, log results/claude_sol_review_verify.log.

## Verdict

A good, careful frontier report. Its epistemics are exactly right (marked vs exact vs
isolated jump; irreducible Frobenius != geometric simplicity; positive-dimensional locus
!= infinitude; bounded search != nonexistence), and every numeric claim I spot-checked is
correct: record-curve reductions #J(F_31)=864, #J(F_37)=1248 (gcd 96, cross-prime checks
valid), P_37/P_73 as printed, C_66 exactness pattern (36/144). The two-prime gcd exactness
argument is clean and worth adopting in the paper as the black-box-free certification.

## The one substantive gap found — now fixed

Sol calls the cyclic [40] (contact-5 family, t=-1/3) and [28] (contact-7 rational-root)
jumps "exact geometrically simple", but the underlying notes recorded only Q-simplicity
(irreducible Frobenius at p=7 resp. p=5). Worse, p=5 for the [28] curve has
chi_5 = x^4+2x^2+25 (a_p=0), so the root-power test cannot even run there.
NEW CERTIFICATES (this session): both curves pass the strict root-power test
([Q(pi^n):Q]=4 for 2<=n<=12):
  [40]: y^2 = -324x^5+1296x^4+1944x^3-5103x^2-4374x+6561, exact [40], STRICT at p=17, 23.
  [28]: y^2 = 4x^5+21x^4-70x^3+79x^2-42x+9, exact [28], STRICT at p=7, 31.
ERRATUM (2026-07-20, same day): an earlier version of this note claimed the [40] curve
beats Elkies's cyclic [32] as a record. WRONG: Elkies's 2002 page ALREADY exhibits an
order-40 curve y^2 = (3x+4)(x^4+5x^3+8x^2+(19/4)x+1), "the highest known (as of 6/2001)
for a simple genus-2 Jacobian", plus curves of orders 34 and 39 (and a rational
31-SUBGROUP with no rational generator). So cyclic [40] is Elkies 2001/2002; our
contact-5 curve is an independent realization of the KNOWN group, now carrying strict
geometric-simplicity certificates. Table 1 of NotesAndTodo.tex correctly credits [40]
and [34] to elkies2002; its [39] row is MISSING the elkies2002 credit and should get it.

## Staleness (not errors)

Predates: curve #2 of (2,2,2,12) (07-18), curve #3 + the multigrade-variety theorem
(07-20, notes/claude_ari_surface_22212.md), the T5-pencil structure theorems, and the
current ranked program (notes/fable_final_top10_2026_07_19.md). Its open-directions list
([49],[60],[63],[2,24],[2,6,6],[2,2,4,8]) is far weaker than the live top-10
([8,8],[2,24],[3,12],[2,2,24],[6,12],[2,2,16],[4,12],...).

## Does it suggest new paths? Mostly no — the repo already tried them

Checked every implicit upgrade against the notes:
- [2,42]-type (contact-7 two-root +3): bounded-negative, hits confined to a 5-adic
  boundary disk (notes/contact7_two_root_plus3_search_2026_07_11.md).
- [46] (Kuru-Sadek [23] + branch point): condition curve has geometric genus 16;
  quotient/boundary analysis open (notes/order23_extra2_quotients.md).
- [45] (contact-9 + 5): direct chart obstructed at p=3, boundary-confined
  (notes/contact9_family.md).
- [2,32] (Elkies-32 + extra root): height-5000 negative, p=13 blocks the extra-root
  condition (notes/elkies32_extra2.md).
- 80-via-[40]: the order-8 half does not divide at t=-1/3 (contact5_order40_family.md).
- 56 = 8x7: contact7_family.md already has the target-56 subfamily section.

## Genuinely actionable items extracted

1. DONE: strict certificates for [40] and [28]; promote [40] as the cyclic record.
2. Port the Elkies 16D=T Riemann-Roch rank-drop machinery (already reconstructed in-repo
   for [32]: code/elkies32_* ) to the Z/48 = 16x3 lane: encode the 16-chain as a rank
   drop + independent 3-contact, instead of iterated halving square-conditions. The Z48
   notes only ever tried halving-chain scans; this is a different algebraic encoding and
   the one fresh strategy this review surfaced.
3. Audit what else the Kuru-Sadek quadratic-order machinery yields at genus 2 beyond
   [23] — the repo only verified their printed member.
4. Standing rule applied: depth-3 Richelot + divisibility audit of curve #3 launched
   (code/claude_curve3_audit.m) to test the T5-component-rigidity conjecture on the
   first post-conjecture hit.
