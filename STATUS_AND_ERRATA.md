# Status and errata index

The notes in this repository are dated snapshots and are intentionally not
retro-edited: a note may assert something that a later note corrects. This
index maps the significant corrected claims to their corrections, so that a
reader who lands on a stale assertion can find its final status. The
authoritative record for everything the paper claims is the verification
suite in `paper/scripts_and_data/` (every script asserts what it claims and
was run green); nothing below affects any claim of the paper.

## Corrections to our own working claims

| Working claim (where you may encounter it) | Final status | Correction recorded in |
|---|---|---|
| Conductors computed via PARI `genus2red` (various early notes) | `genus2red` omits the 2-part of the conductor whenever its factorization matrix carries the sentinel row `[2,-1]`; every affected number in the repository was recomputed with Magma (notably the generic `[2,2,14]` witness conductor, `2^2 * 9550095752925 = 38200383011700`) | `notes/claude_ov_erratum_2026_07_25.md` (full statement + repository sweep) |
| "Two new [2,22] curves" (2026-07-31) | Retracted 2026-08-02: the curves are isomorphic to previously known ones | ERRATUM at the top of `notes/nonrm_222_hunt_2026_07_31.md` |
| The displayed `[2,2,4,4]` split witness (mid-August table drafts) | The table briefly displayed a curve that differed from the verified witness ("witness drift"); fixed, and both tables re-verified row-for-row on 2026-08-22 | `notes/claude_paper_computation_audit_2026_08_22.md`, item D1 |
| Cyclic `[40]` claimed as a new record (July drafts) | Retracted: Elkies (2001–02) already exhibited cyclic orders 34, 39, 40; our curve stands as an independent realization with new strict certificates, and the cyclic record credit is Elkies' | `notes/claude_top10_2026_07_25.md` (corrected framing); the paper's Source column |
| `[6,12]` p-adic (Coleman) computation bundle (2026-07-13) | Withdrawn when review found printed digits at unproven precision ("printed zeros are not certified zeros"); superseded by precision-safe drivers | `notes/m612_review_and_top3_plan_2026_07_13.md` (withdrawal recorded in the same note) |
| Phi38 "genus-7 verdict" (Z/35 lane, 2026-07-18) | Corrected to "genus ≥ 7": two finite-field computations give the lower bound; the exact-over-Q computation never completed (an 11-hour runaway job, terminated) | `notes/claude_prod_04_35.md` (post-review correction) |
| Split-census "new group" announcements (2026-08-13) | Several were new only relative to the production-LMFDB + HLP baseline, not relative to the extended database; corrected in the same note, and the final Table 2 reflects the corrected accounting | `notes/split_census_2026_08_13.md` (correction sections) |
| An early "infinite family" claim from a single seed (the kappa = 0 tale) | Retracted: the seed lay on a spurious perfect-square degeneration branch; codified as the "never claim a family from its seed" rule that later hits were held to | `skills/validate-and-record-a-hit/SKILL.md` §7; `notes/agent_a2_24_d0_derivation.md` |
| The earliest `(2,2,4,4)` candidate bank (`tor2244.txt`, 26,653 rows) | The original bank had a generation bug; later stages regenerated candidates independently and do not rely on it | noted in `notes/target_22412_second_stage_2026_07_18.md` |
| `[31]` novelty phrasing (late July) | Scoped precisely after review: ours is the first *explicit* order-31 realization (and the first on a genus-2 Jacobian); existence on some inexplicit member of the isogeny class follows from Alessandrì–Coppola Thm 4.2/4.3 | `notes/claude_z31_rm_witness_2026_07_30.md`; `reports/z31-order31-writeup/` |
| Early half-divisor derivation (GPT-5.4) | Contained an error; rederived correctly (second version) | `notes/derivations/README.md` (both versions kept deliberately) |

## Corrections to the literature found along the way

* **Bernard–Leprévost–Pohst 2009, row C4**: the printed parameter row
  contains a misprint; the unique consistent fix is `b = -277/243`, which
  yields a genuine (previously unreported) second explicit `[2,22]` curve.
  See `notes/claude_generic_222_2214_plan_2026_07_23.md`.
* **PARI/GP `genus2red` 2-part omission** (software behavior, documented
  for p > 2): see `notes/claude_ov_erratum_2026_07_25.md` for the exact
  failure signature and the workaround used repository-wide.

If you find a stale claim not listed here, check for later notes on the
same route (`ls notes/ | grep <route>`) before trusting it, and treat the
paper plus `paper/scripts_and_data/` as the final word.
