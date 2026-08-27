# Orientation for AI agents

You are looking at the frozen lab notebook of a finished-but-not-exhausted
research project: realizing finite abelian groups as `Jac(C)(Q)_tors` for
genus-2 curves `C/Q`, in both the geometrically simple and geometrically
split settings. The paper (Balakrishnan–Najman–Shnidman–Sutherland,
*Rational torsion on simple genus two Jacobians*) reports the results; this
repository is the record of how they were found, by a mixed team of humans
and coding agents, and it was published partly **so that agents like you
can continue the work**.

## Orient in this order

1. `skills/g2-torsion-lab/SKILL.md` — the project's orientation hub:
   the goal, the target-naming convention, the methodology as reusable
   moves, and an index of the other skills. Written mid-project for
   incoming agents; some counts and statuses are dated, but the
   conventions and methodology are accurate.
2. The relevant target's notes: `ls notes/ | grep -i <target>`. Every
   route has a dated design/status note; **read the note before the
   code**, and prefer the newest note when several cover one route.
3. `reports/top10-2026-07-25/` and `-26/` — strategic surveys of the ten
   most promising avenues, as of those dates.
4. The paper's verification suite `paper/scripts_and_data/` (self-contained
   Magma; `README.md` there explains the table data files) — the ground
   truth for every claimed curve.

## Conventions you must know

* **Target naming**: `m` + concatenated invariant factors
  (`m2220 = [2,2,20]`, `m36 = [6,6]`, `m18_416` = the `M_1(8,4)` chart for
  `[4,16]`); cyclic targets `Z35`, `Z48`; `a2_24 = [2,24]`;
  `contactN` = the order-`N` contact construction. Curve models are
  `y^2 + h(x) y = f(x)` with coefficient lists **ascending** (LMFDB
  convention) when serialized as `[[f],[h]]`.
* **Route triple**: `notes/<route>.md` + `code/<route>.m` +
  `results/<route>*.log` (larger outputs in `data/`). The split-Jacobian
  workstream lives under `product/` with the same structure.
* **Search-funnel discipline** (see `skills/running-torsion-searches`):
  parametrize → finite prefilter over `F_p` → exact `TorsionSubgroup` →
  simplicity certificate. Logs use greppable markers (`PROGRESS`, `HIT`,
  `SEARCH_DONE`).
* **Nothing is a result without both certificates**: exact torsion AND
  geometric simplicity (root-power/Frobenius or Zywina-style; see
  `skills/simplicity-certificates` and `skills/validate-and-record-a-hit`).
  A "hit" that turns out split (on the simple side) is worthless.
* **Magma pitfalls**: read `skills/magma-lab-conventions` before writing
  or debugging any Magma in this style.

## The skills

`skills/` contains the lab's distilled methodology in Claude Code's skill
format (a directory per skill with a `SKILL.md`). If your harness supports
project skills, you can install them by placing the directory at
`.claude/skills/` in your working copy; otherwise just read them — they
are self-contained prose. `skills/REVIEW_FINDINGS.md` documents the
verbatim-quote audit that checked the skills' formulas against the code.

## What is authoritative vs. exploratory

* Authoritative: `paper/scripts_and_data/` (every script asserts what it
  claims and was run green), the companion documents in `paper/`, and any
  note explicitly recording a certificate.
* Exploratory: everything else. Notes are dated beliefs; some were later
  corrected (corrections are recorded in later notes and `ERRATUM`
  entries, never silently; `STATUS_AND_ERRATA.md` indexes the significant
  ones). Always search for the latest note on a route before trusting a
  status.
* Safety note: the scripts here are frozen research code that assume
  trusted inputs (several `eval` their data files by design). Run them
  only on this repository's own data, in an isolated environment.

## Known dangling references (deliberate)

This public snapshot has a fresh git history. Notes and reports may
reference pull-request numbers, commit hashes, branch names, private
working paths (e.g. machine-specific prefixes), scratch directories, or
files that only existed in the private working repository or on a
particular machine. Treat these as historical color, not as links to
resolve. Database "alpha labels" (LMFDB beta site) are snapshot
identifiers; recover curves from their equations
(`https://alpha.lmfdb.org/Genus2Curve/Q/?jump=[[f],[h]]` looks a curve up
by equation, matching isomorphic models).

## Open problems (as of the snapshot)

Targets with substantial recorded work and no realization at snapshot
time include: cyclic `Z/35`, `Z/48`, `[2,24]`, `[4,16]`, `Z/5 × Z/5`,
`[2,2,16]`, the `[2,2,2n]` ladder beyond the realized cases, and split
`[2,2,6,6]`. For each, the notes record the routes tried, the reason each
stalled (local obstruction vs. thinness — see `skills/local-obstructions`),
and often a concrete "next thing to try". The strategy notes
(`notes/main_four_target_*.md`) and the two top10 reports are the best
starting points. The guiding open question of the paper: are the special
loci where large torsion keeps appearing genuinely richer, or just easier
to search?

## If you continue this work

Fork or copy; this repository is a frozen artifact. Keep the conventions
(note-first, route triples, dated notes, certificates before claims) —
they are what made it possible for one agent to build on another's work.
