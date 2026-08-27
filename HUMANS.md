# Orientation for human mathematicians

## The problem

For which finite abelian groups `G` is there a genus-2 curve `C/Q` with
`Jac(C)(Q)_tors ≅ G` — with `Jac(C)` **geometrically simple**, and
separately with `Jac(C)` **geometrically split**? The paper's two census
tables record every group we could either find in the literature and
databases or realize ourselves: 73 groups on the simple side and 76 on the
split side at the time of writing, each row an explicit curve whose exact
torsion (and simplicity or splitness) is certified by a Magma script. The
problem is far from exhausted, and a central hope of publishing this
notebook is that others — human or AI — will push it further.

## What this repository is

A frozen snapshot of the project's working repository: the "lab notebook"
of a loosely coordinated collaboration between four mathematicians and a
rotating cast of AI coding agents (Anthropic's Claude models, including
long-running Claude Code sessions; OpenAI's GPT models — one long-lived
GPT session was nicknamed "Sol" — and the codex review agent, which
reviewed most pull requests). Nothing about the structure was planned in
advance. Each route to each torsion target left behind a triple

    notes/<route>.md   — the design note and running status (the memory)
    code/<route>.m     — the Magma (occasionally Python/Sage/C) scripts
    results/<route>*.log — the captured run logs

and the notes are the part worth reading first: they record what was
tried, what worked, what died, and *why we believed what we believed* at
each dated step — including beliefs that were later corrected (errata are
recorded rather than erased).

## Directory map

| Directory | Contents |
|---|---|
| `notes/` | ~325 dated design/status notes. **The actual notebook.** |
| `code/` | ~1250 scripts: one per experiment; mostly Magma, some Python/Sage/GP/C. |
| `results/` | ~3200 captured run logs and result tables. |
| `data/` | Larger generated data: sieve outputs, certificates, equation dumps. |
| `product/` | The geometrically **split** workstream (gluing constructions, the split census) with its own code/data/logs. |
| `reports/` | Self-contained HTML session reports written by agents for humans — polished summaries of major pushes. |
| `paper/` | The companion documents (realization census, split table, per-group source classification) and `scripts_and_data/`, the paper's verification suite. |
| `skills/` | The distilled lab methodology, written as loadable "skills" for coding agents (see AGENTS.md) but perfectly readable as short methodology surveys. |

## Where to start reading

* `reports/top10-2026-07-25/` and `reports/top10-2026-07-26/` — two
  editions of "the ten most promising avenues", written mid-project as
  strategic surveys. The best single picture of the frontier as we saw it.
* `reports/order96-discovery-session/` — the complete July 17–18 Claude
  Code session, as a readable transcript (`transcript.md` renders on
  GitHub; `index.html` is the styled version), in which the order-96
  group [2,2,2,12] was found overnight: the assignment, the top-10
  ranking, the duel with GPT 5.6's rival list, the launch of the lanes,
  and the morning after. **This is a transcript of one of the many
  human/AI conversations that took place over the project, included as a
  representative example** (the other reports are write-ups, not
  transcripts).
* `reports/z31-order31-writeup/` — the story of the first explicit
  order-31 realization (a real-multiplication curve attached to the
  newform 1830.2.a.q), from observation to certificate.
* `reports/generic-2214/` — how the first generic (`End = Z`) `[2,2,14]`
  curves were found by re-angling a construction that had stalled.
* `reports/split-torsion-table/` and `product/split_torsion_table.md` —
  the split-side census and how its 76 rows were assembled and certified.
* `notes/main_four_target_*.md` — the recurring "state of the union"
  passes over the main open targets.
* `skills/g2-torsion-lab/SKILL.md` — the project's own orientation hub:
  target-naming conventions, the methodology as five reusable moves
  (contact constructions, named charts, halving, 2-rank control, coprime
  composition) plus the search-funnel discipline.

## Unfinished threads (nuggets may be buried here)

The notes document many targets that resisted us, usually with a
diagnosis of *why* (local obstruction vs. global thinness), reusable
charts, and partially explored covers. Examples with substantial recorded
work: cyclic `Z/35` and `Z/48`, `[2,24]`, `[4,16]` (the `m18_416` chart),
`Z/5 × Z/5`, `[2,2,16]`, the `[2,2,2n]` ladder beyond `n = 6`, and the
question of a split `[2,2,6,6]`. The dead ends are documented as
carefully as the successes, and several of the paper's results came from
reviving a route someone else had abandoned.

## Caveats

* Nothing in `notes/` is refereed; claims that matter are the ones that
  made it into the paper, whose computational content is certified by the
  self-contained scripts in `paper/scripts_and_data/`.
* Notes are dated snapshots and are **intentionally not retro-edited**:
  a note may assert something a later note corrects. When in doubt,
  search for later notes on the same route, and consult
  `STATUS_AND_ERRATA.md` at the repository root, which maps the
  significant corrected claims to their corrections.
* Database identifiers from the LMFDB beta site ("alpha labels") that
  appear in notes and logs are not permanent identifiers; curves are
  always recoverable from their displayed equations.
* This is a partial record: it contains what passed through the shared
  repository, not every conversation, session, or private scratch file of
  every collaborator.
