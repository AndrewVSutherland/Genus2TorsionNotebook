# Skill authoring convention (read before writing any SKILL.md)

This directory is a **skill library** for the geometrically-simple genus-2
Jacobian torsion project. Its audience is junior/mid engineers and
Sonnet-class models who must debug, extend, validate, and advance the project
**without the original authors present**. Optimise for *correctness* and
*actionability*, not brevity.

## File format

Each skill is one directory `skills/<skill-name>/SKILL.md`. Optional
supporting files (reference tables, template `.m` scripts) may live alongside
`SKILL.md` and be linked from it.

`SKILL.md` begins with YAML frontmatter, then markdown:

```markdown
---
name: <skill-name>            # kebab-case, matches the directory
description: <one sentence>    # WHAT it covers + WHEN to use it; this is how a
                               # model decides to load the skill, so be specific
                               # and include trigger words (target names, error
                               # symptoms, technique names).
---

# <Title>

## When to use this
...
## <body sections>
...
```

## Quality bar (non-negotiable)

1. **Ground every claim in the repo.** Every formula, chart, function name,
   file path, and numeric result must be traceable to a real file in
   `code/`, `notes/`, `results/`, or `paper/`. Cite the file (e.g.
   `see code/torsion_cover_lab_utils.m`). If you cannot verify a claim from
   the repo, mark it `UNVERIFIED` rather than asserting it.
2. **No invented Magma APIs.** Only use Magma intrinsics that appear in
   existing scripts or that you have confirmed. When in doubt, copy the idiom
   from a real script and cite it.
3. **Show, don't just tell.** Include the actual code idiom / command line a
   reader would run, copied or adapted from a real script, with the file it
   came from.
4. **Cross-link.** Reference sibling skills by name in a `## See also`
   section so a reader can navigate. The hub skill is `g2-torsion-lab`.
5. **Pitfalls section.** Each skill ends with a `## Pitfalls` section listing
   the concrete ways a cheaper model will get this wrong, with the fix.
6. **Self-contained.** A reader who loads only this skill (plus the hub)
   should be able to act. Do not assume they have read the whole `notes/`
   tree; summarise the load-bearing facts and point to the note for depth.

## House facts (true across the project — state them, don't re-derive)

- Everything runs in **Magma** (invoked as `magma` or
  `magma`). No Sage/Pari except the historical
  Sage/Lombardo endomorphism check.
- Curves are genus 2, `y^2 = f(x)` with `deg f in {5,6}`; work on the
  Jacobian `J`. "Simple" always means **geometrically simple** (`End(J_Qbar)`
  has no nontrivial idempotent), which must be *certified*, not assumed.
- Operational limits: **at most 3 concurrent Magma jobs**; the machine has
  OOM'd at 6. Cap memory with `SetMemoryLimit(n*10^9)`. Git commits are made
  locally; **the user pushes** (agent push fails on auth).
- Target naming: `m` followed by the invariant-factor digits, e.g.
  `m2220` = `Z/2 x Z/2 x Z/20`, `m36` = `Z/6 x Z/6`, `m3222` = the `[2,2,16]`
  worker line, `a2_24` / `A(2,24)` = `Z/2 x Z/24`, `Z35` = `Z/35`,
  `Z48` = `Z/48`, `m18_416` = the `M_1(8,4)` chart for `[4,16]`,
  `contactN` = the order-`N` contact construction.

## Process for authors

- Read `AUTHORING.md` (this file) and the hub `g2-torsion-lab/SKILL.md`
  first, then the specific canonical files named in your task.
- **Do not run Magma during authoring** (to avoid exceeding the concurrent-job
  limit while several authors run in parallel). Ground statements by reading
  files. Validation runs happen in the review phase.
- Write the `SKILL.md`, then return a short summary: what you wrote, which
  files you grounded each section in, and anything you could not verify.
