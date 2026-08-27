# The extended-database torsion histogram (finally run)

Date: 2026-07-25.  Run by the orchestrator directly against the LMFDB mirror
once MCP access came back; this closes avenue 4 of
`notes/claude_top10_2026_07_25.md`.

## Why this needed doing

The 2026-07-22 census only *refined witnesses for groups already on the
collaborators' hand-curated `NotesAndTodo.tex` Table 1* — the query was
`... WHERE is_simple_geom AND geom_end_alg='Q' AND torsion_subgroup IN (...)`
with the `IN`-list inherited, not measured.  So "absent from the database" was
an assumption for every open target except the 11-slice.  No histogram over the
extended Booker–Sutherland collection had ever been computed.

## The query

```sql
SELECT torsion_subgroup, count(*), min(cond)
FROM g2c_curves_new WHERE is_simple_geom
GROUP BY torsion_subgroup ORDER BY 2 DESC;
```

**65 rows** (64 nontrivial groups).  Production `g2c_curves` gives 36 rows for
comparison.

## Result 1 — nothing is hiding in the database

Diffed against the 71 rows of `paper/torsion_realizations.tex`:

```text
In the DB but NOT a paper row:  (none)
```

The paper's table is **complete with respect to the extended database**.  Every
open target — `[2,24]`, `[3,12]`, `[4,12]`, `[6,12]`, `[8,8]`, `[2,2,16]`,
`[2,2,18]`, `[5,5]`, `[35]`, `[45]`, `[63]`, `[7,7]`, `[5,10]`, `[2,2,24]`,
`[2,2,4,8]`, `[4,16]`, `[2,2,2,14]`, `[2,2,2,24]`, `[2,2,4,12]` — is now
**confirmed absent by a measured query** rather than by inheritance.  The
epistemic risk flagged in the review is retired.

Exactly seven paper rows are absent from the DB, i.e. are genuinely
construction-only:

```text
[6,6]   [2,4,8]   [2,2,20]   [2,2,2,8]   [2,2,4,4]   [2,2,2,10]   [2,2,2,12]
```

(the first five plus the record are this project's; `[2,2,2,10]` is Elkies 2024).

## Result 2 — an attribution error in the paper

`[2,4,4]` is marked **"new"** in `paper/torsion_realizations.tex`, but the
extended database contains **four geometrically simple curves with exact torsion
`[2,4,4]`, all with `geom_end_alg = 'Q'` (generic)**:

```text
172260.bj3   cond 172260
180180.cg7   cond 180180
727650.qf8   cond 727650
956340.fc7   cond 956340
```

By the paper's own stated conventions — display the smallest-conductor known
generic witness, and cite `\cite{BookerSutherland}` for extended-DB rows — that
row should be re-attributed to Booker–Sutherland and probably re-displayed with
`172260.bj3`.  The project's `M(2,4,4)` pilot curve remains a valid independent
realization; it is the *first-realization* credit that is wrong.

This is worth contrasting with the three groups the 2026-07-18 certificate scan
also "realized from existing banks": `[2,2,4,4]` and `[2,2,2,8]` are genuinely
absent from the DB, so those rows are safe.  Only `[2,4,4]` is affected.

## Result 3 — the RM/CM/QM breakdown (also never run)

```sql
SELECT geom_end_alg, torsion_subgroup, count(*), min(cond) FROM g2c_curves_new
WHERE is_simple_geom AND geom_end_alg <> 'Q' AND torsion_order >= 11
GROUP BY 1,2;
```

22 rows.  The two rows that matter for the current program:

```text
RM  [2,22]    1 curve   cond  19044      <- the sole [2,22], as recorded
RM  [2,2,14]  2 curves  cond 152100      <- BOTH DB [2,2,14] are RM
RM  [22]      3 curves  cond  19044      (of 37 total [22]: 34 are generic)
```

So the DB's `[2,2,14]` entries are *both* RM, which independently confirms that
the five contact-7 curves found on 2026-07-23 are the first generic witnesses.
Other RM seeds of possible interest: `[17]` (4 curves, min cond 194688), `[19]`
(2, cond 5476), `[26]` (3, cond 21600), `[13]` (3, cond 28800), and a surprising
12 RM curves with full `[2,2,2,2]`.

## Result 4 — a stale claim corrected

`notes/claude_sib_D_orthogonal.md` states that the largest Jacobian torsion
order anywhere in the LMFDB is 39.  That is **production-only**.  In the
extended DB the largest torsion order among geometrically simple curves is
**56**, attained by `[2,28]` (cond 28200) and `[2,2,14]` (cond 152100); order 48
occurs as `[2,2,2,6]` (cond 39600) and `[2,2,12]` (cond 12300).  Any argument
of the form "order 96 is far beyond anything catalogued" should be re-derived
against 56, not 39.

## What this does and does not settle

It settles *presence in the database*.  It says nothing about whether a target
is realizable — the database is a small-conductor collection, and the review
already recorded the decisive counterexample: `[2,2,14]` was "RM-only in the
DB" on 2026-07-22 and had five generic witnesses of conductor `10^13`–`10^19`
by 2026-07-23.  Small conductor is populated by Eisenstein/modular factories
that *force* RM (both RM witnesses here have square conductor, `138^2` and
`390^2`), so DB endomorphism statistics remain weak evidence about the generic
world.  Read the `[2,22]` row the same way.
