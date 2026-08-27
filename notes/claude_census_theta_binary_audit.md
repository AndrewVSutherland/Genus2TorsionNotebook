# Binary Theta column audit (2026-08-19)

Task: replace the census Theta column (witness-level, three states) by the
binary all-examples version: `+` iff EVERY known example of the group has a
torsion class of order equal to the largest invariant factor of
doubled-point type 2P-K; `--` otherwise.  Deliverable:
`paper/census_table_theta_binary.tex`.

## Logic and universe

- Any row with witness-level circ/dash verdict is `--` immediately (the
  witness itself is a counterexample to "every").
- The 17 bullet rows were audited against: all LMFDB curves with the exact
  torsion group and is_simple_geom = true (www + alpha, deduped by label;
  data/claude_census_2pk_lmfdb_pool.m, 165 curves), all repo-recorded
  family fibers (contact-7, contact-9, Flynn-11, Daowsud-Schmidt-11,
  M(12)-z, Leprevost-21, Kuru-Sadek-23), the displayed witnesses
  (data/claude_census_2pk_audit.txt), and BFT 2014 Example 19 (with the
  13x^3-105 misprint correction; the printed 12 gives only [3]).
  Detector identical to code/claude_census_2pk_audit.m.
  Script: code/claude_census_2pk_all_audit.m ->
  results/claude_census_2pk_all_audit.log (193 audits).

## Verdicts

`+` (12 rows): [13] [15] [17] [19] [21] [23] [24] [25] [27] [29] [33] [39].

Flips bullet -> `--` (5 rows, with counterexamples):
- [7]: 14/90 fail, e.g. 4815.a.14445.1, 4925.a.4925.1, 10468.a.83744.1.
- [9]: 4/33 fail, e.g. 2768.a.354304.1, 3568.d.913408.1.
- [11]: 1/16 fails: 4489.c.4489.1 (conductor 67^2).
- [12]: 18/22 fail, incl. ALL four M(12)-z family fibers; the bullet
  witness 762.a.3048.1 is one of only 4 passing examples.
- [3,3]: BFT Ex. 19 y^2 = (13x^3-105)^2 - 12(x^2-3x-3)^3 has exact [3,3],
  no maximal-order doubled class (witness 5100.a passes).

Structural remark: for the odd-order contact/CF families the maximal class
g = [R - infty] has 2g = [2R - K] doubled, so every fiber passes
automatically; the audited fibers confirm.  No such protection exists for
even orders -- hence the [12] collapse.

## Caveats

- Not audited (no equations in repo): the LSSV 2024 PQM [3,3]-curve
  (moot: [3,3] is already `--` via BFT Ex. 19), members of Leprevost's
  [13]/[17]/[19] families beyond LMFDB, Platonov-Petrunin examples beyond
  LMFDB.  A future counterexample among these could only flip a current
  `+` to `--`, never the reverse.
- In passing, the [20] row of the v3 table has a missing `&` before
  $\infty$ (cells shift by one); fixed in the delivered file, flagged to
  Filip.
- The caption's Theta legend line was updated to the binary meaning; the
  prose audit paragraph above the table still describes the old
  bullet/circ convention and needs a matching one-sentence edit (proposed
  text given in the session reply).
