# Source-citation audit of the BookerSutherland rows (2026-08-20)

Question (Filip): Table 2 cites BookerSutherland (unpublished, bib entry 8)
as source for [3], but the production LMFDB (BSSVY 2016, entry 7) obviously
contains split [3]-examples -- isn't that incorrect?

Answer: YES for Table 2 (14 rows); NO for Table 1 (all its (8)-citations
are genuinely first).  Method: for every row citing BookerSutherland,
query the production LMFDB API for a curve with the exact torsion group
and is_simple_geom matching the table (true for Table 1, false for
Table 2).  44 queries, all resolved (logs:
$CLAUDE_JOB_DIR/tmp/lmfdb_bs_check.log, lmfdb_bs_retry.log; transcribed
here since job tmp is ephemeral).

## Table 1 (simple census): all 17 BookerSutherland citations CORRECT
[2,16] [2,18] [2,20] [2,22] [2,26] [2,28] [3,6] [4,4] [4,8] [2,2,8]
[2,2,10] [2,2,12] [2,2,14] [2,4,4] [2,2,2,2] [2,2,2,4] [2,2,2,6]:
production has NO geometrically simple example of any of them ->
the extended database is the first documented source.

## Table 2 (split census): 14 rows should cite BSSVY, not BookerSutherland
Production split witnesses found (first API match, not nec. smallest):
[3] 847.d.847.1; [4] 672.a.172032.1; [5] 847.a.847.1; [6] 448.a.448.1;
[7] 2704.b.692224.1; [8] 1152.a.147456.1; [9] 1696.b.434176.1;
[2,2] 784.c.614656.1; [2,4] 720.a.6480.1; [2,6] 600.a.96000.1;
[2,8] 600.b.30000.1; [4,4] 504.a.27216.1; [2,2,2] 1584.a.684288.1;
[2,2,4] 630.a.34020.1.
Fix: in these 14 rows replace \cite{BookerSutherland} by \cite{BSSVY}.

## Table 2 rows correctly citing BookerSutherland (21)
[14] [16] [18] [2,14] [2,16] [2,18] [2,20] [2,30] [2,48] [3,24] [4,8]
[4,12] [4,16] [2,2,10] [2,2,12] [2,2,16] [2,4,4] [2,2,2,2] [2,2,2,4]
[2,2,2,6] [2,2,2,8]: no production split example exists.

## Witness check
For all 14 corrected groups the displayed (extended-DB) witness has
strictly smaller conductor than the production example found -- witness
choices unaffected; only citations change.

## Caveats
- "First" here compares the two databases only.  The audit
  (paper/torsion_sources.tex) is by its own inclusion rule only
  REPRESENTATIVE for split constructions, so pre-2016 literature could
  still predate BSSVY for some of the 14 (e.g. classical split
  constructions with exact small torsion were never systematically
  searched).  If desired, extend the audit before publication.
- Cause of the error: the split-table source column was filled with the
  home database of the displayed witness (provenance), while the
  literature rows use first-source semantics.  The 14 corrections make
  the whole column first-source.
