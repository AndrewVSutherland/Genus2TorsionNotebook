# Lane 5: the conic-bases run was adjudicated and STOPPED by the orchestrator

Date: 2026-07-25, ~13:55 local.  **Do not spend time re-judging whether the
orphaned job was hung — that decision has been made and executed.**

## What was measured

The run was
`magma -b BASEFILE:=data/claude_ov_88conic2_bases400.txt NCH:=8 PSB:=400 NV:=6 USEMW:=1 code/claude_ov_88conic_base2.m`
with the `SetClassGroupBounds("GRH")` fix in place (line 29), so it was the
corrected computation, not a stray.  It was **not hung**: all 8 children were at
~98% CPU and writing.  The problem was the tail.  After 4 h 10 m:

```text
child   done/50   rate/h   ETA to finish   last write
part0     18       4.32       7.4 h         30 min ago
part1     14       3.36      10.7 h         29 min ago
part2     14       3.36      10.7 h         31 min ago
part3     16       3.84       8.9 h         15 min ago
part4      6       1.44      30.6 h        108 min ago   <-- inside one base
part5      6       1.44      30.6 h        101 min ago   <-- inside one base
part6     14       3.36      10.7 h         22 min ago
part7     12       2.88      13.2 h         27 min ago
```

Aggregate throughput 24 bases/h would suggest 12.5 h, but `WaitForAllChildren()`
waits for the slowest child, so the critical path was **~31 h more (~35 h
total)** — and that is a lower bound, since children 4 and 5 had each been inside
a single `MordellWeilGroup` call for over 1.5 h.  Per-base cost is heavily
tailed.

The collaborator needs claude-box back within 8 h, so the job was killed by PID
(all nine PIDs verified against `/proc/<pid>/cmdline` first).

## What you have

Every completed base is intact: the script `Flush()`es after each base, so
nothing was lost in the kill.

```text
data/claude_ov_88conic2_partial/part0..7.txt      the raw child files
results/claude_ov_88conic2_partial_harvest.log    concatenated, 108 base records
```

That is a **25% sample of the 400 bases, chosen by the sharding (idx mod 8), so
it is not biased toward easy bases** — except in one direction you must state
explicitly: bases that are cheap to compute are over-represented, because the
expensive ones are exactly the ones still unfinished.  Any rate you quote from
this sample (e.g. "fraction of fibres with a locally soluble conic") should be
reported with that caveat.

## What to do

Work from these 108 base records rather than restarting the sweep.  Lane 5's
actual question — how many fibres give a locally soluble conic, and what are the
obstruction primes where they do not — is very likely answerable at this sample
size, and the answer is GRH-conditional either way.  Aggregate the 100-part lift
profile, exact-test any lift with the validated `[8,8]`-containment gate, and
write the note.  If you conclude a larger sample is genuinely needed, say so in
your report with a per-base cost estimate and a proposed cap (e.g. a per-base
time limit, or a base subset chosen to avoid the expensive stratum) — do not
simply relaunch the same 400-base run.
