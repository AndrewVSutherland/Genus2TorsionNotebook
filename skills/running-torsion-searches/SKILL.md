---
name: running-torsion-searches
description: Design, launch, monitor, and know when to KILL a Magma parameter search - the staged funnel architecture, NParts sharding, greppable log markers (PROGRESS/HIT/SEARCH_DONE), background launch + watcher patterns, height-enumeration subtleties, rate extrapolation, pipeline validation on a known hit, and prefilter false-positive budgeting. WHEN starting any torsion search, when a run seems stuck or silent, when deciding scan height, or when interpreting a search that found nothing. Trigger words - search funnel, height scan, NParts, Part, background job, monitor, PROGRESS, SEARCH_DONE, kill a scan, smoke test.
---

# Running torsion searches

## When to use this

Load this before launching any parameter search, and again before *believing*
one. A search here is a staged funnel over a parametrized family, sharded into
background Magma jobs, monitored through greppable log markers, and — crucially
— validated on a known positive before its negatives mean anything.

## 1. The funnel architecture

Stage order, strictly by increasing unit cost. Verbatim structure from the
production `[2,24]` search `code/agent_a2_24_wsplit_3tors.m` (its main loop):

```text
enumerate parameters (r, t, beta)                      free
  -> parametrized condition (W-split p-roots)          one Roots() call
  -> degenerate skips:                                  free
       Degree(f) lt 5, Discriminant(f) eq 0,
       even-sextic skip: Coefficient(f,1) eq 0 and
         Coefficient(f,3) eq 0 and Coefficient(f,5) eq 0   // bielliptic => split => worthless
  -> IntModel (clear denominators)                      free
  -> algebraic filter: TwoRank(fInt) ge 2               one Factorization
  -> finite prefilter: 3|#J(F_q), q in PreP             ~14 point counts
  -> marked-class exact order: D8 (try/catch)           a few Jacobian adds
  -> exact TorsionSubgroup(J)                           EXPENSIVE
  -> SimplicityCertificate                              only on actual hits
```

**Instrument every stage with a counter** (`tested`, `rank2`, `pre3`, `tors`,
`hits24`, ...) and print them all in the `PROGRESS` lines — the funnel's shape
in the log is the primary debugging and design tool. See `finite-prefilters`
for the middle stages, `two-rank-and-factor-types` for the algebraic filter,
`simplicity-certificates` for the last.

## 2. Sharding

Shard on the **outermost** loop only, with a running index (verbatim pattern,
all search scripts):

```magma
ridx := 0;
for rv in rvals do
    ridx +:= 1;
    if (ridx mod NParts) ne Part then continue; end if;
    ...
```

Launch `NParts` (house standard: 3, matching the concurrent-job cap) parts as
separate background jobs, each writing its own log
`data/<name>_part<k>.log`. Parts are independent; a crashed part loses only
its own slice.

## 3. Log marker conventions (machine-greppable)

One line per event, fixed leading token, **full parameters plus the integral
`f` on hit lines** so any hit is reproducible from the log alone:

```text
PROGRESS tested=... p=... rank2=... pre3=... tors=... hits24=... hits224=...
HIT24 r=5 p=-5/2 t=-9/2 torsion=[ 24 ] 2rank=1 simple=true (q=17 chi=...)
  f=701706240*x^6 - ...
TARGET_2_24 r=... t=... b=... p=... torsion=... simple=...
SEARCH_DONE tested=... p=... rank2=... pre3=... tors=... hits24=... hits224=...
TORSION_HISTOGRAM
  [ 2, 8 ] : 1708
DONE
```

(Formats from `code/agent_a2_24_wsplit_3tors.m` /
`code/agent_a2_24_ztors_sample.m`; the recorded `HIT24` line above is real —
verbatim in `data/agent_a2_24_composite_h12_part1.log`, and recorded in prose
in `notes/agent_a2_24_composite.md`.) `PROGRESS` should fire every
`progress` iterations (a passed parameter, typically 1e5–3e6 depending on
stage cost). The summary of a finished run is then a one-liner:
`grep -h "SEARCH_DONE" data/<name>_part*.log` and
`grep -h "^HIT" data/<name>_part*.log`. A `TORSION_HISTOGRAM` block at the end
of each part turns even a zero-hit run into data (it is how "the prefilter
survivors are all `[2,8]`" was learned).

## 4. Launch and monitor

Launch (respecting the ≤3-job cap; see `magma-lab-conventions` §4–5 for `-b`
buffering, `SetColumns(0)`, `quit;`, `MemGB`):

```bash
for k in 0 1 2; do
  nohup magma -b H:=26 NParts:=3 Part:=$k MemGB:=3 \
    code/agent_a2_24_wsplit_3tors.m > data/wsplit_h26_part${k}.log 2>&1 &
done
```

Watcher loop — a monitor must match **failure signatures too**, because a
crashed job's silence is indistinguishable from "still running":

```bash
until grep -qE "SEARCH_DONE" data/wsplit_h26_part*.log \
   || grep -qE "User error|Runtime error|Segmentation|Killed|total memory" data/wsplit_h26_part*.log \
   || ! pgrep -f wsplit_3tors >/dev/null; do sleep 60; done
```

Check liveness with `pgrep -f <script-name> | wc -l`. Never poll by killing
the job to look at stdout (`-b` buffering: you would see nothing and lose the
run).

## 5. Height enumeration subtleties

The standard enumerator (inline in each script; also
`TC_HeightRationals` in `code/torsion_cover_lab_utils.m`):

```magma
function HeightRationals(HH)
    vals := [];
    for den in [1..HH] do for num in [-HH..HH] do
        if GCD(num, den) eq 1 then Append(~vals, Q!num/den); end if;
    end for; end for;
    return Sort(Setseq(Seqset(vals)));
end function;
```

**It returns values sorted BY VALUE, not by height.** Consequence: a specific
low-height point of interest can sit deep in scan order — the known `[2,24]`-
adjacent example at `r = 1/3` is index ~7104 of 12175 at `H=100`, so a scan
runs for hours before touching the one parameter you care about. If early
low-height coverage matters (it usually does), either iterate by increasing
height, or smoke-test the known point explicitly first (§7). Size grows
roughly like `H^2` — `#HeightRationals(45) = 2511` values (whence the
`2511^2 ≈ 6.3M` pairs of the H=45 `d=0` scan), and a 3-parameter scan cubes
it.

## 6. Rate measurement and extrapolation

Before any big run:

1. **Syntax smoke test:** `timeout 40 magma -b RH:=2 TH:=2 BH:=2 code/x.m` —
   a full mini-run to `DONE` in seconds catches syntax errors and
   return-arity bugs that would otherwise surface hours in.
2. **Rate test:** run one part a few minutes, read the first `PROGRESS`
   lines, compute tested/sec, extrapolate the full grid.
3. **Apply the high-height correction:** rational arithmetic slows several-
   fold as numerators/denominators grow — observed live during the
   `agent_a2_24_d0_fast.m` H=100 run: a scan that opened at ~25k pairs/sec
   dropped to ~5k/sec deep into the height range (session observation; the
   run was killed and its log not preserved, so treat the exact ratio as a
   rule of thumb and re-measure on your own run). Early rates overestimate;
   budget accordingly.

A search that cannot finish in hours should be redesigned — a stronger
parametrized condition (the W-split move, `two-rank-and-factor-types`) or a
smaller decisive slice — not left running for days.

## 7. Validate the pipeline on a known hit

**A funnel that has never seen a positive is untested.** Before trusting any
negative result, plant a known example inside the scan range and confirm the
funnel flags it end-to-end (all stages, including the hit-line printing). The
model: the `d=0` fast scan (`code/agent_a2_24_d0_fast.m`) at `H=45` returned
exactly one genuine hit — the already-known curve B at `(r,t)=(1/3,-1)` — which
simultaneously validated the machinery and calibrated the density
(`notes/agent_a2_24_d0_derivation.md`). A subtler failure this catches: the
Mumford scaling mistake (`magma-lab-conventions` §6) makes the marked-class
stage silently reject *everything*, which in an unvalidated funnel reads as
"the locus is empty".

## 8. Prefilter false-positive budgeting

From the production `[2,24]` run (`notes/agent_a2_24_composite.md`): ~1M
curves passed the 2-rank stage, 1709 passed the 14-prime 3-torsion prefilter,
**zero** were genuine (`TORSION_HISTOGRAM [2,8]:1708`). Survivor counts times
`TorsionSubgroup` cost is a real budget line — see `finite-prefilters` for the
calibration reasoning. Design the funnel so the exact stage's expected load is
affordable *before* launching.

## 9. When to KILL a scan

Kill it when:

- **The answer is already decided.** A density/rationality question can be
  settled by a complete low-height scan: 1 genuine point in 6.3M pairs
  settles "this cover is not a dense rational family" — more height adds
  nothing (`notes/agent_a2_24_d0_derivation.md`). State the decision rule
  before launching ("if < N hits by height H, conclude X").
- **Measured rates show it cannot reach the decisive region** in acceptable
  time (§6) — redesign instead.
- **A part crashed** (failure signature in the log, or `pgrep` count
  dropped): fix and relaunch that part; the shards are independent.

When you kill, **record the decision and the counts** in the route's note —
a killed scan with recorded counters is data; a killed scan without is wasted
compute. Do not kill expecting partial stdout (buffering).

## Pitfalls

- **Stages out of cost order** — e.g. exact order checks before the finite
  prefilter. The funnel's whole value is rejecting cheaply first.
- **No stage counters** — a zero-hit run with no counters cannot distinguish
  "filter too strong" from "family empty" from "bug".
- **Sharding an inner loop** — parts then overlap or skip work; shard the
  outermost index only.
- **Monitors that only match success** — add the failure alternation
  (`User error|Runtime error|Segmentation|Killed|total memory`); silence is
  not progress.
- **Believing a negative from an unvalidated funnel** (§7) — plant the known
  positive first.
- **Assuming value-sorted = height-sorted** (§5) — your target parameter may
  be 60% of the way through the scan.
- **Leaving a decided scan running** (§9) — the marginal information of more
  height is often exactly zero.
- **Over-cap concurrency** — 3 jobs max, `MemGB` set, or the machine OOMs
  and takes all parts down (`magma-lab-conventions` §5).

## See also

- `magma-lab-conventions` — batch mode, parameter passing, memory, scaling
  rule; read it first.
- `finite-prefilters` — the middle funnel stages and their calibration.
- `two-rank-and-factor-types` — the free algebraic stage; the W-split move
  that replaces filtering with parametrizing.
- `simplicity-certificates`, `validate-and-record-a-hit` — what happens to a
  hit after the funnel.
- `local-obstructions`, `component-boundary-analysis` — what to do when a
  well-validated funnel finds nothing.
- `g2-torsion-lab` — hub.
