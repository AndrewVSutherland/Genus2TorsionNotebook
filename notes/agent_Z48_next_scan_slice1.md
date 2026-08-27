# Agent Z/48 next scan, SliceClass 1

Date: 2026-07-02.

This is the second production partition of the RTHeight=4 search outside the
completed RTHeight=3 box, using the partitioned sign-aware A(8)->A(16) driver.

Command:

```text
magma -b RTHeight:=4 ExcludeRTHeight:=3 SearchBound:=10 \
  PrimeBound:=43 MinGood:=3 SliceMod:=4 SliceClass:=1 ProgressSlices:=10 \
  MaxExact:=50 code/agent_Z48_next_scan.m
```

Output was captured in:

```text
results/Z48_next_slice1.log
```

The command completed successfully.

Final counters:

```text
SEARCH_DONE rawSlices=462 completedSkipped=182 eligibleSlices=280
partitionSkipped=210 maxSliceSkipped=0 runSlices=70 sliceBuildFail=0
tested=1129030 commonRootPairs=4 rationalRoots=4 singular=4
nonsingular=0 pointGateReject=0 pointGatePass=0 exactTried=0
squareReject=0 d8Reject=0 d16Reject=0 certified=0 z48Hits=0
```

Slice categories actually run:

```text
Rnew    : 28
Tnew    : 26
bothnew : 16
```

Root categories:

```text
Tnew    : 3
bothnew : 1
```

There were no nonsingular roots in this partition.  Therefore the point-count
gate was never reached:

```text
FIRST_POINT_GATE_KILLS: empty
FIRST_POINT_GATE_KILLS_BY_CATEGORY: empty
```

There were no gate survivors, no exact Jacobian certifications, and no torsion
invariant computations:

```text
CERTIFIED_TORSION_HISTOGRAM: empty
```

Outcome: no Z/48 hit occurred in `SliceClass:=1`.
