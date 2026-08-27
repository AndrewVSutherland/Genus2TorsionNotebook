# Agent Z/48 next scan: SliceClass 2

Date: 2026-07-02.

Worker: `Z48_PART2`.

This records the requested reproducible partition of the Z/48 A(16)+3
point-count gate scan:

```text
magma -b RTHeight:=4 ExcludeRTHeight:=3 SearchBound:=10 \
  PrimeBound:=43 MinGood:=3 SliceMod:=4 SliceClass:=2 ProgressSlices:=10 \
  MaxExact:=50 code/agent_Z48_next_scan.m
```

The raw terminal transcript is saved in:

```text
results/Z48_next_slice2.log
```

## Final counters

```text
SEARCH_DONE rawSlices=462 completedSkipped=182 eligibleSlices=280
partitionSkipped=210 maxSliceSkipped=0 runSlices=70 sliceBuildFail=0
tested=1129030 commonRootPairs=8 rationalRoots=10 singular=10
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
Rnew : 10
```

There were no nonsingular root categories.

## Point-count gate

All ten rational roots in this partition were singular.  Therefore no
nonsingular candidate reached the point-count gate.

First point-count kills:

```text
none
```

First point-count kills by category:

```text
none
```

## Gate survivors and exact torsion

There were no gate survivors:

```text
pointGatePass=0
exactTried=0
certified=0
```

Consequently, no exact torsion invariant was computed in this partition.  The
certified torsion histogram was empty.

## Outcome

This completes the `RTHeight=4` outside `RTHeight=3` partition with
`SliceMod:=4, SliceClass:=2`.

No Z/48 hit occurred:

```text
z48Hits=0
```
