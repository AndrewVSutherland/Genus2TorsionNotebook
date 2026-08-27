# Agent Z/48: next scan SliceClass 3

Date: 2026-07-02.

This is the fourth reproducible partition of the RTHeight 4 boundary outside
the completed RTHeight 3 box for the partitioned sign-aware A(8)->A(16)
square-root scan.

Command:

```text
magma -b RTHeight:=4 ExcludeRTHeight:=3 SearchBound:=10 \
  PrimeBound:=43 MinGood:=3 SliceMod:=4 SliceClass:=3 ProgressSlices:=10 \
  MaxExact:=50 code/agent_Z48_next_scan.m
```

Raw log:

```text
results/Z48_next_slice3.log
```

Final counters:

```text
SEARCH_DONE rawSlices=462 completedSkipped=182 eligibleSlices=280
partitionSkipped=210 maxSliceSkipped=0 runSlices=70 sliceBuildFail=0
tested=1129030 commonRootPairs=5 rationalRoots=5 singular=2
nonsingular=3 pointGateReject=3 pointGatePass=0 exactTried=0
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
Rnew : 5
```

Nonsingular root categories:

```text
Rnew : 3
```

Point-count first kills among the three nonsingular gate rejects:

```text
p=5 : 1
p=7 : 2
```

By category:

```text
Rnew:5 : 1
Rnew:7 : 2
```

There were no `GATE_PASS_A16_PLUS3` lines in this partition.  Consequently no
exact A16 certification or `TorsionSubgroup` computation was attempted:

```text
pointGatePass=0
exactTried=0
certified=0
CERTIFIED_TORSION_HISTOGRAM is empty
```

Outcome:

```text
RTHeight=4 outside RTHeight=3, SliceClass 3 mod 4.
```

This partition produced no Z/48 hit.  Its only nonsingular A16-equation roots
were three `Rnew` roots, all killed by the finite-field point-count gate before
exact certification.
