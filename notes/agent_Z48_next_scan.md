# Agent Z/48: next partitioned A16 scan

Date: 2026-07-02.

This is the third-pass continuation of the Z/48 lane.  I chose bounded option
A, with a small improvement from C: a partitioned larger A16 square-root scan
that explicitly skips the completed RTHeight=3 box and keeps a sharper
point-count first-kill table.

New driver:

```text
code/agent_Z48_next_scan.m
```

It is the same sign-aware A(8)->A(16) square-root system as
`agent_Z48_simultaneous_A16_plus3.m`, but with these additional controls:

```text
ExcludeRTHeight:=3       skip slices already covered by RTHeight=3
SliceMod, SliceClass     run a reproducible partition of the remaining slices
MaxSlices                optional hard cap on run slices
```

For each rational A16-equation root it applies the necessary gate

```text
48 | #J(F_l)
```

at good primes before exact Jacobian certification.  Gate survivors are then
certified exactly as A16 candidates and tested by `TorsionSubgroup` on an
integral square model.

## Smoke command

```text
magma -b RTHeight:=2 ExcludeRTHeight:=1 SearchBound:=2 \
  PrimeBound:=17 MinGood:=1 SliceMod:=2 SliceClass:=0 ProgressSlices:=1 \
  MaxExact:=2 code/agent_Z48_next_scan.m
```

Result:

```text
SEARCH_DONE rawSlices=30 completedSkipped=2 eligibleSlices=28
partitionSkipped=14 maxSliceSkipped=0 runSlices=14 sliceBuildFail=0
tested=686 commonRootPairs=3 rationalRoots=3 singular=3 nonsingular=0
pointGateReject=0 pointGatePass=0 exactTried=0 squareReject=0
d8Reject=0 d16Reject=0 certified=0 z48Hits=0
```

## Production partition

Command:

```text
magma -b RTHeight:=4 ExcludeRTHeight:=3 SearchBound:=10 \
  PrimeBound:=43 MinGood:=3 SliceMod:=4 SliceClass:=0 ProgressSlices:=10 \
  MaxExact:=50 code/agent_Z48_next_scan.m
```

This covered one quarter of the new RTHeight=4 slices, not the completed
RTHeight=3 box.

Final counters:

```text
SEARCH_DONE rawSlices=462 completedSkipped=182 eligibleSlices=280
partitionSkipped=210 maxSliceSkipped=0 runSlices=70 sliceBuildFail=0
tested=1129030 commonRootPairs=15 rationalRoots=15 singular=9
nonsingular=6 pointGateReject=5 pointGatePass=1 exactTried=1
squareReject=0 d8Reject=0 d16Reject=0 certified=1 z48Hits=0
```

Slice categories actually run:

```text
Rnew    : 28
Tnew    : 26
bothnew : 16
```

Root categories:

```text
Rnew    : 8
Tnew    : 4
bothnew : 3
```

Nonsingular root categories:

```text
Rnew    : 1
Tnew    : 2
bothnew : 3
```

Point-count first kills among the five nonsingular gate rejects:

```text
p=5 : 3
p=7 : 2
```

By category:

```text
Rnew:5    : 1
Tnew:7    : 2
bothnew:5 : 2
```

## One gate survivor

There was one new `bothnew` A16 candidate that passed the 48 point-count gate:

```text
r=-1/4, t=-1/4, mu=-1/2, y=-5/8, p=-41/144,
N=5/8, z=125/96
gate primes: <13,192>, <17,336>, <19,384>
running gcd: 48
```

Exact Jacobian certification succeeded for A16, but the integral square model
had torsion invariants

```text
[16]
```

so this is not a Z/48 hit.  The factor degrees of the defining sextic were

```text
[<1,1>, <2,1>, <3,1>]
```

The curve printed by the run was:

```text
f = -15/16*x^6 + 25/288*x^5 - 475/2304*x^4 - 125/768*x^3
    + 169375/82944*x^2 - 15625/82944*x + 390625/331776
```

## Outcome

The third-pass partition genuinely goes beyond the completed boxes:

```text
RTHeight=4 outside RTHeight=3, SliceClass 0 mod 4.
```

It produced no Z/48 hits.  Compared with the previous RTHeight=3 run, this
adds six nonsingular A16 roots outside the completed box: five are killed by
the point-count gate, and one is an exact A16 false positive for the gate with
torsion `[16]`.

The remaining reproducible partitions are:

```text
SliceClass:=1, 2, 3 with the same RTHeight=4 ExcludeRTHeight=3 parameters.
```
