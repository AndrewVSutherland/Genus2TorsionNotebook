# Agent Z/48 cubic production

Date: 2026-07-02.

Worker: `Z48_CUBIC_PRODUCTION`.

## Purpose

This was a bounded next production scout for the Z/48 lane using the
cubic-contact/exact-3 layer, not the older A16-only gate.  I did not rerun the
completed `RTHeight:=4`, `SearchBound:=10` partition sweep.  The new partition
scout moved outward to the shell

```text
RTHeight:=5 ExcludeRTHeight:=4 SearchBound:=12.
```

I also ran one targeted fixed-slice calibration at the known old point-count
survivor, with `SearchBound:=12`, to exercise the exact 3-part filter on a
slice where a point-count survivor is known to occur.

## Driver

The original route file remains:

```text
code/agent_Z48_cubic_contact_route.m
```

In this Magma invocation, non-fixed mode in that file hit an unassigned
optional `FixedRDen`/`FixedR` initialization path before `RunSearch` was
available.  I therefore made the allowed wrapper copy:

```text
code/agent_Z48_cubic_production_driver.m
```

It is copied from `agent_Z48_cubic_contact_route.m`; the only functional
change is robust initialization of optional fixed-slice parameters.  The
A16 equations, point-count gate, cubic-contact diagnostic, and exact
`TorsionSubgroup` 3-part test are unchanged.

## Commands

New RTHeight 5 shell, four selected residue classes, each capped at 10 slices:

```text
magma -b Mode:=search RTHeight:=5 ExcludeRTHeight:=4 \
  SearchBound:=12 PrimeBound:=43 MinGood:=3 GateMode:=48 SliceMod:=4 \
  SliceClass:=0 MaxSlices:=10 MaxA16Roots:=30 MaxExact:=10 \
  ContactHeight:=2 ProgressSlices:=1 MaxGatePrint:=20 \
  code/agent_Z48_cubic_production_driver.m \
  > results/Z48_cubic_prod_rt5_shell_s12_c0.log 2>&1

magma -b Mode:=search RTHeight:=5 ExcludeRTHeight:=4 \
  SearchBound:=12 PrimeBound:=43 MinGood:=3 GateMode:=48 SliceMod:=4 \
  SliceClass:=1 MaxSlices:=10 MaxA16Roots:=30 MaxExact:=10 \
  ContactHeight:=2 ProgressSlices:=1 MaxGatePrint:=20 \
  code/agent_Z48_cubic_production_driver.m \
  > results/Z48_cubic_prod_rt5_shell_s12_c1.log 2>&1

magma -b Mode:=search RTHeight:=5 ExcludeRTHeight:=4 \
  SearchBound:=12 PrimeBound:=43 MinGood:=3 GateMode:=48 SliceMod:=4 \
  SliceClass:=2 MaxSlices:=10 MaxA16Roots:=30 MaxExact:=10 \
  ContactHeight:=2 ProgressSlices:=1 MaxGatePrint:=20 \
  code/agent_Z48_cubic_production_driver.m \
  > results/Z48_cubic_prod_rt5_shell_s12_c2.log 2>&1

magma -b Mode:=search RTHeight:=5 ExcludeRTHeight:=4 \
  SearchBound:=12 PrimeBound:=43 MinGood:=3 GateMode:=48 SliceMod:=4 \
  SliceClass:=3 MaxSlices:=10 MaxA16Roots:=30 MaxExact:=10 \
  ContactHeight:=2 ProgressSlices:=1 MaxGatePrint:=20 \
  code/agent_Z48_cubic_production_driver.m \
  > results/Z48_cubic_prod_rt5_shell_s12_c3.log 2>&1
```

Targeted fixed-slice calibration through the original route:

```text
magma -b Mode:=search FixedRNum:=-1 FixedRDen:=4 \
  FixedTNum:=-1 FixedTDen:=4 SearchBound:=12 PrimeBound:=43 MinGood:=3 \
  GateMode:=48 MaxSlices:=1 MaxA16Roots:=30 MaxExact:=10 ContactHeight:=2 \
  ProgressSlices:=1 MaxGatePrint:=20 \
  code/agent_Z48_cubic_contact_route.m \
  > results/Z48_cubic_prod_fixed_survivor_s12.log 2>&1
```

## RTHeight 5 Shell Counts

Aggregate over the four RTHeight 5 shell logs:

```text
runSlices=40
tested=1339560
commonRootPairs=4
rationalRoots=4
singular=3
nonsingular=1
pointGateReject=1
pointGatePass=0
threeChecked=0
threeReject=0
threePass=0
contactHits=0
exactTried=0
certified=0
z48Hits=0
```

The only nonsingular root in this shell sample was a `bothnew` slice and was
killed immediately by the point-count gate:

```text
FIRST_POINT_GATE_KILLS
  5 : 1
FIRST_POINT_GATE_KILLS_BY_CATEGORY
  bothnew:5 : 1
```

Slice categories actually sampled:

```text
Rnew    : 22
Tnew    : 2
bothnew : 16
```

Root categories:

```text
Rnew    : 1
bothnew : 3
```

Nonsingular categories:

```text
bothnew : 1
```

## Fixed Survivor Calibration

The targeted fixed slice found the known point-count survivor:

```text
r=-1/4, t=-1/4, mu=-1/2, y=-5/8, p=-41/144, N=5/8, z=125/96
gate used=[ <13, 192, 0, 192>, <17, 336, 0, 48>, <19, 384, 0, 48> ]
```

The new 3-part layer rejected it exactly:

```text
contactHeight=2
contactChecked=49
contactHit=false
exactOK=true
has3=false
torsion=[ 16 ]
order=16
exponent=16
```

Fixed-slice counters:

```text
tested=33489
rationalRoots=1
nonsingular=1
pointGatePass=1
threeChecked=1
threeReject=1
threePass=0
contactHits=0
certified=0
z48Hits=0
```

Since `threePass=0`, no exact A16-with-3 certificate was triggered and no
curve was printed as a Z/48 candidate.

## Verdict

The new RTHeight 5 shell scout is cold in this bounded sample: one nonsingular
A16-equation root appeared and was killed at `p=5` before the exact 3-part
stage.  The targeted point-count survivor still has exact torsion `[16]`, so
it has no rational 3-primary component and is not a Z/48 curve.

No cubic-contact witness, exact 3-part pass, certified A16-with-3 curve, or
Z/48 hit occurred in these production runs.
