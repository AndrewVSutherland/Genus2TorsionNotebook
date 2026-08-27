# Agent Z/48 cubic-contact route

Date: 2026-07-02.

Worker: `Z48_CUBIC`.

## Purpose

The RTHeight 4 partition sweep finished the larger A16 point-count route.  Its
only `48 | #J(F_l)` survivor was an exact A16 curve with torsion `[16]`, so the
next route should add the 3-primary condition itself instead of enlarging the
same A16 prefilter.

New driver:

```text
code/agent_Z48_cubic_contact_route.m
```

It reuses the sign-aware A(8)->A(16) equations from
`agent_Z48_next_scan.m`, but it sends point-count survivors through a new
3-part diagnostic before attempting to call anything a Z/48 candidate.

## Algebraic condition

For a smooth genus-2 model

```text
C: y^2 = f(x)
```

with a verified rational point of order `16`, the Z/48 question is equivalent
to:

```text
J(Q)[3] != 0.
```

Indeed an order-16 point plus any nonzero rational 3-torsion point has order
`48`, and the converse is its 3-primary component.

The explicit contact form used for the diagnostic is:

```text
h(x)^2 - f(x) = Lambda*q(x)^3,
q = x^2 + U*x + V,
h = H3*x^3 + H2*x^2 + H1*x + H0.
```

For sextic A8/A16 models the leading equation is

```text
Lambda = H3^2 - lc(f).
```

For each bounded rational `U,V`, the script solves the generic `H3 != 0`
branch recursively from the top coefficients and checks the remaining three
coefficient equations exactly.  This bounded contact search is only a witness
search.  The decisive 3-part test is exact:

```text
TorsionSubgroup(Jacobian(HyperellipticCurve(square_integral_model(f))))
```

The route therefore rejects a candidate unless the exact torsion exponent is
divisible by `3`.

## Script behavior

Modes:

```text
Mode:=smoke   tests recorded A16 examples
Mode:=search  runs bounded A16 equation generation plus the 3-part diagnostic
Mode:=both    runs both
```

Important search controls:

```text
RTHeight, ExcludeRTHeight       reuse the post-RTHeight-4 partition box
SearchBound                     height for (mu,y)
GateMode, PrimeBound, MinGood   point-count gate before exact work
ContactHeight                   height bound for U,V in the witness search
DoExact3                        exact TorsionSubgroup 3-part check
FixedRNum/FixedRDen,
FixedTNum/FixedTDen             single-slice diagnostic mode
```

The `FixedR/FixedT` controls are there to avoid rerunning a full partition
when we want to exercise the new filter on a known survivor slice.

## Smoke test

Command:

```text
magma -b Mode:=smoke ContactHeight:=2 PrimeBound:=43 \
  MinGood:=3 code/agent_Z48_cubic_contact_route.m \
  > results/Z48_cubic_contact_smoke.log 2>&1
```

Result:

```text
SMOKE_DONE candidates=7 verified16=7 rejectedByExact3=7
contactHeightHits=0 z48Unexpected=0
```

The important new negative is the SliceClass 0 RTHeight 4 gate survivor:

```text
r=-1/4, t=-1/4, mu=-1/2, y=-5/8, p=-41/144, N=5/8, z=125/96
```

It still passes the original point-count gate:

```text
good primes: <13,192>, <17,336>, <19,384>
running gcd: 48
```

but the new 3-part diagnostic gives:

```text
contactHeight=2 contactChecked=49 contactHit=false
exactOK=true has3=false torsion=[16] exponent=16
```

The six simple known A16 examples from the earlier notes also verify as A16
and are rejected by exact 3-part torsion; every one has exact torsion `[16]`.

## Bounded search runs

First tiny RTHeight 4 outside-RTHeight 3 sample:

```text
magma -b Mode:=search RTHeight:=4 ExcludeRTHeight:=3 \
  SearchBound:=2 PrimeBound:=19 MinGood:=1 GateMode:=48 MaxSlices:=12 \
  MaxA16Roots:=20 MaxExact:=5 ContactHeight:=1 ProgressSlices:=2 \
  code/agent_Z48_cubic_contact_route.m \
  > results/Z48_cubic_contact_small_search.log 2>&1
```

Counts:

```text
runSlices=12 tested=588 commonRootPairs=0 rationalRoots=0
nonsingular=0 pointGatePass=0 threeChecked=0 z48Hits=0
```

Fixed survivor slice, just large enough to include `mu=-1/2, y=-5/8`:

```text
magma -b Mode:=search FixedRNum:=-1 FixedRDen:=4 \
  FixedTNum:=-1 FixedTDen:=4 SearchBound:=8 PrimeBound:=43 MinGood:=3 \
  GateMode:=48 MaxSlices:=1 MaxA16Roots:=20 MaxExact:=5 ContactHeight:=1 \
  ProgressSlices:=1 code/agent_Z48_cubic_contact_route.m \
  > results/Z48_cubic_contact_fixed_slice.log 2>&1
```

Counts:

```text
tested=7569 commonRootPairs=1 rationalRoots=1 singular=0 nonsingular=1
pointGatePass=1 threeChecked=1 threeReject=1 threePass=0
contactHits=0 exactTried=0 certified=0 z48Hits=0
```

This is the new route doing the intended thing: the old false positive reaches
the `48 | #J(F_l)` gate, but is killed by exact absence of rational 3-torsion.

## Verdict

The cubic-contact/3-part layer is now wired into the A16 candidate stream.  It
rejects the known `[16]` false positive and all simple A16 examples.  No Z/48
candidate appeared in the bounded smoke or fixed-slice runs.

Next production parameters should stay fixed-slice or partitioned, not blind:

```text
1. Use GateMode:=48, PrimeBound:=43, MinGood:=3 as before.
2. Keep ContactHeight:=1 or 2 for witness diagnostics.
3. Keep DoExact3:=true; it is the decisive test.
4. Run one RTHeight 4 partition at a time only if there is a reason to enlarge
   SearchBound beyond the completed A16 sweep, e.g. SearchBound:=12 or 14.
5. Prefer fixed-slice reruns for any future point-count survivors so the exact
   3-part filter is exercised immediately.
```
