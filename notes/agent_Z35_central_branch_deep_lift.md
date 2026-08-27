# Z/35 central branch deep lift

Date: 2026-07-02.

Worker: `Z35_CENTRAL`.

This continues `agent_Z35_liftable_branch_lift.md` but only in the two
surviving central first directions of the `b=0,r=1` pole chart:

```text
t=1: <1,1,1,0,1,0>
t=2: <1,1,1,2,1,0>
```

Code:

```text
code/agent_Z35_central_branch_deep_lift.m
```

Logs used here:

```text
results/Z35_central_smoke_k4.log
results/Z35_central_deep_k8.log
results/Z35_central_deep_k10.log
```

The local convention is unchanged: after

```text
a=1+3*A, b=3*B, c0=t+3*C0, c1=t+3*C1, c2=t+3*C2, r=1+3*R,
```

the five equations are divided by one power of `3`.  Thus `H=0 mod 3^k`
means the original coefficient equations vanish modulo `3^(k+1)`.

## Next-digit affine systems

At both central directions the mod-`3` Jacobian has rank `3`, so a nonempty
next-digit set is a coset of a three-dimensional kernel.

Common obstruction rows on the five scaled residuals:

```text
[1,0,2,1,0],  [0,1,2,0,1]
```

For `t=1`, the homogeneous correction kernel is:

```text
[1,0,0,2,2,2],
[0,1,0,1,0,0],
[0,0,1,1,1,1].
```

The initial next digits satisfy:

```text
dC2 + 2*dR + 2 = 0
dB  + 2*dC1 + dR + 2 = 0
dA  + 2*dC0 + dR + 1 = 0
```

For `t=2`, the homogeneous correction kernel is:

```text
[1,0,0,1,1,2],
[0,1,0,2,0,0],
[0,0,1,1,1,2].
```

The initial next digits satisfy:

```text
dC2 + dR + 2 = 0
dB  + dC1 + dR = 0
dA  + dC0 + dR + 2 = 0
```

In all printed good-child samples, the liftable next-digit subset kept the
same three-equation affine shape, with only the constants changing.

## Exact smoke transition through `3^4`

Command:

```text
magma -b max_k:=4 direct_depth:=3 recon_height:=80 \
    max_store:=100000 sample_parent_limit:=20 \
    code/agent_Z35_central_branch_deep_lift.m \
    > results/Z35_central_smoke_k4.log
```

For each central branch:

```text
k  modulus  lifts
1  3        1
2  9        27
3  27       729
4  81       6561
```

The exact `k=3 -> 4` transition was:

```text
729 depth-3 residues:
  243 have 27 next digits
  486 have 0 next digits
```

So the total two-central-branch table through `3^4` is:

```text
k  modulus  total lifts  nonzero central branches
1  3        2            2
2  9        54           2
3  27       1458         2
4  81       13122        2
```

The smoke reconstruction probe tried three sample lifts per branch at
height `80`.  All six reconstructed exact checks returned
`reconstructed_exact_zero false`.  The samples did satisfy the original
coefficient equations modulo the recorded precision and had nonzero
discriminant.

## Interrupted `k=8` attempt

Command:

```text
magma -b max_k:=8 direct_depth:=4 recon_height:=200 \
    max_store:=100000 \
    code/agent_Z35_central_branch_deep_lift.m \
    > results/Z35_central_deep_k8.log
```

This was interrupted after the `t=1` branch reached and stored depth `4`.
Recorded rows:

```text
k  modulus  lifts
1  3        1
2  9        27
3  27       729
4  81       6561
```

No `t=2` rows or reconstruction rows were recorded in this interrupted log.

## Interrupted `k=10` diagnostic

Command:

```text
magma -b max_k:=10 direct_depth:=6 recon_height:=1000 \
    max_store:=2500000 sample_parent_limit:=1000 \
    code/agent_Z35_central_branch_deep_lift.m \
    > results/Z35_central_deep_k10.log
```

This was stopped after the `t=1` branch reached depth `6`, ran one exact
lookahead to depth `7`, and printed representative recurrence diagnostics.

Exact `t=1` direct rows:

```text
k  modulus  lifts      correction histogram
1  3        1
2  9        27         [<27,1>]
3  27       729        [<27,27>]
4  81       6561       [<0,486>, <27,243>]
5  243      177147     [<27,6561>]
6  729      1594323    [<0,118098>, <27,59049>]
```

Exact one-step lookahead:

```text
k=7, modulus 2187: 43046721 lifts, histogram [<27,1594323>].
```

The sampled `k=7 -> 8` obstruction over 1000 depth-6 parents showed an
all-or-none pattern among the 27 depth-7 children:

```text
887 sampled parents: 0 liftable children
113 sampled parents: 27 liftable children
```

The child obstruction vectors in that sample were:

```text
[0,0] : 3051
[1,2] : 11961
[2,1] : 11988
```

The representative-parent diagnostic then gave the following inferred rows
for `t=1`:

```text
k   modulus  inferred lifts  input parents  good parents  obstruction histogram
8   6561     129140163       1594323        177147        [0,0]:177147, [1,2]:708588, [2,1]:708588
9   19683    3486784401      177147         177147        [0,0]:177147
10  59049    31381059609     177147         59049         [0,0]:59049, [1,2]:118098
```

These `k=8..10` rows are not raw enumerations of all leaves.  They are
diagnostic recurrence rows based on representative descendants and the
observed stable affine coset behavior.  The run was interrupted at
`sample_lifts_and_reconstruction`, so no `k=10` rational reconstruction row
was recorded in this log.

## Verdict

The central branches are definitely not dead through scaled modulus `3^7`:
the `t=1` branch alone has `43,046,721` lifts at `3^7`, and the earlier
all-branch run recorded the same behavior for the `t=2` central branch.

Raw enumeration past this point is the wrong default.  The next obstruction is
not a smooth constant `27`-fold lift, and it is not captured safely by a naive
parity multiplier.  The useful structure is:

```text
linear rank-3 correction coset
  plus two left-obstruction residuals
  plus an all-or-none affine next-digit subset in sampled children.
```

Recommended next step: build a true compressed-state automaton keyed by the
two obstruction residuals and the changing affine constants, then prove or
disprove the all-or-none child-coset rule symbolically.  Only after that
should `3^8`, `3^9`, and `3^10` counts be treated as exact rather than
representative-recurrence diagnostics.  No rational lift candidate has appeared
in the recorded reconstruction probes.
