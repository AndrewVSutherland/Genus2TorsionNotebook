# Z/35 liftable b=0 pole branch lift

Date: 2026-07-02.

This continues `agent_Z35_b0_pole_blowup.md` and pushes the eighteen
third-pass liftable first directions beyond the mod `27` table.

Code:

```text
code/agent_Z35_liftable_branch_lift.m
```

Run:

```text
magma -b max_k:=7 recon_height:=80 \
    code/agent_Z35_liftable_branch_lift.m \
    > results/Z35_lift_branch_k7.log
```

The local convention is the same as in the previous note: after

```text
a=1+3*A, b=3*B, c0=t+3*C0, c1=t+3*C1, c2=t+3*C2, r=1+3*R,
```

the transformed equations are divided by one power of `3`; the table below is
for the scaled equations `H=0 mod 3^k`.  Thus the original coefficient
equations vanish modulo `3^(k+1)`.

## First directions

The nine liftable first directions for each center form affine planes over
`F_3`.

For `t=1`:

```text
B=2-A, R=C0-A, C1=C0+A+1, C2=C0+1-A.
```

For `t=2`:

```text
B=2-A, R=2-C0-A, C1=C0+2-A, C2=C0+A+2.
```

After blowing up an individual first direction once more, the next-digit
systems all have `27` mod-`3` solutions with rank histogram

```text
[ <3, 27> ].
```

This rank-3 first impression is not enough: higher Hensel steps have real
obstructions.

## Exact lift table

Across all eighteen first directions:

```text
k   modulus  total lifts  nonzero branches
1   3        18           18
2   9        486          18
3   27       13122        18
4   81       39366        6
5   243      1062882      6
6   729      15943230     6
7   2187     86093442     2
```

Per center, the counts are exactly half of these.

Branch classification:

```text
t=1 dies at mod 81:
  <0,2,0,1,1,0>
  <0,2,2,0,0,2>
  <1,1,0,2,0,2>
  <1,1,2,1,2,1>
  <2,0,0,0,2,1>
  <2,0,2,2,1,0>

t=1 lifts to mod 729 but dies at mod 2187:
  <0,2,1,2,2,1>
  <2,0,1,1,0,2>

t=1 survives through mod 2187:
  <1,1,1,0,1,0>

t=2 dies at mod 81:
  <0,2,0,2,2,2>
  <0,2,2,1,1,0>
  <1,1,0,1,0,1>
  <1,1,2,0,2,2>
  <2,0,0,0,1,0>
  <2,0,2,2,0,1>

t=2 lifts to mod 729 but dies at mod 2187:
  <0,2,1,0,0,1>
  <2,0,1,1,2,2>

t=2 survives through mod 2187:
  <1,1,1,2,1,0>
```

The two surviving directions are the central directions in the two affine
planes.

## Smoothness and reconstruction

For each branch the script stores a sample at the deepest retained precision
and verifies the original coefficient equations at that precision.  The sample
curves have nonzero discriminant; the observed `3`-adic discriminant valuations
are positive, as expected in this bad `p=3` chart, but no sampled branch is
forced singular by these checks.

The point-contact-5 equations are exactly the equations being lifted, so any
exact smooth rational lift would carry the necessary 5-contact.  A small
height-80 rational reconstruction probe was run on the stored samples.  Every
reported `reconstructed_exact_zero` value is `false`; no exact rational point
or clean rational pattern appeared.

## Verdict

The earlier mod `27` table was overly optimistic.  The `b=0,r=1` pole chart is
not uniformly a smooth 3-dimensional 3-adic family:

```text
12/18 branches die immediately at mod 81,
4/18 more die at mod 2187,
2/18 survive through mod 2187.
```

So this chart is not dead, but it has narrowed sharply.  The next useful Z/35
target is to isolate the two surviving central branches and push their
obstruction equations/parameterization beyond `3^7`, rather than spending more
time on the twelve mod-`81` dead directions.
