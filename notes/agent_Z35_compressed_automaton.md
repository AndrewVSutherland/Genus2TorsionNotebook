# Z/35 compressed obstruction automaton

Date: 2026-07-02.

Worker: `Z35_COMPRESSED_AUTOMATON`.

Code:

```text
code/agent_Z35_compressed_automaton.m
```

Main logs:

```text
results/Z35_compressed_smoke_k4.log
results/Z35_compressed_markov_smoke_k4.log
results/Z35_compressed_k10.log
```

The local convention is unchanged from the previous Z/35 notes.  After

```text
a=1+3*A, b=3*B, c0=t+3*C0, c1=t+3*C1, c2=t+3*C2, r=1+3*R,
```

the five equations are divided by one power of `3`.  Thus `H=0 mod 3^k`
means the original coefficient equations vanish modulo `3^(k+1)`.

## State

The visible obstruction state is

```text
S_k(x) = (t, c1,c2,c3, o1,o2) in {1,2} x F_3^3 x F_3^2.
```

Here `c=(c1,c2,c3)` are the constants in the affine next-digit system

```text
E_i(d) + c_i = 0,  i=1,2,3,
```

and `o=(o1,o2)` are the two left-obstruction residuals.  The common
left-obstruction rows on the five scaled residuals are

```text
[1,0,2,1,0],  [0,1,2,0,1].
```

For `t=1`, the next-digit equations use

```text
E1 = dC2 + 2*dR
E2 = dB + 2*dC1 + dR
E3 = dA + 2*dC0 + dR
```

with initial state

```text
<t=1, c=[2,2,1], obs=[0,0]>.
```

For `t=2`, the equations use

```text
E1 = dC2 + dR
E2 = dB + dC1 + dR
E3 = dA + dC0 + dR
```

with initial state

```text
<t=2, c=[2,0,2], obs=[0,0]>.
```

When `obs=[0,0]`, the next digits are evaluated as

```text
d0 + span(kernel_rows)
```

with the branch-specific three-row kernels recorded in the log.  The script
does not enumerate all raw descendants for the certified `k=6 -> k=8`
lookahead; it groups parents by the exact mod-`9` tail

```text
(x mod 9, H(x)/3^k mod 9).
```

## Commands

Smoke validation through scaled `3^4`:

```text
magma -b max_k:=4 seed_depth:=2 print_conflict_limit:=4 \
    code/agent_Z35_compressed_automaton.m \
    > results/Z35_compressed_smoke_k4.log
```

Projected-state Markov test:

```text
magma -b max_k:=4 seed_depth:=2 use_finite_tail:=1 \
    print_conflict_limit:=4 \
    code/agent_Z35_compressed_automaton.m \
    > results/Z35_compressed_markov_smoke_k4.log
```

Main bounded run:

```text
magma -b max_k:=10 print_conflict_limit:=4 \
    code/agent_Z35_compressed_automaton.m \
    > results/Z35_compressed_k10.log
```

## Exact validation

The smoke run matches the exact direct `3^4` table for the two central
branches:

```text
k  modulus  total lifts
1  3        2
2  9        54
3  27       1458
4  81       13122
```

For `t=1`, the main run also matches the recorded exact rows through scaled
`3^7`:

```text
k  modulus  lifts
1  3        1
2  9        27
3  27       729
4  81       6561
5  243      177147
6  729      1594323
7  2187     43046721
```

The same direct seed and certified `k=7` row were obtained for `t=2`.

## Certified counts

For each branch, the `k=6` parent table has `1,594,323` residues.  Grouping
by the exact mod-`9` tail gives only `13,122` state classes for the certified
two-step lookahead.

Per central branch:

```text
k  modulus  lifts       status
1  3        1           certified direct seed
2  9        27          certified direct seed
3  27       729         certified direct seed
4  81       6561        certified direct seed
5  243      177147      certified direct seed
6  729      1594323     certified direct seed
7  2187     43046721    certified grouped mod-9 cosets
8  6561     129140163   certified grouped mod-9 child obstruction
```

For the two central branches together:

```text
k  modulus  total lifts
1  3        2
2  9        54
3  27       1458
4  81       13122
5  243      354294
6  729      3188646
7  2187     86093442
8  6561     258280326
```

The obstruction code in the logs is `o1 + 3*o2`, so code `0` is `[0,0]`,
code `5` is `[2,1]`, and code `7` is `[1,2]`.

## Representative diagnostics

The previous representative recurrence was reproduced for both central
branches, but these rows remain non-certified:

Per branch:

```text
k   modulus  inferred lifts  status
9   19683    3486784401      representative_not_certified
10  59049    31381059609     representative_not_certified
```

Two-branch totals:

```text
k   modulus  inferred total
9   19683    6973568802
10  59049    62762119218
```

## All-or-none and Markov verdict

No counterexample to the all-or-none child-coset rule was found in the
certified scans: every encountered liftable parent coset had either `0` or
`27` children with zero next obstruction.  In the main `k=10` run the
certified `k=6 -> k=8` child-good histogram per branch was

```text
0 good children:  1417176 parents
27 good children: 177147 parents
```

However, the projected five-symbol state `(t,c,obs)` is **not** itself a
Markov state.  The finite-tail smoke run found transition conflicts already
at `k=2`; for example, the same visible state

```text
<t=1, c=[2,2,1], obs=[0,0]>
```

has exact tail representatives with different child-good counts (`27` versus
`0`).  Thus the small symbolic state is useful for reporting obstructions,
but certified propagation needs at least the finite tail used by the mod-`9`
grouped scan.

Verdict: the two central branches are certified through scaled `3^8`.
The sampled `3^9` and `3^10` rows are reproduced for both branches, but they
should remain labelled representative until a deeper finite-tail or symbolic
closure argument replaces the representative recurrence.
