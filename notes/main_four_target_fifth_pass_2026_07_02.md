# Four target fifth pass

Date: 2026-07-02.

This pass continues `main_four_target_fourth_pass_2026_07_02.md`.  The user
asked to pursue all next directions with subagents.

## Active workers

```text
Hilbert   Z/35 compressed-state obstruction automaton
Jason     Z/5 x Z/5 b2=0 global/local sieve from smooth F7/F11 charts
Popper    Z/48 bounded cubic-contact/exact-3 production scout
Dalton    A(2,24) alternate/broader scout beyond the cold height-5 split scan
```

## Local Z/35 sanity check

I ran the existing experimental compressed path in

```text
code/agent_Z35_central_branch_deep_lift.m
```

with:

```text
magma -b max_k:=5 direct_depth:=3 use_compressed:=1 \
  max_store:=100000 recon_height:=20 sample_parent_limit:=10 \
  code/agent_Z35_central_branch_deep_lift.m
```

The grouped transition agrees with the exact smoke counts for both central
branches.  For each of `t=1` and `t=2`:

```text
k=4: 6561 lifts; 243 parents have 27 children and 486 have 0.
k=5: 177147 lifts; all 6561 depth-4 states have 27 children.
```

The compressed rows per branch were:

```text
<4, 81, 6561, 721, 729, 243, 486,
    [<"0",486>, <"27",243>],
    [<"[ 0, 0 ]",243>, <"[ 1, 2 ]",486>], 27>

<5, 243, 177147, 27, 6561, 6561, 0,
    [<"27",6561>],
    [<"[ 0, 0 ]",6561>], 0>
```

This validates the compressed machinery through one nontrivial obstructed
transition and one full-lift transition.  Hilbert's task is to turn this into
a true compressed-state automaton beyond the existing `direct_depth >=
max_k-direct_depth` limitation and to separate certified counts from
representative recurrence diagnostics.

## Current expected priority

```text
1. Z/35 compressed-state automaton: highest priority.
2. Z/5 x Z/5 b2=0 global/local sieve: second priority.
3. Z/48 cubic-contact production: bounded background.
4. A(2,24) alternate scout: low-priority sanity check.
```

## Results landed in this pass

### Z/35 compressed automaton

New files:

```text
code/agent_Z35_compressed_automaton.m
notes/agent_Z35_compressed_automaton.md
results/Z35_compressed_smoke_k4.log
results/Z35_compressed_markov_smoke_k4.log
results/Z35_compressed_finite_tail_k6.log
results/Z35_compressed_finite_tail_k7.log
results/Z35_compressed_finite_tail_k8.log
results/Z35_compressed_k10.log
```

The compressed Z/35 lane is now certified through scaled `3^8` for both
central `b=0,r=1` branches.  Per branch:

```text
k  modulus  lifts
1  3        1
2  9        27
3  27       729
4  81       6561
5  243      177147
6  729      1594323
7  2187     43046721
8  6561     129140163
```

Two-branch totals:

```text
k  modulus  total_lifts
1  3        2
2  9        54
3  27       1458
4  81       13122
5  243      354294
6  729      3188646
7  2187     86093442
8  6561     258280326
```

The finite-tail runs show that the visible five-symbol projected state
`(t,c,obs)` is not Markov by itself.  The all-or-none child-coset rule,
however, survived every certified scan through `k=8`: each encountered
liftable parent coset had either `0` or `27` children with zero next
obstruction.

The `k=9` and `k=10` rows in `Z35_compressed_k10.log` are explicitly
representative diagnostics, not certified counts:

```text
per branch k=9:  3486784401
per branch k=10: 31381059609
```

So the next Z/35 move is not more raw depth.  It is to prove the observed
all-or-none child-coset rule symbolically, or to enlarge the finite state
enough to make the transition genuinely Markov.

### Z/5 x Z/5, b2=0 global/local sieve

New files:

```text
code/agent_z5x5_b2zero_global_sieve.m
notes/agent_z5x5_b2zero_global_sieve.md
results/z5x5_b2zero_global_default.log
results/z5x5_b2zero_global_frontier1M.log
```

The selected smooth `F_7/F_11` charts remain locally alive, but the combined
CRT constraints are strong.  The low CRT search over `7^4*11^3` tested
`122544` tuple products across twelve bounded boxes/sign choices; the first
added good-prime sieve killed every tuple, so `exact_hits=0`.

The high-power CRT frontier uses unique lifts modulo

```text
7^7 * 11^6 = 1458956660623.
```

For all four `F_7/F_11` sign combinations, the tuple coordinate height is
forced above `1000000`; in particular `K`, `s`, and `b1` have no rational
representative of naive height at most `1000000`.  Adjacent probes
`F7 h1=2,h2=0` and `F11 h1=1,h2=1` also lift uniquely, so they are smooth
local points rather than cheap local obstructions.

### Z/48 cubic-contact/exact-3 production

New files:

```text
code/agent_Z48_cubic_production_driver.m
notes/agent_Z48_cubic_production.md
results/Z48_cubic_prod_rt5_shell_s12_c0.log
results/Z48_cubic_prod_rt5_shell_s12_c1.log
results/Z48_cubic_prod_rt5_shell_s12_c2.log
results/Z48_cubic_prod_rt5_shell_s12_c3.log
results/Z48_cubic_prod_fixed_survivor_s12.log
```

The RTHeight `5` shell sample with `SearchBound:=12` was cold:

```text
runSlices=40
tested=1339560
rationalRoots=4
nonsingular=1
pointGatePass=0
threePass=0
certified=0
z48Hits=0
```

The known old point-count survivor

```text
r=-1/4, t=-1/4, mu=-1/2, y=-5/8, p=-41/144, N=5/8, z=125/96
```

again passes the point-count gate, but the exact 3-part layer rejects it:

```text
torsion=[16], exponent=16, has3=false
```

So no Z/48 hit was produced, and the false survivor is now cleanly explained
as a pure 2-primary torsion case.

### A(2,24) alternate scout

New files:

```text
code/agent_A2_24_alt_scout.m
notes/agent_A2_24_alt_scout.md
results/A2_24_alt_syntax.log
results/A2_24_alt_h6_shell_full.log
```

The full height-6 shell beyond the closed height-5 box was small enough to
run completely:

```text
active_checked=43136
split_fibers=4
order12_split_fibers=4
translated_order12_rows=32
low_rows_le_4=0
rational_M_root_rows=0
exact_divisible_rows=0
torsion_cert_rows=0
```

The only new split/order-12 fibers were

```text
(-1/5, +/-6/5,  5/4)
( 1/5, +/-6/5, -5/4)
```

All four have residual quartic type `[1,1,2]` and eight rational
2-torsion translations.  The 32 translated rows had saturated affine
degree distribution:

```text
[16]  : 24 rows
[8,8] : 8 rows
```

No branch had degree `<=4`, no rational `M`-root appeared, and no exact
halving/torsion certification triggered.

## Updated priority after this pass

```text
1. Z/35: prove the all-or-none child-coset rule or build the enlarged Markov state.
2. Z/5 x Z/5: convert the selected-chart CRT height obstruction into a broader chart cover/sieve.
3. Z/48: continue only as bounded background; exact 3-part filtering is working.
4. A(2,24): stop naive split height-box expansion unless used as background coverage.
```
