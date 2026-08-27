# agent_m18_416_den7 notes

Goal: complementary `[4,16]` search in the `M_1(8,4)` family for rational
parameters whose `R` or `w` denominator is divisible by `7`.  These are p=7
projective-boundary cases skipped by strict affine CRT residue vectors.

Driver:

```text
agent_m18_416_den7_boundary_search.m
```

Logs:

```text
agent_m18_416_den7_h20.log
agent_m18_416_den7_h30.log
```

The requested `agent_m18_416_search_notes.md` and
`agent_m18_416_search_crt.m` existed as unmaterialized placeholders in this checkout:
`ls -l` showed nonzero document sizes, but `wc -l` and `xxd` saw no materialized
bytes.  I used Filip's readable notes and scripts as the source of formulas.

## Method

The p=7 denominator charts are treated as projective boundary charts:

```text
Rinf_w      R = 1/z, z = 0, w finite mod 7
r_Winf      w = 1/z, z = 0, R finite mod 7
Rinf_Winf   R = 1/zR, w = 1/zW, zR = zW = 0
```

For each chart the script rebuilds the cleared FIRST_COVER and TARGET_416
coefficient equations over `F_7`, evaluates the special fiber, and brute-forces
the auxiliary variables.  This is a closure test, not a Hensel proof, but it is
the same cleared-boundary style as `m18_m14_second_halving_boundary_equations.m`.

The only live p=7 denominator stratum found is:

```text
r_Winf:1    R finite with R == 1 mod 7, w = infinity
```

All `R = infinity` finite-`w` strata and the simultaneous infinity corner are
empty already for FIRST_COVER closure.

Auxiliary good-prime filters were applied at `11,13,17`: if both denominators
are units and the reduction is good-open, require `Tx/2` and `P_R/2` over the
finite Jacobian; boundary/bad reductions at those primes are allowed through.

For exact checking, `T_x` is checked by the exact first-cover equations over
`Q`.  I initially used Magma's `IsDivisibleBy(Tx,2)`, but it hit an internal
`saturation.m` uninitialized-variable error on the first smooth denominator
candidate.  `P_R` is still checked by exact Magma Jacobian divisibility, with a
rational-model fallback; no `P_R` divisibility errors occurred.

## Runs

Height 20 command:

```text
magma -b height:=20 aux_primes:="11,13,17" progress_interval:=100 max_hits:=20 agent_m18_416_den7_boundary_search.m > agent_m18_416_den7_h20.log 2>&1
```

Height 20 summary:

```text
parameters 511
total_pairs 261121
den7_pairs 52272
p7_chart_pass 3510
aux_pass 1066
smooth 1012
exact_tests 1012
Tx_halves 0
PR_halves 0
both_halves 0
hits 0
```

Height 30 command:

```text
magma -b height:=30 aux_primes:="11,13,17" progress_interval:=500 max_hits:=20 agent_m18_416_den7_boundary_search.m > agent_m18_416_den7_h30.log 2>&1
```

Height 30 summary:

```text
parameters 1111
total_pairs 1234321
den7_pairs 287592
p7_chart_pass 18906
aux_pass 6052
smooth 5914
exact_tests 5914
Tx_halves 0
PR_halves 0
both_halves 0
hits 0
```

Height 30 chart counts:

```text
Rinf_w:*      all p7_pass 0
Rinf_Winf     p7_pass 0
r_Winf:1      seen 18906, p7_pass 18906, exact 5914, Tx_half 0, PR_half 0
r_Winf:other  p7_pass 0
```

Conclusion: through height 30, the denominator-at-7 complement has only one
live p=7 projective closure stratum, `R == 1 mod 7, w = infinity`; after
auxiliary filters and exact checks on all 5914 smooth candidates, no `T_x`
halves, no `P_R` halves, and no `[4,16]` hits occurred.
