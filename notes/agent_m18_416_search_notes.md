# FRONT 1B: p=7-boundary CRT search for `[4,16]`

This is a rational-search continuation inside the `M_1(8,4)` family.  The
finite diagnostic says the good affine `[4,16]` locus is empty over `F_7`, so
the search starts from the cleared `p=7` boundary closure rather than from a
blind `(R,w)` height box.

## Script

`agent_m18_416_search_crt.m`

The script:

1. Recomputes the cleared `p=7` boundary closure for `TARGET_416`.
2. Keeps only affine rational parameters whose reductions `(R,w) mod 7` are in
   that closure.
3. Optionally applies larger-prime good-open filters: if a candidate is good
   and nonboundary modulo an auxiliary prime, its residue must be an open
   `[4,16]` residue at that prime; boundary/bad residues are allowed through.
4. Runs exact rational checks only after the residue filters:
   - exact first halving of `T_x=[x,0]`;
   - exact halving of `P_R=(-R,Y_R)`;
   - exact torsion computation only if `P_R` halves.

The current script deliberately does not enumerate parameters with denominator
divisible by `7`; those are affine-infinity/boundary-chart cases for the
separate `p=7` boundary-lift/blowup worker.

## p=7 closure

The script reproduces the known closure counts:

```text
p7_boundary_closure boundary 41 first_closure 24 target416_closure 19 allowed_pairs 19
```

Allowed affine residue pairs are:

```text
<0,0>, <0,2>, <0,3>, <0,4>, <0,5>,
<1,0>, <1,1>, <1,2>, <1,3>, <1,4>, <1,5>, <1,6>,
<3,0>, <3,3>, <3,4>,
<5,2>, <5,5>,
<6,1>, <6,6>
```

The live strata are the `R`, `R-1`, `R-1&w`, `Q&w`, `Q&Quartic&R&R+w&R-w&w`,
the `Lplus/Lminus` intersections with `R` or `R+-w`, and the four
`Q&Quartic` corner intersections at `(R,w)=(1,1),(1,6),(6,1),(6,6)`.

## Logged runs

Commands run:

```text
magma -b height:=20 max_exact_tests:=20 aux_primes:="" progress_interval:=5000 agent_m18_416_search_crt.m > agent_m18_416_search_h20_p7only.log

magma -b height:=20 max_exact_tests:=20 aux_primes:="11,13" progress_interval:=5000 agent_m18_416_search_crt.m > agent_m18_416_search_h20_aux11_13.log
```

Results:

```text
agent_m18_416_search_h20_p7only.log:
  parameter_count 511
  p7_affine_parameters 457
  p7_denominator_bad_parameters 54
  p7_crt_pair_budget 78376
  candidates_after_aux 70671
  family_smooth 68803
  first_possible 73
  tangent_bases 20
  exact_tests 20
  first_verified 20
  pr_halved 0
  hits 0

agent_m18_416_search_h20_aux11_13.log:
  aux_prime_profile 11: good_open 28, target416 4
  aux_prime_profile 13: good_open 56, target416 4
  AUX_KILL_COUNTS: 11 -> 13972, 13 -> 17960
  candidates_after_aux 46444
  family_smooth 44246
  first_possible 40
  tangent_bases 5
  exact_tests 5
  first_verified 5
  pr_halved 0
  hits 0
```

No `PR_HALF` or `HIT_416` lines occurred in either run.  Thus this pass found
filtered rational candidates and verified exact first-halving points, but no
candidate where `P_R` is divisible by `2`, hence no `[4,16]` hit.
