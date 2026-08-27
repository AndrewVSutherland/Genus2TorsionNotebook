# Seven-torsion hybrid attempts

This records quick necessary-condition tests for combining existing even-torsion
bases with rational `7`-torsion.

For rational `7`-torsion, every good reduction prime `p != 7` must satisfy

```text
7 | #J(F_p).
```

The computations below use this as a filter before any exact torsion test.


## `[4,4] + 7` on the A(2,2,4,4) tuple list

The script is

```text
code/a2244_plus7_tuple_sieve.m
```

It tests the existing tuple file

```text
paper/scripts_and_data/tor2244.txt
```

where a tuple `[a,b,c,d]` represents

```text
y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2).
```

The reduction-only run found 8 survivors among 26653 tuples:

```text
checked 26653
smooth 26653
survivors 8
```

The survivors were then tested exactly:

```text
data/a2244_plus7_tor2244_exact.txt
```

All eight have exact torsion

```text
[2,2,4,4]
```

and none has rational `7`-torsion.  Thus this known A(2,2,4,4) tuple list gives no `[4,4]+7` example.

The first-kill distribution in the reduction sieve was:

```text
13: 2795
17: 4241
19: 3232
23: 4834
29: 5315
31: 3153
37: 1900
41: 825
43: 164
47: 64
53: 94
59: 19
61: 9
```


## `12 + 7` on the M(12) family

The script is

```text
code/m12_plus7_search.m
```

The finite-field scan of the M(12) `(z,r)` chart is saved in

```text
data/m12_plus7_finite.txt
```

It already shows a boundary obstruction at small primes:

```text
p=3: total 2,  good 0, seven_possible 0
p=5: total 12, good 6, seven_possible 0
```

For larger primes there are good affine points with `7 | #J(F_p)`, so this is
not a global finite-field emptiness statement.  But any rational point in this
chart must reduce to bad or boundary behavior at `p=5`.

The bounded rational search

```text
magma -b mode:="search" height:=15 max_hits:=20 max_tests:=100 code/m12_plus7_search.m
```

gave

```text
checked 81510
smooth 80940
split_T 80940
seven_survivors 0
torsion_tests 0
hits 0
```

A height-30 run was interrupted after about two minutes with no output beyond
the header, so it is not recorded as a completed search.  A height-15
prime-by-prime diagnostic was also interrupted because it was slower than the
search itself.

Conclusion: `12+7` looks locally constrained in the current M(12) chart, with
the first visible obstruction at `p=5`.  It is not disproved by this test, but
it would require boundary analysis rather than height extension.


## `[2,2,2,2] + 7` on the full-split genus-2 chart

The broader full rational 2-torsion test uses the normalized full-split model

```text
y^2 = x*(x-1)*(x-a)*(x-b)*(x-c),
```

with distinct `a,b,c` away from `0,1`.  This is the correct test for
`[2,2,2,2]+7`; it is much less restrictive than the previous A(2,2,4,4)
tuple list.

The script is

```text
code/m2222_plus7_search.m
```

The finite good-chart test is saved in

```text
data/m2222_plus7_finite.txt
```

It shows immediate local restrictions:

```text
p=5:  total 1,   good 1,   seven_possible 0
p=13: total 165, good 165, seven_possible 0
```

Thus any rational full-split example with rational `7`-torsion must reduce to
bad/boundary behavior at both `p=5` and `p=13` in this normalized chart.  Other
larger primes have many good residues with `7 | #J(F_p)`, so this is a
boundary-forcing obstruction, not a proof of global emptiness.

The rational height searches were:

```text
magma -b mode:="search" height:=8  code/m2222_plus7_search.m
magma -b mode:="search" height:=10 code/m2222_plus7_search.m
```

with results:

```text
height 8:
  params 85
  checked 98770
  smooth 98770
  reduction_survivors 0
  exact_tests 0
  hits 0

height 10:
  params 125
  checked 317750
  smooth 317750
  reduction_survivors 0
  exact_tests 0
  hits 0
```

The saved output files are

```text
data/m2222_plus7_h8.txt
data/m2222_plus7_h10.txt
```

Conclusion: `[2,2,2,2]+7` was the right broader test, but the normalized
full-split chart is strongly boundary-constrained by `p=5` and `p=13`.  A
height-10 rational search found no candidate even reaching exact torsion
testing.  The next serious version would need to parameterize the simultaneous
bad/boundary reductions at `5` and `13`, rather than simply increasing height.


### Boundary-focused `[2,2,2,2] + 7` search at `p=5,13`

Because the open full-split chart has no good `7`-possible points at `p=5`
and `p=13`, I added the boundary-focused script

```text
code/m2222_plus7_boundary_search.m
```

For a rational parameter `q`, it records whether `q` reduces to `0`, `1`,
`infinity`, or collides with another one of `a,b,c` modulo `p`.  It then keeps
only triples that are boundary at both `5` and `13`, and tests the remaining
primes for the necessary condition `7 | #J(F_p)`.

The height-10 run

```text
magma -b height:=10 code/m2222_plus7_boundary_search.m
```

saved to

```text
data/m2222_plus7_boundary_h10.txt
```

found

```text
checked 317750
smooth 317750
boundary_5 310150
boundary_5_13 116458
reduction_survivors 0
exact_tests 0
hits 0
```

The larger height-12 run

```text
magma -b height:=12 combo_min_count:=1000000 code/m2222_plus7_boundary_search.m
```

saved to

```text
data/m2222_plus7_boundary_h12.txt
```

found

```text
checked 971970
smooth 971970
boundary_5 942274
boundary_5_13 361774
reduction_survivors 0
exact_tests 0
hits 0
```

Among the `361774` height-12 triples on the required `5` and `13` boundaries,
the first remaining good prime killing them was distributed as

```text
11: 70174
17: 188110
19: 69384
23: 24707
29: 7919
31: 1289
37: 148
41: 43
```

So the simultaneous `5`/`13` boundary condition alone is not enough.  Once it
is imposed, the auxiliary good-prime obstruction, especially at `17`, still
kills every tested candidate before exact torsion is reached.


### Refining the `[2,2,2,2] + 7` boundary search at `11` and `17`

The height-12 `5/13` boundary run showed many first kills at `11` and `17`, so
I separated those primes explicitly in

```text
code/m2222_plus7_boundary_refine_11_17.m
```

The run

```text
magma -b height:=12 code/m2222_plus7_boundary_refine_11_17.m
```

saved to

```text
data/m2222_plus7_boundary_refine_11_17_h12.txt
```

found

```text
checked 971970
smooth 971970
boundary_5_13 361774
compatible_11 291600
compatible_11_17 103490
post_survivors 0
exact_tests 0
hits 0
```

Here `compatible_11` means the triple is either bad/boundary at `11` or has
good reduction with `7 | #J(F_11)`, and similarly for `17`.

The status counts after the forced `5/13` boundary were:

```text
11: boundary    269654
11: good_kill7  70174
11: good_pass7  21946

17: boundary     73348
17: good_kill7  188110
17: good_pass7   30142
```

Among the `291600` triples compatible at `11`, the joint `11/17` statuses were:

```text
11 boundary,    17 boundary:    67450
11 boundary,    17 good_kill7: 174410
11 boundary,    17 good_pass7:  27794
11 good_pass7,  17 boundary:     5898
11 good_pass7,  17 good_kill7:  13700
11 good_pass7,  17 good_pass7:   2348
```

Thus `103490` triples survived both `11` and `17`, but all were killed by later
primes:

```text
19: 69384
23: 24707
29: 7919
31: 1289
37: 148
41: 43
```

So the obstruction moves from `5/13` boundary to `11/17` compatibility and then
to `19/23`.  Even after retaining every bad/boundary case at `5,13,11,17`, no
height-12 triple reaches exact torsion testing.


## Practical conclusion

The broader `[2,2,2,2]+7` test is more appropriate than `[4,4]+7`, but it is
also locally constrained: good full-split reductions at `p=5` and `p=13` are
impossible for rational `7`-torsion.  The existing A(2,2,4,4) tuple list gives
no `[4,4]+7` example after exact torsion, and the `12+7` route is also forced
to boundary at small primes.  None of these routes currently looks promising by
plain height extension.
