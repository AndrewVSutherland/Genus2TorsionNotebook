# M_1(8,4) from the M_1(8,2^w) odd model

This records the computation for halving the rational 2-torsion class
`W_0 - infinity` in the `M_1(8,2^w)` family.

## Family

Use the odd model from `NotesAndTodo.tex`:

```text
C: y^2 = x*A(x)*B(x)
```

where

```text
A = n^4*x^2
  + (m^3*n + 4*m^2*t + m*n^3 - 8*m*n*t + 4*n^2*t)*x
  + m^4

B = (m*n + 2*n^2 + 4*t)*x^2
  + (m^2 + 4*m*n + n^2 + 8*t)*x
  + (2*m^2 + m*n + 4*t).
```

The class to halve is represented by `[x,0]`.

## Tangent Equation

Write `h=A*B`, with leading coefficient `c4`. A half exists if there are
rational `U,V,M,N` such that

```text
h(x) - x*(M*x+N)^2 = c4*(x^2 + U*x + V)^2.
```

The constant equation is

```text
h(0) = c4*V^2.
```

On the affine chart `n=1`, put `R=m/n`. Introduce `w` by the first square
condition:

```text
(2*R^2 + R + 4*t)/(R + 2 + 4*t) = w^2.
```

Solving for `t` gives

```text
t = (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1)).
```

Then `V = +/- R^2*w`. After eliminating `M,N` by squaring, the remaining
equation factors into two quadratics in `U`. Their discriminants are:

```text
Delta_plus =
  -4*(w-R)*(R-1)^2*(R+1)*(R+w)*(w+1)
    *(R*w - 3*R + 3*w - 1)

Delta_minus =
   4*(w-1)*(R+1)*(R*w + 3*R + 3*w + 1)
    *(R^4 - 2*R^3 + R^2*w^2 - R^2 + 2*R*w^2 - w^2).
```

These discriminant-square conditions are necessary but not sufficient by
themselves: after choosing a rational `U`, one still has to require
`M^2 = c3 - 2*c4*U` and `N^2 = c1 - 2*c4*U*V` to be rational squares.

## Scripts

```text
code/m18_m14_square_condition_test.m
```

shows that the first square condition alone is not sufficient.

```text
code/m18_m14_tangent_search.m
```

solves the full tangent equations and verifies exact divisibility with Magma.

```text
code/m18_m14_cover_search.m
```

is the earlier exploratory discriminant-cover search; it can have false
positives because it does not impose the final `M^2,N^2` square conditions
before verification.

## Search Results

The direct tangent search:

```text
magma -b height:=20 max_hits:=80 code/m18_m14_tangent_search.m
```

stopped at the requested hit cap:

```text
checked 103688
tangent_points 80
verified 80
hits 80
```

Every verified hit had `half_order 4` and torsion subgroup `[4,8]`. Many
also had irreducible Frobenius certificates for simplicity.

One simple example is:

```text
R = -8,  w = 6,  t = 12/5
```

giving

```text
y^2 =
  18/5*x^5 + 24489/25*x^4 + 708048/25*x^3
  + 6179904/25*x^2 + 2654208/5*x.
```

Equivalently, after scaling `y`,

```text
y^2 =
  90*x^5 + 24489*x^4 + 708048*x^3
  + 6179904*x^2 + 13271040*x.
```

Magma computes torsion `[4,8]`. A simplicity certificate is given by the
irreducible Frobenius polynomial at `p=43`:

```text
1849*T^4 + 38*T^2 + 1.
```


## One-split follow-up: `[2,4,8]`

The full-split follow-up below asks both quadratics `A` and `B` in

```text
C: y^2 = x*A(x)*B(x)
```

to split, which targets the stronger `[2,2,4,8]` situation.  A weaker and much
more natural test is to require only one of `A,B` to split.  This should add
one independent rational `2`-torsion class to the infinite `[4,8]` family and
therefore target exact torsion `[2,4,8]`.

I added

```text
code/m18_m14_one_split_search.m
```

The height-20 run

```text
magma -b height:=20 max_hits:=20 code/m18_m14_one_split_search.m \
    > data/m18_m14_one_split_h20.txt
```

found `8` verified hits, all with exact torsion

```text
J(Q)_tors = [2,4,8].
```

The summary was:

```text
checked        259080
one_split        1132
both_split          0
tangent_points      8
verified            8
hits                8
```

A geometrically simple example is obtained from

```text
R = -16/11,  w = 14/11,  t = 42/55.
```

One integral odd model used for the exact torsion computation is

```text
y^2 =
  7061463847622250*x^5
  + 104632219276049025*x^4
  + 135735215960638800*x^3
  + 188573481843278400*x^2
  + 51200550567936000*x.
```

Here `A` splits and `B` does not:

```text
A = x^2 + (8392/605)*x + 65536/14641,
B = (18/5)*x^2 + (2061/605)*x + 3528/605,
disc(A) = 63872064/366025,
disc(B) = -5297643/73205.
```

Magma gives

```text
J(Q)_tors = [2,4,8].
```

It also gives a reduced Weierstrass model for the same curve:

```text
y^2 + (x^2 + x)*y =
  60*x^5 + 3982*x^4 + 1100*x^3 + 130682*x^2 - 239760*x.
```

At `p=47`, the Frobenius characteristic polynomial is

```text
T^4 - 4*T^3 + 30*T^2 - 188*T + 2209.
```

It is irreducible, and its `12`th-power transform is also irreducible:

```text
Y^4 - 1253509316*Y^3 - 58850623665731800826*Y^2
    - 145647106516923685790048393156*Y
    + 13500460747057082764996435506735298654081.
```

Thus this gives a geometrically simple Jacobian with rational torsion
`[2,4,8]`.  This is a useful positive result: the `[2,4,8]` target is easy in
the `[4,8]` family, while the stronger full-split `[2,2,4,8]` target remains
the hard one.

## Full-Split Follow-Up

`code/m18_m14_full_split_search.m` adds the direct full-rational-2 test inside the same tangent locus. It requires both quadratics `A` and `B` in

```text
C: y^2 = x*A(x)*B(x)
```

to split over `Q`, then re-runs the exact tangent verification. The first run found no full-split points:

```text
magma -b height:=20 max_hits:=20 code/m18_m14_full_split_search.m

checked 259080
split_quadratics 0
tangent_points 0
verified 0
hits 0
```

The splitting conditions after the first square parameter `w` are already quite restrictive. On the chart `n=1`, the two quadratic discriminants are:

```text
disc(A) =
  -4*(w-R)*(R-1)^2*(R+w)
    *(R^4 - 2*R^3 + R^2*w^2 - R^2 + 2*R*w^2 - w^2)
    /((w-1)^2*(w+1)^2)

disc(B) =
  (R-1)^2*(R*w - 3*R + 3*w - 1)*(R*w + 3*R + 3*w + 1)
    /((w-1)*(w+1)).
```

A small finite-field check found no nonboundary split points over `F_3` or `F_5`; for larger primes the split locus is nonempty. Thus this chart is not globally locally empty, but any rational full-split point found this way is forced into small-prime boundary behavior.

`code/m18_m2228_curve_extra_halving_search.m` checks the known rational curve on the full-split `M_1(8,2,2,2)` surface for an additional independent halving. The first run

```text
magma -b height:=20 max_hits:=10 code/m18_m2228_curve_extra_halving_search.m

checked 508
smooth 508
qsplit 508
base_torsion 508
hits 0
```

found only the expected torsion `[2,2,2,8]` on that one-parameter slice, not the desired `[2,2,4,8]`.

The broader repo already has a stronger direct search for this same target in `code/m2248_surface_cover_search.m`, `code/m2248_sieve.m`, and `notes/a2244_search_summary.txt`. Those computations show good-reduction obstructions at `11` and `23`, so the useful search is a forced-boundary height enumeration with the signed depth-4 real-compatible filter. Continuing that search at height `20000` gave:

```text
python3 code/run_surface_tuple_chunks.py --B 20000 \
  --output data/surface_tuples_B20000_N_depth4_real_earlypart_11_23_ax1_20.txt \
  --boundary 11:N,23:N --local-depth 4 \
  --chunk-size 10 --start 11 --end 20 --resume

ax=11..20 return=0 elapsed=44.5s rows=0
```

Together with the existing `ax=1..10` pilot, this gives no signed depth-4 real-compatible survivors for `ax <= 20` at height `20000`.


The explicit `B`-split reparameterization is implemented in `code/m18_m14_bsplit_cover_search.m`.  It writes

```text
R = (-3*w^2*z + 2*w*z^2 + 2*w - z)/(w^2*z - 2*w*z^2 - 2*w + 3*z)
```

from the condition that the normalized `B` quadratic has roots with product `w^2`.  The first run again found no full-split points:

```text
magma -b height:=20 max_hits:=20 code/m18_m14_bsplit_cover_search.m

checked 259080
valid_chart 259076
a_split 0
tangent_points 0
verified 0
hits 0
```

So the failure of the direct `(R,w)` scan is not just an artifact of missing small `B`-split parametrizations.


Continuing the chunked K3 boundary search with the common output

```text
data/surface_tuples_B20000_N_depth4_real_earlypart_11_23.txt
```

and common chunk directory

```text
data/surface_tuples_B20000_N_depth4_real_earlypart_11_23.txt.chunks
```

keeps the same settings: height `B=20000`, boundary `11:N,23:N`, signed depth-4 real-compatible filter, and 10-wide `ax` chunks.  The completed chunks through `ax=100` are all empty:

```text
ax=1..10    rows=0
ax=11..20   rows=0
ax=21..30   rows=0
ax=31..40   rows=0
ax=41..50   rows=0
ax=51..60   rows=0
ax=61..70   rows=0
ax=71..80   rows=0
ax=81..90   rows=0
ax=91..100  rows=0
```

The merged output has `0` unique rows, so there were no candidates to send to the exact `m2248` sieve in this range.

## Attempt to add 5-torsion to the `[4,8]` family

The next hybrid attempt was to use the `M_1(8,4)` tangent-cover family, which already gives many simple Jacobians with torsion `[4,8]`, and filter for possible rational `5`-torsion.  The implementation is the `m14_search` mode of

```text
code/elkies5_4_6_hybrid_search.m
```

The necessary finite-field condition is that for every good reduction prime `p`, the Jacobian order over `F_p` must be divisible by `5`.  The script lets bad or boundary reduction pass; it only kills a candidate when it has good reduction and `5` does not divide `#J(F_p)`.

The finite-density check was

```text
magma -b mode:="finite_all" code/elkies5_4_6_hybrid_search.m
```

For the `M_1(8,4)` chart it found good affine residues with `5 | #J(F_p)` at some primes, for example `p=13,23,29,31,37`, but none at `p=7,11,17,19` in the finite good affine chart.  Thus any rational hit in this chart must reduce to bad or boundary behavior at these obstructing primes.

Exact searches:

```text
magma -b mode:="m14_search" height:=18 max_hits:=20 max_tests:=200 code/elkies5_4_6_hybrid_search.m
magma -b mode:="m14_search" height:=30 max_hits:=20 max_tests:=200 code/elkies5_4_6_hybrid_search.m
magma -b mode:="m14_search" height:=50 max_hits:=20 max_tests:=200 code/elkies5_4_6_hybrid_search.m
```

Outputs:

```text
height 18:
  checked 164024
  cover 2016
  smooth 152
  five_survivors 0
  verified 0
  torsion_tests 0
  hits 0

height 30:
  checked 1229880
  cover 5536
  smooth 384
  five_survivors 0
  verified 0
  torsion_tests 0
  hits 0

height 50:
  checked 9566648
  cover 15322
  smooth 922
  five_survivors 0
  verified 0
  torsion_tests 0
  hits 0
```

I added a prime-by-prime diagnostic:

```text
code/elkies5_m14_plus5_prime_diagnostic.m
```

At height `50`, among the `922` smooth tangent-cover candidates, no candidate passed all necessary 5-divisibility reductions.  The strongest killing primes were:

```text
p=7:  good 220, pass5 0, kill5 220
p=11: good 338, pass5 0, kill5 338
p=17: good 482, pass5 10, kill5 472
p=19: good 550, pass5 2, kill5 548
```

The saved data files are:

```text
data/elkies5_4_6_finite_all.txt
data/elkies5_m14_plus5_h18.txt
data/elkies5_m14_plus5_h30.txt
data/elkies5_m14_plus5_h50.txt
data/elkies5_m14_plus5_prime_diag_h30.txt
data/elkies5_m14_plus5_prime_diag_h50.txt
```

Conclusion: this specific `[4,8]` plus `5` chart is not promising by blind height extension.  A rational example would have to lie on simultaneous bad/boundary branches at the small obstructing primes, especially `7`, `11`, and `19`; the height-50 boundary-compatible search still found no candidate even before exact torsion tests.

## Boundary branches for the `[4,8] + 5` attempt

After the exact `[4,8] + 5` search found no `5`-survivors, I analyzed the forced boundary behavior at the small obstructing primes `7`, `11`, and `19`.

The coarse classifier is

```text
code/elkies5_m14_plus5_boundary_branches.m
```

and the factor-level summary is

```text
code/elkies5_m14_plus5_boundary_factor_summary.m
```

For the `n=1` chart with parameters `(R,w)`, after clearing denominators the discriminant factors as

```text
R^8*w^4*(R-1)^14*(R+1)^2*(w-1)^19*(w+1)^19
*(w-R)*(R+w)
*(R*w - 3*R + 3*w - 1)
*(R*w + 3*R + 3*w + 1)
*(R^4 - 2*R^3 + R^2*w^2 - R^2 + 2*R*w^2 - w^2)
*(-2*R^2 + R*w^2 - R + 2*w^2)^8.
```

I label these nontrivial factors as

```text
w-R, R+w, Lplus, Lminus, Qminus, Eresultant.
```

The height-50 factor-level run

```text
magma -b height:=50 code/elkies5_m14_plus5_boundary_factor_summary.m
```

saved to

```text
data/elkies5_m14_plus5_boundary_factor_h50.txt
```

reproduced the same population as the exact search:

```text
checked 9566648
cover 15322
smooth 922
```

Among these `922` smooth rational `[4,8]` candidates, exactly `178` are not killed by the necessary `5 | #J(F_p)` test at all three obstructing primes `7,11,19`; that means they are bad/boundary there, except for two `p=19` boundary-chart points with good reduction and `5 | #J(F_19)`.

However, no candidate survives the full filter through the usual primes:

```text
branch_ok_7_11_19 178
global_ok_all_filter_primes 0
```

The first good prime killing those `178` boundary-compatible candidates is distributed as

```text
13: 52
17: 74
23: 28
29: 12
31: 10
37: 2
```

So the obstruction is not only at `7,11,19`; once those boundary conditions are met, the next primes, especially `13` and `17`, eliminate all candidates in this height-50 sample.

The most common named boundary branches among the `178` candidates are:

```text
p=7:
  R+w & Lplus: 46
  w=inf: 23
  R=-1, w=1, R+w, Lminus, Qminus, Eresultant: 14
  R=inf: 12
  w=0 & Qminus: 11

p=11:
  R=1: 30
  Qminus: 28
  R=inf: 18
  w=inf: 16
  R=0: 13

p=19:
  Eresultant: 44
  R=1: 26
  Lminus: 18
  Qminus: 16
  R=0: 11
  R=inf: 11
```

Conclusion: the boundary analysis confirms that this `[4,8] + 5` route is not failing because of a single missed congruence family.  There are many possible bad/boundary branches at `7,11,19`, but after imposing them the remaining candidates are systematically killed by other good primes.  A blind height extension of the same chart is therefore low value.  If this route is pursued further, the next step should include `13` and `17` in the boundary analysis from the start, or move to a different model for combining even torsion with `5`-torsion.

## Attempt to add 3-torsion to the `[4,8]` family

I also tested the same `M_1(8,4)` tangent-cover family for possible rational `3`-torsion.  If the Jacobian has rational `3`-torsion, then for every good reduction prime `p != 3`, one must have

```text
3 | #J(F_p).
```

The dedicated script is

```text
code/elkies3_m14_plus3_search.m
```

with modes `finite`, `search`, and `prime_diag`.

The finite-density check

```text
magma -b mode:="finite" code/elkies3_m14_plus3_search.m
```

showed immediate small-prime restrictions on the good affine tangent-cover chart:

```text
p=5:  checked 8,   cover 8,   good 0,   three_possible 0
p=7:  checked 24,  cover 22,  good 6,   three_possible 0
p=11: checked 80,  cover 68,  good 16,  three_possible 0
p=13: checked 120, cover 108, good 44,  three_possible 0
p=17: checked 224, cover 196, good 108, three_possible 58
p=19: checked 288, cover 252, good 128, three_possible 64
```

Thus any rational hit must be on bad/boundary behavior at least at `p=5,7,11,13`.

Exact searches were then run with the reduction filter allowing bad reduction to pass:

```text
magma -b mode:="search" height:=30 max_hits:=20 max_tests:=200 code/elkies3_m14_plus3_search.m
magma -b mode:="search" height:=50 max_hits:=20 max_tests:=200 code/elkies3_m14_plus3_search.m
```

Outputs:

```text
height 30:
  checked 1229880
  cover 5536
  smooth 384
  three_survivors 0
  verified 0
  torsion_tests 0
  hits 0

height 50:
  checked 9566648
  cover 15322
  smooth 922
  three_survivors 0
  verified 0
  torsion_tests 0
  hits 0
```

The height-50 prime diagnostic

```text
magma -b mode:="prime_diag" height:=50 code/elkies3_m14_plus3_search.m
```

saved to

```text
data/elkies3_m14_plus3_prime_diag_h50.txt
```

showed the main first killers:

```text
p=5:  good 120, pass3 0,  kill3 120, first_kill 120
p=7:  good 220, pass3 6,  kill3 214, first_kill 194
p=11: good 338, pass3 14, kill3 324, first_kill 206
p=13: good 398, pass3 12, kill3 386, first_kill 174
```

So `p=5` is especially severe: every good-reduction height-50 candidate at `p=5` fails the necessary `3 | #J(F_5)` condition.

I then ran the factor-level boundary summary

```text
code/elkies3_m14_plus3_boundary_factor_summary.m
```

at height `50`, saved to

```text
data/elkies3_m14_plus3_boundary_factor_h50.txt
```

It reproduced the same search population and found

```text
checked 9566648
cover 15322
smooth 922
branch_ok_5_7_11_13 228
global_ok_all_filter_primes 0
```

Thus `228` candidates are not killed at the forced branch primes `5,7,11,13`, but all are killed by later good primes.  The first later killers are

```text
17: 42
19: 68
23: 68
29: 20
31: 16
37: 6
41: 6
43: 2
```

Conclusion: the `[4,8] + 3` version is at least as obstructed as `[4,8] + 5` in this chart.  It has a very strong good-reduction obstruction at `p=5,7,11,13`; after forcing those boundary behaviors, every height-50 candidate is still killed by another good prime, mostly `19` and `23`.  This route does not look promising by height extension.

## Attempt to add 7-torsion to the `[4,8]` family

Since the `M_1(8,4)` tangent-cover family is an infinite source of simple
Jacobians with torsion `[4,8]`, I also tested the fresh target `[4,8] + 7`,
which would give torsion order divisible by `224`.

The implementation is

```text
code/elkies7_m14_plus7_search.m
```

with modes `finite`, `search`, and `prime_diag`.  It uses the necessary
good-reduction condition

```text
7 | #J(F_p)
```

for every good prime `p != 7`.

The finite-density check

```text
magma -b mode:="finite" code/elkies7_m14_plus7_search.m \
    > data/elkies7_m14_plus7_finite.txt
```

showed that the good tangent-cover chart has no `+7` residues at `p=11` and
`p=13`:

```text
p=11: checked 80,  cover 68,  good 16,  seven_possible 0
p=13: checked 120, cover 108, good 44,  seven_possible 0
```

The small primes `p=3` and `p=5` are also degenerate for this chart: at `p=3`
there are no affine parameter pairs after excluding `w=0,+-1`, and at `p=5`
there are no good smooth tangent-cover points.  For larger primes, `+7`
residues do occur, for example at `p=17,19,23,37,43`.

The height-50 search

```text
magma -b mode:="search" height:=50 max_hits:=20 max_tests:=200 \
    code/elkies7_m14_plus7_search.m > data/elkies7_m14_plus7_h50.txt
```

reused the same rational tangent-cover population as the `[4,8]+3` and
`[4,8]+5` searches:

```text
checked         9566648
cover           15322
smooth            922
seven_survivors     0
verified            0
torsion_tests       0
hits                0
```

The prime diagnostic

```text
magma -b mode:="prime_diag" height:=50 \
    code/elkies7_m14_plus7_search.m \
    > data/elkies7_m14_plus7_prime_diag_h50.txt
```

showed the main first killers:

```text
p=5:  good 120, pass7 0,   kill7 120, first_kill 120
p=11: good 338, pass7 0,   kill7 338, first_kill 296
p=13: good 398, pass7 0,   kill7 398, first_kill 232
p=17: good 482, pass7 46,  kill7 436, first_kill 140
p=19: good 550, pass7 100, kill7 450, first_kill 80
```

Thus `[4,8]+7` behaves like the `[4,8]+3` and `[4,8]+5` attempts: the
infinite `[4,8]` family is real and often geometrically simple, but adding odd
torsion in this tangent-cover chart is forced onto simultaneous small-prime
boundary behavior.  The height-50 boundary-compatible sample still has no
global reduction survivor, so a blind height extension is not promising.

## Search for `[8,8]` or `[4,16]` in the `[4,8]` family

The natural next 2-primary question is whether the same `M_1(8,4)`
tangent-cover family contains examples with exact torsion `[8,8]` or
`[4,16]`.  I added

```text
code/m18_m14_torsion_structure_search.m
```

This scans rational `R,w` of bounded height, applies the two discriminant-cover
conditions and the final tangent square conditions, verifies the halving of
`[x,0]` with Magma exact `IsDivisibleBy`, and then computes the full exact
group `J(Q)_tors`.  The current target test is containment-based: it flags any
2-primary structure containing `[8,8]` or `[4,16]`, and also any torsion subgroup
of exponent at least `16`.  The height-50 and height-80 runs below were made
before that cleanup with the exact `[8,8]`/`[4,16]` plus exponent-`16` predicate,
but their `max_order = 64` summaries show that no larger exponent-`8` group
containing `[8,8]` was missed.

The height-50 run

```text
magma -b height:=50 max_tests:=3000 progress_interval:=1000000 \
    code/m18_m14_torsion_structure_search.m \
    > data/m18_m14_torsion_structure_h50.txt
```

gave

```text
checked        9566648
cover            15322
smooth             922
tangent_points     178
verified           178
exact_tests        178
targets              0
max_order           64
max_exponent         8

TORSION_COUNTS
  [2,4,8] 8
  [4,8]   170
```

The deeper height-80 run

```text
magma -b height:=80 max_tests:=10000 progress_interval:=5000000 \
    code/m18_m14_torsion_structure_search.m \
    > data/m18_m14_torsion_structure_h80.txt
```

finished with

```text
checked        61795320
cover             38424
smooth             1784
tangent_points      300
verified            300
exact_tests         300
targets               0
max_order            64
max_exponent          8

TORSION_COUNTS
  [2,4,8] 8
  [4,8]   292
```

Conclusion: in the low-height part of this infinite `[4,8]` family, the only
observed 2-primary enlargement is the already-understood one-split
`[2,4,8]` phenomenon.  No `[8,8]`, no `[4,16]`, and no exponent-`16` torsion
appeared through height `80`.  Since exact torsion tests are sparse and cheap
after the tangent-cover filter, the next improvement should be a
target-specific finite-prime 2-primary compatibility sieve, not simply a blind
height increase.

## Algebraic conditions for `[8,8]` and `[4,16]`

After the height-80 exact torsion scan found no `[8,8]` or `[4,16]`, I switched
from height search to the actual second-halving equations.

For a specialization

```text
C: y^2 = f(x) = x*A(x)*B(x)
```

there are three rational `2`-torsion classes

```text
T_x = [x,0],      T_A = [A,0],      T_B = [B,0],
T_x + T_A + T_B = 0.
```

The diagnostic script

```text
code/m18_m14_twopower_diagnostic.m
```

shows the generic pattern in the `[4,8]` family: all three classes are
`2`-divisible, but only the `T_A` chain is already `4`-divisible.  Thus:

- `[4,16]` means pushing the `T_A` chain one more step.
- `[8,8]` means making the independent half of `T_x` or equivalently `T_B`
  divisible one more time.

The `T_A` chain has a clean rational point.  On the `R,w` chart,

```text
t = (2R^2 + (1-w^2)R - 2w^2)/(4(w^2-1)),
Q = R^2 - (1/2)R w^2 + (1/2)R - w^2,
Y_R = -2R(R-1)^2 Q/(w^2-1).
```

Then

```text
P_R = (-R, Y_R) in C(Q),      f(-R) = Y_R^2,
4(P_R - infinity) = T_A.
```

For example, at `R=-8,w=6`, the diagnostic gives

```text
P_R = (8, 31104/5),
P_R has order 8,
P_R is not divisible by 2,
J(Q)_tors = [4,8].
```

### `[4,16]` condition

Let

```text
q(x) = x^2 + a*x + b,
ell(x) = c*x^2 + d*x + e,
c4 = lc(A*B) = lc(f).
```

Then `P_R - infinity` is divisible by `2` if and only if there are rational
`a,b,c,d,e` such that

```text
f(x) - ell(x)^2 = c4*(x+R)*q(x)^2.          (416)
```

The sign of `ell(-R)=+-Y_R` is irrelevant, since `P_R` is divisible by `2` if
and only if `-P_R` is divisible by `2`.  Therefore a `[4,16]` point in the
`[4,8]` family is obtained by imposing the existing first-cover equations for
halving `T_x`, together with `(416)`.

There is a smaller open-chart form of `(416)`.  Write

```text
h = A*B = c4*x^4 + h3*x^3 + h2*x^2 + h1*x + h0.
```

The constant coefficient of `(416)` is

```text
e^2 + c4*R*b^2 = 0.
```

On the good open chart `b*c4*R != 0`, introduce

```text
ss^2 = -c4*R,       e = eta*ss*b,       eta = +-1.
```

The `x^4` and `x` coefficients then solve for `a` and `d`:

```text
a = (h3 - c^2 - c4*R)/(2*c4),
d = (h0 - c4*b^2 - 2*c4*R*a*b)/(2*eta*ss*b).
```

Thus the reduced `[4,16]` cover over the `[4,8]` tangent cover is cut out by:

```text
FIRST_COVER,
ss^2 + c4*R = 0,
Coeff_x2(416) = 0,
Coeff_x3(416) = 0,
```

with the displayed substitutions for `a,d,e`.  This reduction is generated by

```text
code/m18_m14_416_reduced_conditions.m
```

and the outputs are

```text
data/m18_m14_416_reduced_conditions_summary.txt
data/m18_m14_416_reduced_conditions_full.txt
```

For each sign `eta=+-1`, the generated summary is:

```text
FIRST_COVER:          4 equations, term counts 21, 24, 21, 0
SECOND_REDUCED_416:  3 equations, term counts 4, 375, 135
```

Here the first reduced equation is the square relation
`ss^2 + c4*R = 0`; the other two are the genuine remaining halving equations.

I then ran a targeted search using this reduced necessary condition before exact
Jacobian arithmetic:

```text
code/m18_m14_416_reduced_search.m
data/m18_m14_416_reduced_search_h50.txt
```

The height-50 run checked the same `R,w` box used in the older torsion scans,
but first imposed `ss^2=-c4*R`, then the `[4,8]` tangent-cover test, and only
then exact-tested whether `P_R` is divisible by `2`:

```text
checked        9566648
square416        10984
cover             6296
smooth              12
tangent_points       2
tx_halves            2
pr_tests             2
pr_halves            0
torsion_tests        0
simple_hits          0
```

So the reduced square relation is very restrictive in low height; among the
few smooth `[4,8]` tangent points that satisfy it, none actually halves `P_R`.
This is consistent with the earlier full height-80 torsion scan, which found
`300` exact `[4,8]` tangent points and no target with exponent `16`.

### `[8,8]` condition

The first `[4,8]` cover is represented by variables `U,V,M,N` satisfying

```text
h(x) - x*(M*x+N)^2 = c4*(x^2 + U*x + V)^2,      h=A*B.
```

Put

```text
a0(x) = x^2 + U*x + V,
ell0(x) = x*(M*x+N),
v0(x) = -ell0(x) mod a0(x) = (M*U-N)*x + M*V.
```

Then the half of `T_x` is the Mumford class

```text
H_x = [a0, v0].
```

This class is divisible by `2` if and only if there are rational
`a,b,rho,sigma` such that, with

```text
q(x) = x^2 + a*x + b,
ell1(x) = -v0(x) + a0(x)*(rho*x + sigma),
```

one has

```text
f(x) - ell1(x)^2 = -rho^2*a0(x)*q(x)^2.       (88)
```

Thus a `[8,8]` point is obtained by imposing the first-cover equations and
`(88)`.  This is the right algebraic replacement for blind height search: it is
a finite second-halving cover over the already-derived `[4,8]` cover.

The equation-generation script is

```text
code/m18_m14_second_halving_equations.m
```

and the generated outputs are

```text
data/m18_m14_second_halving_equations_summary.txt
data/m18_m14_second_halving_equations_full.txt
```

The cleared coefficient systems are moderate in size:

```text
FIRST_COVER: 5 equations, one identically zero after clearing.
TARGET_416: 5 equations, term counts 4, 8, 32, 28, 20.
TARGET_88:  6 equations, term counts 4, 20, 62, 52, 41, 8.
```

So the next algebraic move should be to study rational points on these two
second-halving covers, starting with finite-prime and boundary analysis of the
systems `(416)` and `(88)`, rather than increasing the naive `R,w` height box.

## Finite-prime and boundary diagnostics for the second-halving covers

I implemented the finite-prime and boundary analysis requested after deriving
the algebraic second-halving equations.  The scripts are

```text
code/m18_m14_boundary_factors.m
code/m18_m14_second_halving_finite_boundary.m
code/m18_m14_second_halving_boundary_equations.m
```

and the output files are

```text
data/m18_m14_boundary_factors.txt
data/m18_m14_second_halving_finite.txt
data/m18_m14_second_halving_boundary_3_5_7_11_13.txt
data/m18_m14_second_halving_boundary_equations_p3.txt
data/m18_m14_second_halving_boundary_equations_p5_p7.txt
```

### Base boundary factors

For the `R,w` chart, the discriminant numerator factors as

```text
w^4 * R^8 * (R-1)^14 * (R+1)^2 * (R-w) * (R+w)
* (R*w - 3R + 3w - 1)
* (R*w + 3R + 3w + 1)
* (R^2 - (1/2)R*w^2 + (1/2)R - w^2)^8
* (R^4 - 2R^3 + R^2*w^2 - R^2 + 2R*w^2 - w^2).
```

The denominator is `(w-1)^13(w+1)^13`, and the leading coefficient of `A*B`
has numerator `(R-1)(R+1)`.  These are the boundary components used in the
finite boundary classifier.

### Good-open finite-field test

On the good affine chart, I used the finite Jacobian rather than the expanded
coefficient equations:

- `[4,16]` is tested by `T_x` divisible by `2` and `P_R` divisible by `2` in
  `J(F_p)`.
- `[8,8]` is tested by `T_x` divisible by `4` in `J(F_p)`.

The key counts are:

```text
p   good_open  first_Tx_half  target416  target88
3       0             0             0         0
5       0             0             0         0
7       8             8             0         4
11     28            20             4        12
13     56            28             4         8
17    136            64            12        28
19    164           104            16        28
23    288           184            32        44
29    528           256            32       116
31    640           376            56       132
37    968           516            64       228
41   1240           644           112       260
43   1372           792           128       244
```

Conclusions from the good-open test:

- At `p=3` and `p=5`, the entire affine chart is boundary/bad.  Any rational
  `[4,8]`, `[4,16]`, or `[8,8]` example in this chart necessarily has bad or
  boundary reduction at both primes.
- At `p=7`, `[8,8]` has good-open residues, but `[4,16]` has none.  Therefore
  a rational `[4,16]` point in this chart must reduce to the boundary modulo
  `7`.
- From `p=11` onward, both `[4,16]` and `[8,8]` have good-open residues.  There
  is no immediate finite-field obstruction away from the small-prime boundary.

### Boundary classification

The base boundary classifier gives, for the small primes:

```text
p=3:  open_good 0,  boundary_total 9
p=5:  open_good 0,  boundary_total 25
p=7:  open_good 8,  boundary_total 41, open_target416 0, open_target88 4
p=11: open_good 28, boundary_total 93, open_target416 4, open_target88 12
p=13: open_good 56, boundary_total 113, open_target416 4, open_target88 8
```

For `p=3`, every residue is on a multiple boundary intersection.  For `p=5`,
there are also single-factor boundary components such as `R=+-1`, `R=+-w`, and
`w=+-1`, but still no good open chart.  For `p=7`, the `[4,16]` problem is the
interesting one: it has no good-open residue and is forced onto the 41 boundary
base residues.

### Cleared boundary-equation closure

I also evaluated the cleared polynomial systems on boundary base residues for
`p=3,5,7`.  This is not a proof of liftability, since these are singular closed
strata, but it distinguishes a genuine empty boundary from a boundary that
requires blowups/local lifting.

The summaries are:

```text
p=3: boundary 9,  first_closure 6,  target416_closure 6,  target88_closure 6
p=5: boundary 25, first_closure 15, target416_closure 14, target88_closure 9
p=7: boundary 41, first_closure 24, target416_closure 19, target88_closure 18
```

Thus the small-prime obstruction is not that the closed boundary systems are
empty.  Instead:

- `[8,8]` is locally plausible on the good open chart already at `p=7`.
- `[4,16]` is forced to the `p=7` boundary, but the boundary closure still has
  many residues.  A real obstruction, if present, would have to appear after
  blowup/Hensel-depth analysis of those boundary strata.

Practical conclusion: algebraically, `[8,8]` is the more promising of the two
inside this family.  `[4,16]` is more constrained because it has no good-open
`F_7` point; continuing it sensibly means analyzing the `p=7` boundary closure,
not extending rational height search.

## Reduced `[8,8]` halving conditions

I then focused on `[8,8]` alone and rewrote the second-halving equations on the
open first-cover chart.  The starting point is

```text
h - x*(M*x+N)^2 = c4*(x^2+U*x+V)^2,       h=A*B,
```

with

```text
a0 = x^2 + U*x + V,
ell0 = x*(M*x+N).
```

The first half of `T_x=[x,0]` is the Mumford class `[a0,v0]`, where
`v0=-ell0 mod a0`.  A second half exists if

```text
ell1 = ell0 + a0*(rho*x + lambda)
```

and

```text
f - ell1^2 = -rho^2*a0*q^2,       q=x^2+a*x+b.
```

Using `f = ell0^2 + x*c4*a0^2`, this is equivalent to

```text
x*c4*a0 - 2*ell0*(rho*x+lambda)
  - a0*(rho*x+lambda)^2 + rho^2*q^2 = 0.
```

The constant coefficient is

```text
rho^2*b^2 = V*lambda^2.
```

On the good open chart this forces `V` to be a square.  Since the first cover
has

```text
V = eps*R^2*w,       eps = +-1,
```

the rational `[8,8]` locus is contained in the two square subcovers

```text
w = s^2       with eps=+1,
w = -s^2      with eps=-1.
```

After setting `z=1/rho`, `lambda=tau/z`, and `b=eta*R*s*tau` with
`eta=+-1`, the `x^3` coefficient gives

```text
a = U/2 + M*z + tau - (c4/2)*z^2.
```

The remaining second-halving conditions are just the `x` and `x^2`
coefficients of

```text
z^2*x*c4*a0 - 2*z*ell0*(x+tau)
  - a0*(x+tau)^2 + (x^2+a*x+b)^2 = 0.
```

This gives a compact reduced system: the first-cover equations in
`R,s,U,M,N`, plus two second-halving equations in `z,tau`, and the sign
`eta=+-1`.

The implementation is

```text
code/m18_m14_88_reduced_conditions.m
```

with generated outputs

```text
data/m18_m14_88_reduced_conditions_summary.txt
data/m18_m14_88_reduced_conditions_full.txt
data/m18_m14_88_reduced_conditions_factors.txt
```

The two square branches `w=s^2` and `w=-s^2` collapse to the same reduced
equations after pairing with the correct first-cover sign: the curve depends on
`w^2`, and in both cases `V=R^2*s^2`.

The reduced equations have size

```text
FIRST_REDUCED:      3 nonzero equations, term counts 21, 24, 21
SECOND_REDUCED_88: 2 equations, term counts 16, 36
```

No individual equation factors nontrivially in these coordinates.

I also added a finite-field verifier

```text
code/m18_m14_88_reduced_finite.m
```

which solves the reduced equations over finite fields and compares them with
the intrinsic test `T_x in 4J(F_p)`.  For `p=7,11`, the reduced equations match
the intrinsic condition exactly on the square-cover points:

```text
p=7:  good_square_cover_bases 16, intrinsic_Tx_fourdiv 8,  reduced_88 8
p=11: good_square_cover_bases 56, intrinsic_Tx_fourdiv 24, reduced_88 24
```

At `p=13`, the intrinsic finite test has extra points because `-1` is a square
modulo `13`; these correspond to the sign branch `V=-R^2*s^2`, which is not a
rational square branch over `Q`.  The rationally relevant reduced system still
finds `16` square-cover points modulo `13`.

Thus the useful algebraic target for `[8,8]` is no longer the raw
`R,w,U,V,M,N,a,b,rho,sigma` system.  It is the square-subcover system in
`R,s,U,M,N,z,tau` with only two genuine second-halving equations.



## Reduced `[8,8]` exact search

I implemented a target-specific rational search for the reduced square-subcover
system:

```text
code/m18_m14_88_reduced_search.m
```

It enumerates rational `R,s`, puts `w=s^2`, applies the small-prime residue
filter at `p=7,11,13`, solves the first-cover equations exactly, then solves the
two reduced second-halving equations for `z,tau`.  Any solution is verified
again by exact Magma Jacobian arithmetic before computing `J(Q)_tors`.

The unfiltered height-20 run

```text
magma -b height:=20 max_hits:=10 progress_interval:=500000 \
    use_filter:=false code/m18_m14_88_reduced_search.m \
    > data/m18_m14_88_reduced_search_h20.txt
```

gave

```text
checked          261121
open_good        258004
first_bases        1020
first_solutions     4080
second_solutions       0
exact_verified         0
torsion_tests          0
hits                   0
```

The filtered height-40 run

```text
magma -b height:=40 max_hits:=10 progress_interval:=1000000 \
    code/m18_m14_88_reduced_search.m \
    > data/m18_m14_88_reduced_search_h40_filtered.txt
```

gave

```text
checked             3837681
residue_filtered    3813949
open_good             23732
first_bases             150
first_solutions         600
second_solutions          0
exact_verified            0
torsion_tests             0
hits                      0
```

The filtered height-80 run

```text
magma -b height:=80 max_hits:=10 progress_interval:=5000000 \
    code/m18_m14_88_reduced_search.m \
    > data/m18_m14_88_reduced_search_h80_filtered.txt
```

finished with

```text
checked            61826769
residue_filtered   61270575
open_good            556194
first_bases             664
first_solutions        2656
second_solutions          0
exact_verified            0
torsion_tests             0
hits                      0
```

Conclusion: this was a genuine `[8,8]` search, not a blind torsion scan.  The
finite-prime square-cover tests show local points, and the rational search does
find many first-cover points, but none satisfy the reduced second-halving
equations through height `80`.  The next logical move is therefore not another
height extension.  It is to eliminate `M,N,z,tau` from the reduced system, or to
perform a deeper local/boundary analysis on the reduced square subcover, looking
for a global obstruction to rational points.


## Eliminating the reduced `[8,8]` system

I next eliminated the auxiliary variables in the reduced square-subcover system
as far as seemed computationally useful.  The reproducible scripts are

```text
code/m18_m14_88_elimination_diagnostic.py
code/m18_m14_88_residue_profile.m
```

with outputs

```text
data/m18_m14_88_elimination_diagnostic.txt
data/m18_m14_88_residue_profile_7_11_13.txt
data/m18_m14_88_residue_profile_17_43.txt
```

After solving the `x^3` coefficient, the two reduced second-halving equations
can be written in the compact form

```text
E2 = (P(z) + 4*D(z)*tau)/4,
E1 = (2*mu-U)*tau^2 + B(z)*tau + mu^2*c*z^2,
```

where `mu=eta*R*s`, `c=c4`, and

```text
D(z) = 2*mu - U - c*z^2,
P(z) = c^2*z^4 - 4*M*c*z^3
       + (4*M^2 + 2*U*c)*z^2
       + (4*M*U - 8*N)*z + U^2 - 4*mu^2,
B(z) = 2*z*(M*mu-N) - mu*(2*mu-U) - c*mu*z^2.
```

Eliminating `tau` without dividing by `D` gives

```text
G = (2*mu-U)*P^2 - 4*B*P*D + 16*mu^2*c*z^2*D^2.
```

Modulo the first-cover relations

```text
M^2 = c3 - 2*c*U,
N^2 = c1 - 2*c*U*mu^2,
2*M*N = c2 - c*(U^2 + 2*mu^2),
```

the polynomial `G` becomes linear in `M,N`.  Eliminating the tangent signs gives
a single sign-eliminated condition `H_eta(R,s,U,z)=0`.

The first-cover equation on `w=s^2` factors into boundary factors and two
nonboundary branches:

```text
boundary factors: R-1, R+1, s-1, s+1, s^2+1
two U-branches: total degree 9, degree_U 2, 24 terms each
```

The substituted sign-eliminated second-cover condition is much larger:

```text
H_eta_+ : total degree 52, degree_U 10, degree_z 16, 10634 terms
H_eta_- : total degree 52, degree_U 10, degree_z 16, 10634 terms
```

Both `H_eta` polynomials are irreducible over `Q` in these coordinates.  A
direct resultant `Res_U(branch,H_eta)` was tested, but it is too large for
routine interactive use.  Thus this elimination does not reveal a hidden small
factor, boundary component, or low-degree obstruction.

I also profiled the reduced equations over finite fields using the linear
`tau` formula.  The counts below are on the good open rational square branch
`w=s^2`:

```text
p   open_good  first_pairs  reduced_pairs
7       8            8             4
11     28           20            12
13     48           24             8
17    120           56            28
19    164           72            28
23    288          144            44
29    512          208           144
31    640          280           132
37    944          376           260
41   1168          464           276
43   1372          616           244
```

This confirms the earlier `p=7,11,13` filters exactly, but it also shows that
there is no simple finite-prime obstruction through `p=43`: every tested prime
has good-open reduced `[8,8]` residues.

Conclusion: the failure through rational height `80` is not explained by a
visible factor or by local emptiness at small primes.  The reduced `[8,8]`
locus looks like a genuinely high-degree cover of the `[4,8]` family.  The next
useful step, if continuing this exact route, is either a much stronger
multi-prime CRT search using the new residue profiles, or a more geometric
study of the two degree-2 first-cover branches intersected with the irreducible
`H_eta` cover.


## CRT-guided reduced `[8,8]` search

I implemented the multi-prime CRT search suggested by the residue profiles:

```text
code/m18_m14_88_crt_search.m
```

For each chosen prime it recomputes the good-open reduced `[8,8]` residue pairs
`(R,s) mod p`.  It then buckets rational parameters by their full residue vector
and only runs exact arithmetic on rational pairs satisfying all chosen
congruence conditions.  This is a strict good-open search: rational parameters
whose denominators are divisible by one of the CRT primes are excluded, since
those are boundary cases at that prime.

The smoke test with the original three primes was

```text
magma -b height:=40 primes:=7,11,13 max_hits:=5 progress_interval:=100 \
    code/m18_m14_88_crt_search.m \
    > data/m18_m14_88_crt_search_h40_p7_11_13.txt
```

and gave

```text
pair_candidates 528
open_good       528
first_bases      52
first_solutions 208
second_solutions  0
hits              0
```

The six-prime height-500 run

```text
magma -b height:=500 primes:=7,11,13,17,19,23 \
    max_hits:=10 progress_interval:=250 \
    code/m18_m14_88_crt_search.m \
    > data/m18_m14_88_crt_search_h500_p7_23.txt
```

gave

```text
parameter_count        304463
strict_good_parameters 197219
pair_candidates          9336
first_bases               264
first_solutions          1056
second_solutions            0
hits                        0
```

The eight-prime height-1000 run

```text
magma -b height:=1000 primes:=7,11,13,17,19,23,29,31 \
    max_hits:=10 progress_interval:=250 \
    code/m18_m14_88_crt_search.m \
    > data/m18_m14_88_crt_search_h1000_p7_31.txt
```

gave

```text
parameter_count        1216767
strict_good_parameters  743091
pair_candidates           3604
first_bases                248
first_solutions            992
second_solutions             0
hits                         0
```

Finally, the eight-prime height-2000 run

```text
magma -b height:=2000 primes:=7,11,13,17,19,23,29,31 \
    max_hits:=10 progress_interval:=1000 \
    code/m18_m14_88_crt_search.m \
    > data/m18_m14_88_crt_search_h2000_p7_31.txt
```

finished with

```text
parameter_count        4866351
strict_good_parameters 2952341
pair_candidates          51896
open_good                51896
first_bases               1028
first_solutions           4112
second_solutions             0
exact_verified              0
torsion_tests               0
hits                        0
```

Conclusion: the CRT method works computationally and gives a much stronger
negative search than the original height-80 scan, but it still finds no
rational second-halving point.  In particular, among more than fifty thousand
eight-prime-compatible rational pairs of height at most `2000`, over one
thousand exact first-cover bases were reached, yet none satisfied the reduced
second-halving equations.  Further height extension in the same good-open
square branch now looks low value unless paired with a new geometric idea or a
separate boundary-prime search.
