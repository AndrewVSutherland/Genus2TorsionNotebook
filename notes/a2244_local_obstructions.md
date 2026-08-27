# Local obstructions to A(2,2,4,4) on S = M(2,2,2,8)

This note summarizes the local obstruction analysis for the direct K3-surface search.
The input data are the 105 primitive curve tuples in `data/surface_tuples_B2000.txt`.

## Good-reduction criterion

Let the squared branch residues modulo an odd prime be

```text
A = a^2, B = b^2, C = c^2, D = d^2.
```

For an ordering `(A,B,C,D)`, the A(2,2,4,4) square conditions are

```text
(A-C)(A-D), (B-C)(B-D), (A-C)(B-C), (A-D)(B-D) in k^{*2}.
```

At good reduction, write `chi` for the quadratic character. If

```text
x = chi(A-C), y = chi(A-D), z = chi(B-C), w = chi(B-D),
```

then the four conditions are

```text
xy = zw = xz = yw = 1.
```

Thus `x = y = z = w`. Equivalently, after choosing a 2+2 partition of the four branch residues, all four cross-differences must have the same quadratic character. There are only three partitions to test:

```text
12|34, 13|24, 14|23.
```

This is the diagnostic implemented in `code/analyze_a2244_local.py`.

## Height <= 2000 tuple search

For the 105 direct K3 tuples of height <= 2000, all fail before the full M(2,2,4,8) cover: none even pass the local A(2,2,4,4) sieve for the tested primes.

First obstruction distribution, with primes tested in increasing order:

```text
5   23
7   27
11  38
13  3
17  7
19  2
23  4
29  1
```

A small set-cover check shows that, on this finite sample, the three primes

```text
11, 29, 31
```

already obstruct every tuple.

## Structural good-reduction obstructions

The finite-field enumeration is more informative than the first-obstruction list. For nondegenerate points of the K3 surface over finite fields, the counts are:

```text
p   signed_good  signed_good_A2244  square_classes_good  square_classes_A2244
11  240          0                  5                    0
13  288          288                3                    3
17  1920         768                24                   8
19  1728         432                36                   9
23  5280         0                  110                  0
29  12768        4032               196                  42
31  16560        5040               345                  105
```

So `p=11` and `p=23` are genuine good-reduction obstructions: no nondegenerate point of `S(F_p)` lifts locally to A(2,2,4,4) at either prime. By contrast, `p=13` has no good-reduction obstruction in this sense; its obstructions in the height search occur at degenerate reductions.

The small primes `5` and `7` are useful as filters, but they cannot support four distinct nonzero squared branch residues, so their obstruction behavior is boundary/degeneracy behavior rather than a good-reduction obstruction.

## Consequence for future searches

Any rational point on the desired cover must have bad reduction at both `11` and `23`. For an integral representative `[a,b,c,d]`, this means it must satisfy

```text
11 | a b c d * product_{i<j}(a_i^2 - a_j^2)
23 | a b c d * product_{i<j}(a_i^2 - a_j^2).
```

Equivalently, modulo each of `11` and `23`, at least one branch point is zero or two squared branch points collide.

The next search should not be a larger blind height search. It should enumerate points on the K3 surface subject to these forced boundary congruences at `11` and `23`, then apply the exact p-adic/rational square tests. This should cut the search space and target the only residue classes that can survive the structural local obstructions.


## Boundary-stratum search results

I implemented the forced-boundary search in `code/enumerate_surface_tuples.cpp`. The enumerator now accepts required bad primes and, optionally, raw K3 boundary components such as

```text
11:Ead,23:Ecd
```

where the raw K3 coordinates are `[a,b,c,d]` and the components are

```text
Za, Zb, Zc, Zd, Eab, Eac, Ead, Ebc, Ebd, Ecd.
```

The stratum prefilter uses the K3 quadratic in the solved coordinate:

```text
v^2 d^2 + (2uv - 4abc)d + u^2 = 0,
where u = ab + ac + bc and v = a + b + c.
```

Degenerate affine cases are treated conservatively, so the prefilter can admit false positives but should not reject a rational point merely because the chosen affine chart is bad modulo the prime.

Validation and search results:

```text
B=2000, 11:any,23:any:
  38 tuples, exactly matching the previous congruence-filtered file.

B=5000, 11:any,23:any:
  72 tuples.
  Exact all-permutation M(2,2,4,8) sieve: 0 intermediate, 0 full.

B=10000, d-collision family
  11 in {Ead,Ebd,Ecd}, 23 in {Ead,Ebd,Ecd}:
  114 tuples.
  Contains all 72 B=5000 forced-bad tuples and adds 42 new tuples.
  Exact all-permutation M(2,2,4,8) sieve: 0 intermediate, 0 full.
```

This makes the current best interpretation sharper: the forced bad-reduction condition is real, and the `d`-collision boundary family is high-yield, but even this targeted family gives no A(2,2,4,4) source tuple up to height `10000`. The next natural split is to run the remaining boundary families separately, especially zero components and collisions not involving `d`, rather than increasing height uniformly.


## Non-d/zero boundary search results

I next ran the complementary high-priority boundary block at height `10000`:

```text
11 in {Za,Zb,Zc,Zd,Eab,Eac,Ebc},
23 in {Za,Zb,Zc,Zd,Eab,Eac,Ebc}.
```

This covers zero components and collisions not involving the solved coordinate `d`, at both obstruction primes. The merged result was:

```text
B=10000, non-d/zero family:
  118 tuples.
  Intersection with B=10000 d-collision family: 114 tuples.
  New relative to d-collision family: 4 tuples.
  Contains all 72 B=5000 forced-bad tuples.
  Exact all-permutation M(2,2,4,8) sieve: 0 intermediate, 0 full.
```

I also attempted a direct `11:any,23:any` height-`10000` coverage run using the improved exact K3 residue table. It remained too slow as a monolithic coverage check and was interrupted before any chunk completed; all partial files were empty and removed. This reinforces that the useful search direction is boundary-stratum refinement, not an undifferentiated `any` filter.

The practical conclusion after the two completed height-`10000` boundary blocks is that actual surface tuples have strong overlap among boundary components: the non-d/zero family already contains the entire d-collision family at this height. The remaining untested cases are mixed one-prime strata where one prime is only seen through a d-collision component and the other through a non-d/zero component; these should be tested as a smaller targeted complement only if a full coverage certificate at height `10000` is needed.


## Mixed boundary search results

I then tested the remaining mixed complement at height `10000`:

```text
D = {Ead,Ebd,Ecd},
N = {Za,Zb,Zc,Zd,Eab,Eac,Ebc},

11 in D, 23 in N,
and
11 in N, 23 in D.
```

A grouped `11:D,23:N` / `11:N,23:D` run was too dense per chunk to be useful at height `10000`, so I ran the same complement as the 42 sparse component-pair strata. This completed and merged to

```text
B=10000, mixed component-pair family:
  118 tuples.
  This file is byte-identical to the B=10000 non-d/zero family file.
  New relative to non-d/zero family: 0 tuples.
  New relative to d-collision family: 4 tuples.
```

Therefore the two completed height-`10000` searches now give a stable 118-tuple forced-boundary set across the d-collision, non-d/zero, and mixed decompositions that were computationally tractable. Since the mixed file is identical to the non-d/zero file, the previous exact sieve result applies: 0 intermediate and 0 full M(2,2,4,8) witnesses.


## Boundary squareclass mechanism for A(2,2,4,4)

I analyzed whether the boundary strata themselves force failure of the A(2,2,4,4) square conditions. The diagnostic is now reproducible with

```text
code/analyze_a2244_boundary_strata.py
```

and the current output is

```text
data/a2244_boundary_strata_diagnostics_B10000.txt
```

For a fixed 2+2 partition `{i,j}|{k,l}`, the A(2,2,4,4) conditions are equivalent over `Q_p` to the four cross-differences

```text
a_i^2-a_k^2, a_i^2-a_l^2, a_j^2-a_k^2, a_j^2-a_l^2
```

having the same `Q_p^*/Q_p^{*2}` squareclass, allowing zero only in the obvious limiting sense. Thus a simple boundary collision does not automatically kill every partition: if the colliding pair lies on the same side of the 2+2 partition, the odd valuation may not appear among the cross-differences. The obstruction is the combination of valuation parity, unit squareclass, and compatibility of the same partition at both obstruction primes.

On the 118 retained height-`10000` forced-boundary tuples, the result is:

```text
Q11 local possible: 28 / 118
Q23 local possible: 27 / 118
same ordering locally possible at both Q11 and Q23: 1 / 118
same ordering locally possible at Q11, Q23, and R: 0 / 118
```

The surviving partition-set distribution is especially telling:

```text
66 tuples: no partition survives at 11 or 23
13 tuples: only 14|23 survives at 11, none at 23
13 tuples: none at 11, only 14|23 survives at 23
7 tuples: none at 11, only 12|34 survives at 23
6 tuples: only 12|34 survives at 11, none at 23
6 tuples: only 13|24 survives at 11, none at 23
```

Only one tuple has a common finite-local partition at both obstruction primes:

```text
[169,192,468,7644]
```

At `11` it has boundary `E23` and only partition `14|23` survives. At `23` all squared coordinates collide, so all three partitions survive. But every common finite-local ordering has two negative A(2,2,4,4) products, so the real place kills it; the same orderings are also obstructed at many auxiliary finite primes.

Conclusion: a bare boundary component label does not, by itself, force failure. However, on the relevant forced-boundary tuples found up to height `10000`, the boundary squareclass data at `11` and `23`, together with the real sign condition, already forces A(2,2,4,4) failure for every tuple. The mechanism is partition incompatibility and squareclass mismatch on the boundary, not merely the absence of enough height search.

## p-adic residue-class refinement on the boundary

I added a finite residue-class analyzer for the boundary squareclass problem:

```text
code/analyze_a2244_padic_residues.py
```

The current combined output is:

```text
data/a2244_padic_residue_report.txt
```

This works directly on the K3 surface

```text
(ab+ac+ad+bc+bd+cd)^2 = 4abcd
```

modulo `p` and then lifts only the mod-`p` A(2,2,4,4)-ambiguous boundary classes to `p^2`.  The all-zero mod-`p` root is excluded by default, because the tuple search uses primitive projective representatives.

For each 2+2 partition, the test is that the four oriented cross-differences `a_i^2-a_j^2` have the same `Q_p^*/Q_p^{*2}` squareclass.  A cross-difference that is `0 mod p^2` is not declared good; it is recorded as a deep collision requiring higher-order analysis.

The mod-`p` counts are:

```text
p=11:
  Fp roots: 1470
  boundary roots: 1230
  good-reduction roots: 240
  singular roots: 70
  good-reduction roots with any surviving A2244 partition: 0
  per partition: killed 1120, ambiguous 310, resolved_ok 40

p=23:
  Fp roots: 12166
  boundary roots: 6886
  good-reduction roots: 5280
  singular roots: 154
  good-reduction roots with any surviving A2244 partition: 0
  per partition: killed 10120, ambiguous 1210, resolved_ok 836
```

Thus the finite-field obstruction is exactly boundary-supported: no nonboundary mod-`p` class has an A(2,2,4,4) partition surviving at either obstruction prime.

The `p^2` refinement of the mod-`p` ambiguous classes gives:

```text
p=11, per partition among ambiguous root-partitions:
  all_lifts_killed: 0
  deep_only: 280
  resolved_and_deep: 30
  has_resolved_ok_lift: 30
  has_deep_lift: 310

p=23, per partition among ambiguous root-partitions:
  all_lifts_killed: 0
  deep_only: 1144
  resolved_and_deep: 66
  has_resolved_ok_lift: 66
  has_deep_lift: 1210
```

So `p^2` does not prove that the relevant boundary strata themselves force A(2,2,4,4) failure.  In fact every mod-`p` ambiguous root-partition has at least one deep `p^2` lift.  The only `p^2` resolved lifts in this analysis come from the full squared-coordinate collision signature

```text
E12+E13+E14+E23+E24+E34,
```

while the simple and zero-component boundary classes mostly remain deep rather than being killed uniformly.

Conclusion: the boundary labels alone are too coarse for a local impossibility proof.  The evidence now points to a combined mechanism: good reduction is impossible at `11` and `23`, but once forced onto the boundary, failure is caused by partition incompatibility across primes, real signs, auxiliary local conditions, or the full M(2,2,4,8) equations.  The next useful finite computation is not a broader height search; it is a normalized higher-order analysis of the deep `p^2` masks, or equivalently a `p^2`/`p^3` local sieve that keeps the same partition simultaneously at `11` and `23` before feeding residue classes back into the height enumerator.

## Sampled `p^3` refinement of the deep boundary charts

I added a second local diagnostic:

```text
code/analyze_a2244_padic_deep_charts.py
```

The current output is:

```text
data/a2244_padic_deep_chart_p3_report.txt
```

This takes the `p^2`-deep A(2,2,4,4) boundary classes from the previous analyzer, compresses them into broad normalized chart keys, and then lifts several `p^3`-control representatives in each chart.  The control representatives record the next K3 lift obstruction and the `p^2` coefficients of the deep cross-differences.  This is deliberately a chart diagnostic, not a full enumeration of every `p^2` residue class at `p^3`.

The broad chart compression is stable across the two obstruction primes:

```text
p=11:
  broad p^2 chart keys: 273
  sampled p^2 representatives: 618
  p^2 deep root-partitions represented: 2,518,230

p=23:
  broad p^2 chart keys: 273
  sampled p^2 representatives: 666
  p^2 deep root-partitions represented: 102,542,418
```

The sampled `p^3` behavior is mixed, not uniformly obstructed.  For example, at `p=23` the sampled `p^3` lift statuses are:

```text
12|34: killed=7,670,454  deep=3,682,898  resolved_ok=983,986
13|24: killed=7,925,961  deep=3,682,898  resolved_ok=728,479
14|23: killed=7,901,627  deep=3,682,898  resolved_ok=752,813
```

At `p=11` the same qualitative pattern appears:

```text
12|34: killed=192,720  deep=188,276  resolved_ok=63,558
13|24: killed=208,692  deep=188,276  resolved_ok=47,586
14|23: killed=210,023  deep=188,276  resolved_ok=46,255
```

Thus a `p^2` deep mask such as `E13+E23`, `E14`, or the full cross mask does not have a single outcome at `p^3`.  The same broad mask can contain representatives that do not lift to `p^3`, representatives that lift but are killed by squareclass mismatch, representatives that resolve to a valid local A(2,2,4,4) squareclass, and representatives that remain deep at order `p^3`.

The recurring top charts have a simple interpretation.  The largest charts are the three-zero/one-unit boundary charts, e.g.

```text
(a,b,c,d) = (pA,pB,pC,D),  D a unit.
```

For these, the K3 equation at the next order already depends on the linear condition

```text
A+B+C = 0 mod p
```

(up to permutation).  That condition is invisible to the coarser boundary label and explains why the same broad chart contains both no-lift and lifting subcharts.  Once it lifts, the A(2,2,4,4) squareclass outcome depends on the `p^2` coefficients of the deep cross-differences.

Conclusion: the boundary-stratum mechanism is now clearer.  Good reduction is impossible at `11` and `23`, so any rational tuple is forced onto the boundary.  But the boundary itself does not force failure at order `p^2` or in the sampled `p^3` controls.  The obstruction must use finer normalized chart data, and probably the same partition simultaneously at both `11` and `23`, plus the real sign/full M(2,2,4,8) conditions.

The next computationally meaningful step is therefore to build an exact local sieve on the normalized `p^3`-control signatures, rather than pushing all `p^2` residues blindly.  For a height search, this means filtering candidates by the finite data

```text
partition, mod-p boundary signature, p^2 deep mask,
K3 p^3 lift condition, deep-cross p^2 coefficient squareclasses
```

at both `11` and `23` before running the exact A(2,2,4,4)/M(2,2,4,8) sieve.

## Candidate-level signed p-adic signature sieve

I added the candidate-level finite local sieve:

```text
code/a2244_padic_signature_sieve.py
```

The main height-`10000` outputs are:

```text
data/a2244_padic_signature_sieve_B10000.txt
data/a2244_padic_signature_sieve_B10000_p3_survivors.txt
data/a2244_padic_signature_sieve_B10000_depth6.txt
data/a2244_padic_signature_sieve_B10000_depth6_survivors.txt
```

This sieve uses the oriented cross-differences for each 2+2 partition and keeps the sign of the `Q_p` unit.  This matters at `p=11` and `p=23`, where `-1` is nonsquare.  I also regenerated

```text
data/a2244_boundary_strata_diagnostics_B10000.txt
```

with this signed convention.  The corrected exact signed local survivor is now

```text
[605,845,2925,6877]
```

with common finite partition `13|24`; it is not real-compatible and is obstructed at many auxiliary primes.

On the 118 height-`10000` forced-boundary tuples, the signed finite-depth sieve gives:

```text
k=1: finite common at 11 and 23 = 28, real-compatible = 9
k=2: finite common at 11 and 23 = 21, real-compatible = 5
k=3: finite common at 11 and 23 = 14, real-compatible = 4
k=4: finite common at 11 and 23 = 1,  real-compatible = 0
k=5: finite common at 11 and 23 = 1,  real-compatible = 0
k=6: finite common at 11 and 23 = 1,  real-compatible = 0
```

Thus the requested `p^3` signature sieve is useful but not terminal: it reduces the 118 boundary tuples to 14 finite-local candidates, 4 of which are still real-compatible at this finite precision.  One more digit eliminates all real-compatible candidates and leaves only the exact signed finite-local survivor above.

I then fed the 14 `p^3` survivors to the exact all-permutation M(2,2,4,8) sieve:

```text
data/m2248_p3_signature_survivors_intermediate_B10000.txt
data/m2248_p3_signature_survivors_full_B10000.txt
```

Both files are empty:

```text
intermediate witnesses: 0
full witnesses: 0
```

Practical conclusion: for future height search, the useful prefilter is the signed finite-depth partition sieve at `11` and `23`.  At height `10000`, `p^3` already cuts the exact M-sieve input from 118 to 14, and depth `4` plus the real sign condition cuts it to 0 before the M(2,2,4,8) stage.  The K3 next-lift flag should only be used when a signed surface representative is retained; the positive sorted curve tuple does not necessarily satisfy the symmetric K3 equation with all positive signs.

## Integrated depth-filtered height-search pilot

I extended the C++ enumerator so it can apply the signed finite-depth real-compatible local filter at tuple creation.  The optional final argument is now:

```text
REAL_LOCAL_DEPTH
```

For the sorted positive curve tuples, the only real-compatible A(2,2,4,4) partition is `12|34`, so the integrated filter keeps a tuple only if that partition is still locally possible at every required obstruction prime to the requested depth.  The intended setting is therefore:

```text
required primes: 11,23
REAL_LOCAL_DEPTH: 4
```

I also added a resumable chunk driver:

```text
code/run_surface_tuple_chunks.py
```

It runs `code/enumerate_surface_tuples` on consecutive `a`-ranges, writes successful chunks atomically via temporary files, and merges completed chunks.  This is the right way to push beyond height `10000`; monolithic runs give no progress information and are too awkward to interrupt.

A first `B=20000` pilot was run on the non-d/zero boundary block with the depth-4 real-compatible filter:

```text
B=20000
boundary=11:N,23:N
REAL_LOCAL_DEPTH=4
ax=1..10
```

Output:

```text
data/surface_tuples_B20000_N_depth4_real_11_23_ax1_50.txt
```

The completed chunk was:

```text
ax=1..10: 0 tuples, elapsed 115.957s
```

This is useful benchmark data: even with the strong local filter, the cost is dominated by enumerating the K3 triples before candidate tuples are formed.  The filter is excellent for keeping output and Magma input tiny, but it does not by itself make a full `B=20000` sweep cheap.  A full higher-height search should therefore be run as resumable chunks, preferably after either improving the residue prefilter to incorporate the real-compatible partition condition earlier or distributing the `a`-range chunks.

## Early common-partition residue prefilter

I improved the C++ enumerator again by moving a necessary A(2,2,4,4) condition into the residue prefilter.  Previously, the prefilter allowed signed residue triples `(a,b,c) mod 11*23` whenever some `d mod p` put the K3 point on the requested boundary components.  The new version also requires, already modulo each obstruction prime, that the same raw 2+2 partition can survive the oriented A(2,2,4,4) cross-difference squareclass test.

This is a cheap mod-`p` necessary condition, not the full depth-4 test.  It is safe as a prefilter because any eventual A(2,2,4,4) tuple has one global 2+2 partition that must be locally possible at both `11` and `23`.  The stronger depth-4 real-compatible `12|34` test still runs on the final sorted positive tuple.

For the `11:N,23:N` boundary block, the residue-table reduction is:

```text
without early common-partition prefilter:
  allowed signed residue triples modulo 253: 3,133,151 / 16,194,277
  allowed positive c residue entries:        6,266,302 / 32,388,554

with early common-partition prefilter:
  allowed signed residue triples modulo 253:   964,463 / 16,194,277
  allowed positive c residue entries:        1,928,926 / 32,388,554
```

So the residue loop is reduced by about a factor of `3.25` before the expensive integer search begins.

I reran the same `B=20000`, `ax=1..10`, depth-4 real-compatible pilot:

```text
before this prefilter: 115.957s, 0 tuples
after this prefilter:   35.7s,   0 tuples
```

The new pilot output is:

```text
data/surface_tuples_B20000_N_depth4_real_earlypart_11_23_ax1_10.txt
```

This makes a chunked `B=20000` sweep substantially more realistic.  It is still not instant, but the bottleneck has moved in the right direction: the search spends far less time on residue classes that could never support a common A(2,2,4,4) partition at the obstruction primes.


## Component-wise 11/23-adic boundary analysis

I added a component-wise finite local analyzer:

```text
code/a2244_component_adic_analysis.py
```

The report is:

```text
data/a2244_component_adic_analysis.txt
```

It separates the ten raw K3 boundary components

```text
Z1, Z2, Z3, Z4, E12, E13, E14, E23, E24, E34
```

and the three A(2,2,4,4) partitions.  For each of `p=11` and `p=23`, it enumerates all mod-`p` K3 roots, exhaustively lifts every mod-`p` ambiguous component/partition class to `p^2`, and then searches for smooth `p^3` resolved witnesses.

The main outcome is negative for a component-wise obstruction:

```text
p=11:
  Fp_roots=1470
  boundary_roots=1230
  smooth_boundary_roots=1160

p=23:
  Fp_roots=12166
  boundary_roots=6886
  smooth_boundary_roots=6732
```

For both primes, every one of the 30 component/partition pairs has verdict

```text
smooth_p3_resolved
```

in the report.  Thus no raw component `Zi` or `Eij`, even together with a chosen A(2,2,4,4) partition, can be ruled out locally at `11` or at `23` by this component-wise analysis.  In particular, the forced bad-reduction condition is real, but the boundary components themselves are not the obstruction.

The interpretation is now sharper: any obstruction to rational examples must use data not visible from one component at one prime.  The relevant mechanisms are the same global partition having to survive simultaneously at `11` and `23`, real compatibility of that partition, and finally the full `M(2,2,4,8)` cover equations.  This matches the height-`10000` signed-depth sieve, where depth 4 leaves one finite-local tuple but no real-compatible tuple, and the exact M-sieve gives no intermediate or full witnesses.
