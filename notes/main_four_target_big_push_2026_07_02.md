# Four target big push ledger

Date: 2026-07-02.

This ledger tracks the large-budget campaign following
`main_four_target_third_pass_2026_07_02.md`.

## Workers

```text
Zeno      A(2,24) quartic saturated halving branch extraction
Halley    Z/35 b=0,r=1 liftable pole directions beyond mod 27
Rawls     Z/48 RTHeight=4 outside RTHeight=3, SliceClass 1
Newton    Z/48 RTHeight=4 outside RTHeight=3, SliceClass 2
Nash      Z/48 RTHeight=4 outside RTHeight=3, SliceClass 3
Plato     Z/5 x Z/5 full norm, b2=0 elimination/modular saturation
```

## Local integration

Added:

```text
code/torsion_cover_lab_utils.m
```

The helper file is intentionally small: height-rational enumeration,
factor-degree summaries, good-reduction conversion, point-count torsion gates,
and torsion invariant normalization.  It is meant for new scripts, not as a
forced refactor of the existing successful one-off drivers.

## Priority order

```text
1. A(2,24): explicit quartic branches and rational-point/cover analysis.
2. Z/35: branch lifting and parameterization for the 18 liftable directions.
3. Z/48: complete remaining partitions and read off whether the A16 gate is
   fully cold at RTHeight 4.
4. Z/5 x Z/5: b2=0 branch first, then full norm modular components.
5. Infrastructure: only promote helpers after at least two lanes use them.
```

## Operating Pattern From Previous Hits

The useful past examples in this folder follow a repeated pattern:

```text
visible contact torsion
  -> add one extra torsion condition as a cover equation
  -> saturate away fake boundary components
  -> look for unexpectedly small factors/components
  -> point-count gate before expensive exact Jacobians
  -> exact TorsionSubgroup/divisor certification only at the end
```

For `[2,2,20]`, the decisive move was to factor the residual quartic on a
one-parameter contact-5 family, then parametrize the first extra-2 locus.  For
`[6,6]`, the productive route was not the first broad contact-6 scan, but the
extra-root/core-cover refinement plus finite filters before exact torsion.

This is why the current campaign emphasizes the A(2,24) quartic factors, the
Z/35 pole branches, and the Z/48 cubic-contact route rather than enlarging
undifferentiated height boxes.

## Running status

### A(2,24) Quartic Branch Extraction

Files:

```text
code/agent_A2_24_quartic_extract.m
notes/agent_A2_24_quartic_extract.md
results/A2_24_quartic_extract_default.log
```

The four best split fibers

```text
(-1/3,-1, 4/3),  (-1/3,1, 4/3),
( 1/3,-1,-4/3),  ( 1/3,1,-4/3)
```

all reduce to the same genus-2 curve data.  The two extra rational 2-torsion
translations are

```text
T1 = [x^2 + x + 1/3, 0],
T2 = [x^2 + 8/3*x + 16/3, 0].
```

For each extra translation, after removing the boundary

```text
s4 = M^2 - 9/125,
```

the saturated affine factor degrees are `[4,4,8]`.  The quartic components
collapse to three distinct irreducible quartics:

```text
qA = M^4 + 12/13*M^3 + 1728/4225*M^2
     + 1512/21125*M + 324/105625

qB = M^4 + 12/13*M^3 + 2214/4225*M^2
     - 3348/21125*M - 1863/105625

qC = M^4 - 84/143*M^3 + 56646/511225*M^2
     - 14364/2556125*M - 3807/12780625
```

All three are irreducible over `Q` and have no rational roots.  On each
quartic branch, `N` is forced linearly over `Q(M)`, so the rational point
question on these zero-dimensional branches is exactly the rational-root
question for `qA,qB,qC`.  No exact torsion certification was triggered.

Follow-up spawned: close the remaining O/TR degree-16 pieces, degree-8
residual factors, and the `s4=0` boundary for the same common fiber.

### A(2,24) Four-Fiber Closure

Files:

```text
code/agent_A2_24_branch_closure.m
notes/agent_A2_24_branch_closure.md
results/A2_24_branch_closure_compact.log
results/A2_24_branch_closure_default.log
```

The closure pass checked all 16 translated exact order-12 classes occurring in
the four best fibers.  After removing `s4=0`, every saturated affine
`M`-factor has no rational root:

```text
D_A, D_B       O/TR classes:       [16], [16]
D_minus,D_plus extra translations: [4,4,8], [4,4,8]
```

The residual degree-8 factors and the O/TR degree-16 factors are irreducible
over `Q` in the recorded runs and have no rational roots.  The boundary is

```text
s4 = M^2 - 9/125,
```

so it has no rational point; over `Q(sqrt(5))` it degenerates to the expected
fake degree-3 square-quartic boundary.  Exact integral-model checks also give
`IsDivisibleBy(D,2)=false` for all 16 classes.

Verdict: the four best A(2,12) fibers are closed for A(2,24).  The next A2
move, if desired, is a higher-height split-fiber scan for new low-degree
components, not more work on these four fibers.

### A(2,24) Height-5 Expansion

Files:

```text
code/agent_A2_24_height5_lowbranch_scan.m
notes/agent_A2_24_height5_lowbranch_scan.md
results/A2_24_height5_full.log
```

The full height-5 split-fiber scan, skipping the four already closed fibers,
found:

```text
checked=56316
split_fibers=84
translated_order12_rows=352
[16]  : 280 rows
[8,8] : 72 rows
low_rows_le_4=0
rational_M_root_rows=0
errors=0
```

A height-4 sanity run with the closed fibers included confirms that the only
degree-4 rows in the height-5 box are exactly the already closed four fibers,
two extra translations each.  Thus the height-5 shell adds no new A(2,24)
candidate and no exact divisibility target.

Verdict: deprioritize A(2,24) for now unless we want a much broader,
partitioned height search or a different A(2,12) chart.

### Z/5 x Z/5 Linear-B Full Norm Branch

Files:

```text
code/agent_z5x5_b2zero_elim.m
notes/agent_z5x5_b2zero_elim.md
results/z5x5_b2zero_count_p7.log
results/z5x5_b2zero_slice_sat_p7.log
results/z5x5_b2zero_count_p11.log
```

On the `b2=0` branch, the residual equations drop from degrees
`30,27,24,21,18` to degrees `10,9,8,7,6`.  The open product used for the
bounded algebra was:

```text
K * b1 * b0 * disc(U) * Res(B,U) * disc(f).
```

Finite-field counts:

```text
F_3:  raw residual points 103,   open points 0
F_7:  raw residual points 2521,  open points 12
F_11: raw residual points 15281, open points 240
```

The `F_7` and `F_11` sampled open points have Jacobian rank `5` for the five
residual equations, and Magma verifies the contact class and `[U,V]` are
independent 5-torsion classes.  The best local chart found is the `F_7` slice

```text
h1=1, h2=0,
```

where raw slice dimension is `2`, and saturation by the open product gives a
zero-dimensional chart with `11` primary components and two open rational
points over `F_7`.

Verdict: `b2=0` is not killed by the open-chart cleanup.  It is now a finite
field lift / rational reconstruction target.

### Z/5 x Z/5 Linear-B Lift Attempt

Files:

```text
code/agent_z5x5_b2zero_lift.m
notes/agent_z5x5_b2zero_lift.md
results/z5x5_b2zero_lift_default.log
results/z5x5_b2zero_lift_wide.log
```

On the `h1=1,h2=0` slice, the two open `F_7` points and one comparison `F_11`
point are nonsingular in the five slice variables `(K,s,t,b0,b1)`:

```text
F7_slice_A: rank 5, determinant 4 mod 7, unique lift through 7^4
F7_slice_B: rank 5, determinant 4 mod 7, unique lift through 7^4
F11_slice_A: rank 5, determinant 4 mod 11, unique lift through 11^3
```

The wider rational reconstruction search used

```text
|num| <= 200, den <= 100
```

constrained to the final p-adic residue classes.  It found no exact rational
tuple:

```text
F7_slice_A:  90000 tuples tested, 0 exact hits
F7_slice_B:  90000 tuples tested, 0 exact hits
F11_slice_A: 1591812 tuples tested, 0 exact hits
CRT(F7,F11): no coordinate candidates in the same box
```

Verdict: these are genuine smooth p-adic points of the saturated open chart,
not local obstructions.  A rational point, if present in these charts, is
outside the bounded reconstruction box or must be approached by a stronger
global method.

### Z/35 Pole-Branch Lift

Files:

```text
code/agent_Z35_liftable_branch_lift.m
notes/agent_Z35_liftable_branch_lift.md
results/Z35_lift_branch_k7.log
```

The previous mod-27 viability table was too optimistic.  Pushing the 18
liftable first directions for the `b=0,r=1` pole chart gives:

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

Branch taxonomy:

```text
12/18 die at mod 81.
4/18 lift through mod 729 but die at mod 2187.
2/18 survive through mod 2187.
```

The two survivors are the central first directions:

```text
t=1: <1,1,1,0,1,0>
t=2: <1,1,1,2,1,0>
```

No height-80 rational reconstruction produced an exact solution.  Stored
samples verify the original coefficient equations at their recorded precision
and have nonzero discriminant, so the remaining obstruction is a higher
3-adic issue rather than immediate singularity.

Verdict: Z/35 is no longer an 18-branch problem.  It is now a two-central-branch
problem.

### Z/35 Central-Branch Deep Lift

Files:

```text
code/agent_Z35_central_branch_deep_lift.m
notes/agent_Z35_central_branch_deep_lift.md
results/Z35_central_smoke_k4.log
results/Z35_central_deep_k8.log
results/Z35_central_deep_k10.log
```

For both central directions the mod-`3` Jacobian has rank `3`, so nonempty
next-digit sets are affine three-dimensional cosets.  Common left-obstruction
rows on the five scaled residuals are:

```text
[1,0,2,1,0],  [0,1,2,0,1].
```

The exact smoke transition through scaled `3^4` gives, for the two central
branches together:

```text
k  modulus  total lifts  nonzero central branches
1  3        2            2
2  9        54           2
3  27       1458         2
4  81       13122        2
```

For the `t=1` central branch, direct enumeration and one-step lookahead in the
deep run record exact rows through scaled `3^7`:

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

Representative diagnostics for `3^8..3^10` show sampled all-or-none child
coset behavior, but those rows are not certified full counts.  No rational
reconstruction candidate appeared in the recorded probes.

Verdict: raw enumeration is now the wrong tool.  The next Z/35 step should be
a compressed-state automaton or symbolic obstruction calculation keyed by the
two left-obstruction residuals and the changing affine constants.

### Z/48 RTHeight 4 Partition Sweep

Files:

```text
notes/agent_Z48_next_scan.md
notes/agent_Z48_next_scan_slice1.md
notes/agent_Z48_next_scan_slice2.md
notes/agent_Z48_next_scan_slice3.md
results/Z48_next_slice1.log
results/Z48_next_slice2.log
results/Z48_next_slice3.log
```

The remaining partitions `SliceClass:=1,2,3` all completed successfully.
Together with the earlier `SliceClass:=0` run, this finishes the reproducible
RTHeight 4 search outside the completed RTHeight 3 box.

Aggregate over `SliceClass:=0,1,2,3`:

```text
runSlices=280
tested=4516120
commonRootPairs=32
rationalRoots=34
singular=25
nonsingular=9
pointGateReject=8
pointGatePass=1
exactTried=1
certified=1
z48Hits=0
```

The only point-count gate survivor was in `SliceClass:=0`; exact A16
certification succeeded, but the integral square model had torsion invariants
`[16]`, not `[48]`.  The other eight nonsingular roots were killed by the
finite-field gate:

```text
p=5 : 4
p=7 : 4
```

Verdict: the bounded A16 square-root route is cold through RTHeight 4 outside
RTHeight 3.  The next genuinely new Z/48 route should add the cubic-contact
equations for the 3-part, not just keep enlarging the same A8/A16 prefilter.

### Z/48 Cubic-Contact / Exact 3-Part Route

Files:

```text
code/agent_Z48_cubic_contact_route.m
notes/agent_Z48_cubic_contact_route.md
results/Z48_cubic_contact_smoke.log
results/Z48_cubic_contact_small_search.log
results/Z48_cubic_contact_fixed_slice.log
```

The new diagnostic sends A16 point-count survivors through an explicit
3-primary layer.  Algebraically it looks for cubic-contact witnesses

```text
h(x)^2 - f(x) = Lambda*q(x)^3,
q = x^2 + U*x + V,
h = H3*x^3 + H2*x^2 + H1*x + H0,
```

and decisively checks exact rational 3-torsion via `TorsionSubgroup` on an
integral square model.

The seven smoke A16 candidates all verified as A16 and all were rejected by
the exact 3-part test:

```text
SMOKE_DONE candidates=7 verified16=7 rejectedByExact3=7
contactHeightHits=0 z48Unexpected=0
```

Most importantly, the old RTHeight 4 `SliceClass:=0` point-count survivor

```text
r=-1/4, t=-1/4, mu=-1/2, y=-5/8, p=-41/144, N=5/8, z=125/96
```

still passes the old `48 | #J(F_l)` gate with running gcd `48`, but the new
diagnostic gives:

```text
torsion=[16], has3=false, contactHit=false.
```

The fixed-slice rerun finds exactly this survivor and kills it before any
Z/48 certification.  Verdict: the cubic-contact/exact-3 layer is now the right
filter for future Z/48 production searches.
