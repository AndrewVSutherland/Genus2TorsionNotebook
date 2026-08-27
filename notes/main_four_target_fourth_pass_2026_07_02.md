# Four target fourth pass

Date: 2026-07-02.

This summarizes the large-budget campaign recorded in
`main_four_target_big_push_2026_07_02.md`.

## Executive summary

```text
Z/2 x Z/24:
  The four best A(2,12) fibers are closed.  Height 5 gives no new low-degree
  branch or rational M-root.  Deprioritize this lane for now.

Z/35:
  The b=0,r=1 pole chart narrowed from 18 first directions to two central
  3-adic branches.  They survive through scaled 3^7, but raw enumeration past
  that is not viable.  Next move: compressed obstruction automaton.

Z/48:
  RTHeight 4 outside RTHeight 3 is complete and cold.  The only point-count
  survivor is exact torsion [16].  A new cubic-contact/exact-3 filter rejects
  it and is the right route for future Z/48 runs.

Z/5 x Z/5:
  The b2=0 full-norm branch is live.  Open F_7/F_11 points are smooth
  p-adic points with unique lifts, but bounded rational reconstruction found
  no Q-point.  Next move: global algebra or stronger local/global sieve.
```

## Best next directions

1. `Z/35`: build the compressed-state automaton for the two central branches.
   The exact data shows a rank-3 affine correction coset plus two
   left-obstruction residuals.  Proving the sampled all-or-none child-coset
   rule is the most natural way to push beyond `3^7`.

2. `Z/5 x Z/5`: use the smooth `b2=0` p-adic charts as local conditions in a
   global search or elimination.  The finite-field points are real; the
   rational obstruction, if any, is not a naive local obstruction.

3. `Z/48`: keep the new cubic-contact/exact-3 filter, but avoid rerunning the
   same A16 boxes.  Future production should either enlarge `SearchBound`
   deliberately or use fixed-slice reruns when new point-count survivors appear.

4. `Z/2 x Z/24`: pause the current split-fiber approach unless moving to a
   much broader partitioned height search or a different A(2,12) chart.  The
   four degree-4 fibers and the height-5 shell are closed for rational halving.

## Main artifacts

```text
notes/main_four_target_big_push_2026_07_02.md
code/torsion_cover_lab_utils.m

notes/agent_A2_24_quartic_extract.md
notes/agent_A2_24_branch_closure.md
notes/agent_A2_24_height5_lowbranch_scan.md

notes/agent_Z35_liftable_branch_lift.md
notes/agent_Z35_central_branch_deep_lift.md

notes/agent_Z48_next_scan_slice1.md
notes/agent_Z48_next_scan_slice2.md
notes/agent_Z48_next_scan_slice3.md
notes/agent_Z48_cubic_contact_route.md

notes/agent_z5x5_b2zero_elim.md
notes/agent_z5x5_b2zero_lift.md
```
