# Four-task torsion continuation, 2026-07-10

This note records the four follow-ups to
`notes/six_task_followup_2026_07_10.md`. Each calculation was delegated,
bounded, recorded, and then independently reviewed. No new rational torsion
example was found in this pass.

## 1. Degree-2 `[5,5]` Hensel reconstruction

Files:

```text
code/z5x5_degree2_fiber_branches.sage
code/z5x5_degree2_hensel_reconstruct.sage
data/z5x5_degree2_fiber_branches_p11_p19.txt
data/z5x5_degree2_hensel_h30_p7_all_branches.txt
notes/z5x5_degree2_hensel_reconstruction.md
```

The fixed base fibers contain 100 smooth open branches modulo 11 and 16
modulo 19. The search tested all 1,600 cross-prime branch pairs above all 120
rational base triples of coordinate height at most 30, through precisions
`11^7` and `19^7`: 192,000 lifted pairs in total. Of 280 fully reconstructed
fiber tuples, none satisfied the ten equations over `QQ`.

This stops only these two base residue cells. The construction remains live;
the next search must vary the base residue triples rather than lift these cells
more deeply.

## 2. Nonautomorphic `[2,6,6]` surface

Files:

```text
code/contact6_m36_266_surface_decompose.sage
data/contact6_m36_266_surface_decompose_p13_r4_full_open.txt
data/contact6_m36_266_surface_decompose_p19_r4_full_open.txt
notes/contact6_m36_266_surface_decomposition.md
```

On the tractable `r=4` slice, structural saturation removes a two-dimensional
boundary component and leaves a curve with Hilbert polynomial `95*t-355` and
total affine degree 95 over both `F_13` and `F_19`. Smoothness, coprimality,
and both automorphism saturations leave this curve unchanged. Thus a genuine
nonautomorphic curve closure survives; it is not an isolated finite-field
point or a detachable automorphism component.

Primary decomposition, the full five-variable saturation, and exact-`QQ`
saturation reached their 240-second limits. No individual component degree,
normalization, or rational point is claimed.

## 3. Degree-16 `[64]` cover

Files:

```text
code/task3_cover_different_bound.sage
data/task3_cover_different_bound.txt
notes/task3_cover_different_bound.md
```

The bad-fiber support has 124 geometric points. Squarefree residue
certificates prove that the 112 points over the irreducible degree-8 and
degree-104 factors are unramified in the degree-16 cover. The remaining 12
points contribute at most 15 each in characteristic zero, so

```text
deg Different <= 180,
g <= 75.
```

This is the unconditional improvement from the old Newton bound `g<=1785`.
Local Puiseux experiments predict `g=33`, with an intermediate `g<=57`, but
those values remain exploratory until the local equisingularity or integral
bases are certified over characteristic zero.

## 4. Boundary audit for `[3,18]` and `[6,12]`

Files:

```text
code/task4_boundary_component_audit.sage
data/task4_boundary_component_audit.txt
notes/task4_boundary_component_audit.md
```

For contact 9, put `t=eps*s`. The ordinary-node charts `t=3,5 mod 7` have
insufficient 3-primary order and are locally impossible. The chart
`t=4 mod 7` has elliptic normalization `[9]`, a split torus, and the required
generalized-Jacobian 3-primary rank; it is the primary live `[3,18]` chart.
The pole `t=-1 mod 7` also passes the necessary capacity test and is secondary.
These are necessary local tests, not rational examples.

For `[6,12]`, all 42 rank-five residues and their 1,050 lifts modulo 25 were
enumerated. None is fully transverse. The only 160 curve-smoothing lifts were
followed through all 4,000 continuations modulo 125, and every one stays on
the two torsion-degeneracy divisors. The current affine boundary route is
therefore stopped.

## Priority after these four tasks

1. Vary the `[5,5]` base residue cells using small rational bases, two-prime
   25-divisibility screens, and simple-Frobenius screens before branch lifting.
2. Resolve the contact-9 `t=4 mod 7` nodal chart explicitly, including its
   component-group extension, before another `[3,18]` height search.
3. Compute associated primes and normalizations of the degree-95 `r=4`
   `[2,6,6]` curve, matching the components through the known `p=13,19` points.
4. Compute local integral bases at the 12 unresolved `[64]` boundary points
   only if certifying the predicted genus is a priority.
5. Do not spend more compute on the present `[6,12]` chart without a different
   geometric construction.
