# Six-task torsion follow-up, 2026-07-10

This note consolidates the six bounded follow-ups launched from the repository
audit. Each task has its own reproducible code, data, and detailed note. The
notation `[n1,...,nr]` means the corresponding product of cyclic groups.

## 1. Package the exact `[6,6]` example

Files:

```text
code/contact6_m36_66_package.sage
code/contact6_m36_66_package.m
data/contact6_m36_66_package_sage.txt
notes/contact6_m36_66_example.md
```

The reduced model

```text
y^2 = 1872*x^5 - 3000*x^4 + 6969*x^3 - 1691*x^2 + 4875*x
```

has two explicit order-6 generators spanning 36 points. Good reductions give
`#J(F_7)=36` and `#J(F_11)=144`, proving exact torsion `[6,6]`. Sage/Lombardo
returns geometric endomorphism ring `ZZ`, and the 12th-power Frobenius transform
at `p=37` is irreducible. This task is complete.

## 2. Target a non-automorphic `[2,6,6]` lift

Files:

```text
code/contact6_m36_266_targeted_lift_sage.py
code/contact6_m36_266_targeted_lift.m
data/contact6_m36_266_targeted_lift.txt
notes/contact6_m36_266_targeted_lift.md
```

The old `p=19` residue `(r,b)=(2,3)` is on both explicit automorphism loci,
despite its irreducible finite-field Frobenius quartic. A genuinely
non-automorphic residue `(r,b)=(4,5)` has smooth roots at `p=13,19`; all lift
uniquely through `p^6`. The four CRT pairings, with balanced reconstruction
bound `10655549`, give no exact rational point. The next justified step is
full-surface saturation and modular component decomposition, not deeper lifting
of these fixed fibers.

## 3. Scout the contact-9 route to `[3,18]`

Files:

```text
code/contact9_318_finite_scout.sage
code/contact9_318_finite_scout.m
data/contact9_318_finite_p7.txt
notes/contact9_318_finite_scout.md
```

At `p=7`, the only two good affine root-chart residues have
`J(F_7)=[2,36]`, which does not contain `[3,18]`. Every other affine residue is
singular or has `r=0`. The generic affine route is stopped; only a separate
7-adic boundary project could reopen it.

## 4. Test `[6,12]` on the `[6,6]` core

Files:

```text
code/contact6_m36_612_local_feasibility.sage
data/contact6_m36_612_local_feasibility_p5_p11.txt
notes/contact6_m36_612_local_feasibility.md
```

All 476 intersection points modulo 5 are degenerate or have bad special fiber,
so the affine open has no `Z_5` point. Modulo 11 there are 56 fully open target
points with invariants `[12,12]`, but all have Frobenius `(T^2+11)^2` and are
split supersingular. This is a no-go for a broad affine search, not a theorem
excluding points on blown-up 5-adic boundary charts.

## 5. Triage the `[8,8]` and `[64]` covers

Files:

```text
code/task5_cover_geometry.sage
data/task5_cover_geometry.txt
notes/task5_cover_geometry.md
```

The reconstructed `[32]` base has the clean parameterization

```text
r = 1/(t^4-1),
z = 2*(t^4+t^3+t^2+t+1)/(t^2*(t^2+t+1)).
```

The corrected `[64]` cover is an irreducible degree-16 extension of `Q(t)`.
Its exact genus remains uncomputed; the integral-basis and discriminant steps
hit their time caps. The reduced `[8,8]` target is a surface, not a curve. Its
compact second-stage eliminant has degree 8, and the generic cover degree is
128 over `(R,s)`. The next useful work is a controlled one-coordinate fiber or
function-field different calculation, not a larger height search.

## 6. Replace the failed point-contact `[5,5]` model

Files:

```text
code/agent_z5x5_degree2_contact_probe.sage
data/agent_z5x5_degree2_contact_probe_p11_p19.txt
notes/agent_z5x5_degree2_contact.md
```

For `f=(1+a*x+b*x^2)^2-k*x^5`, a general degree-2 order-5 class is governed by

```text
H^2 - f*R^2 = q^5,
deg(H)=5, deg(R)<=2, deg(q)=2.
```

The ten coefficient equations in thirteen variables have smooth, nonboundary,
expected-dimensional points at both `p=11` and `p=19`. At both points the two
order-5 classes are independent, the fiber-variable Jacobian is invertible, and
the 12th-power Frobenius transform is irreducible. This is the strongest live
new-construction route. The next bounded step is Hensel lifting followed by
rational reconstruction in `(a,b,k)`.

## Priority after the six tasks

1. Hensel-lift and reconstruct the degree-2 `[5,5]` cover.
2. Decompose the non-automorphic `[2,6,6]` surface after saturation.
3. Compute the different/genus of the degree-16 `[64]` function field or take a
   controlled `[8,8]` fiber.
4. Leave `[3,18]` and affine `[6,12]` stopped unless boundary geometry becomes
   a separate objective.

Magma was not installed on this system during this run. All new
executed computations therefore use Sage/Singular. The packaged `[6,6]` result
has an independent exact Sage proof and agrees with the repository's existing
captured Magma output. Any future rational hit should still receive a final
exact Magma torsion computation on a licensed host.

## Continuation

The four priorities above were subsequently executed. Their bounded results
and revised priority order are recorded in
`notes/four_task_followup_2026_07_10.md`.
