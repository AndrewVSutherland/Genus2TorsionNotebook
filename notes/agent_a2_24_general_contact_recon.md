# General non-d=0 cubic-contact reconnaissance for [2,24]

Date: 2026-07-06.

Purpose: follow the open door left by the `d=0` rank-0 analysis.  The
contact identity is

```text
f = h3^2 + kappa*(x^2 + U*x + V)^3.
```

The `d=0` split-through-zero slice is now a clean negative, but curve A has
nonzero `d` and irreducible `q3`, so the general contact cover still had to
be tested directly.

## Files

New:

```text
code/agent_a2_24_general_contact_recon.m
results/a2_24_general_contact_recon_p7.log
results/a2_24_general_contact_recon_p11_13_17.log
```

Also patched:

```text
code/agent_a2_24_locus_geometry.m
results/a2_24_locus_geometry_p7_11_13_corrected.log
```

The patch fixes the finite-field D8 representative.  In finite fields the
unscaled A(8) model should use

```text
v8 = -ellBase mod g8
```

not the rational-search integer-model scaling.  After this correction the
finite-field geometry counts are nonzero.

## Exact Q anchors

The script validates the two known simple cyclic Z/24 curves.

Curve A:

```text
r=5, p=-5/2, t=-9/2
d=-57/4
torsion=[24]
q3=x^2 - 435/73*x + 2529/292
disc(q3)=4608/5329, nonsquare over Q
contact identity true
```

Curve B:

```text
r=1/3, p=-1/9, t=-1
d=0
torsion=[24]
q3=x^2 - 2/3*x
disc(q3)=4/9, split over Q
contact identity true
```

So the non-`d=0`, irreducible-`q3` contact branch is not merely formal: it
contains the known simple cyclic Z/24 anchor A.

## Finite-field W-split contact survey

The survey loops through A(8) over `F_ell`, keeps smooth non-even curves,
checks the visible order-8 class, keeps W-split / 2-rank at least 2 curves,
then enumerates actual Mumford representatives `[u,v]` with

```text
u = x^2 + U*x + V,     u | f-v^2,     3*[u,v] = 0.
```

The count `contact_classes` is deduplicated by the contact conic `u=q3`, so
it counts conics rather than the signed pair `D,-D`.

```text
ell  smooth  wsplit  wsplit_j3  contact_bases  contact_conics  split  irreducible  non_d0_bases  non_d0_irred_bases  d0_bases
7    67      18      3          1              1               1      0            0             0                   1
11   472     175     67         50             50              13     37           36            29                  14
13   914     376     92         69             78              29     49           44            33                  25
17   2538    1150    385        331            337             172    165          271           135                 60
```

Examples from the non-`d=0`, irreducible branch:

```text
F_11: <r,p,t,d,2rank,q3,split?> =
      <2,1,3,4,2, x^2 + 3*x + 10, false>

F_13: <2,1,9,5,2, x^2 + 4*x + 2, false>

F_17: <2,3,3,2,2, x^2 + 7*x + 1, false>
```

## Lightweight corrected geometry count

After fixing the finite-field D8 representative, the older geometry helper
reports:

```text
ell  N_o8  N_o8_r2  N_full
7    67    18       3
11   472   175      67
13   914   376      92
```

Here `N_full` uses the necessary condition `3 | #J(F_ell)`, while the new
reconnaissance script enumerates the actual rational contact conics.

## Interpretation

The general non-`d=0` contact door is real.

The d=0 slice is special: at `ell=7`, the only W-split contact point found
is d=0 and split.  But by `ell=11,13,17`, the W-split locus has many contact
points with `d != 0`, and many have irreducible `q3`.

This means the negative d=0 result should not be read as a local obstruction
to [2,24].  It is only a negative for the friendliest slice.

The harder lesson is geometric: rational 3-torsion is a finite cover, not a
hypersurface condition in the A(8) base.  Intersecting the cover with W-split
does not give an obvious rational curve in `(r,p,t)`.  It gives a cover of
the W-split surface, and rational points on that cover appear globally thin.

## Next move

The next serious build should not be another local-count probe.  It should
construct the W-split contact cover with variables

```text
(r,t,beta; U,V,h0,h1,h2,h3,kappa)
```

where `beta` parametrizes the W-split condition, then eliminate or specialize
to find low-dimensional rational subcovers.  The known Q anchor A should be
used as the non-d0 irreducible-contact sanity check throughout.

Practical first specialization: fix or parametrize a small invariant of
`q3`, such as its discriminant squareclass, and see whether the resulting
W-split contact cover drops to genus 0 or genus 1 curves that can be searched
or descended.
