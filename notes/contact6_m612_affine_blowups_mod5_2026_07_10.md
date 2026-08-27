# Contact-6 `[6,12]`: corrected affine first blow-ups modulo 5

Date: 2026-07-10.

## Scope and corrections

This audits

```text
code/contact6_m612_boundary_first_blowups_mod5.m
```

for the five affine base supports

```text
b+3, a+3, a+b+2, DB, DC
```

paired with the four deeper contact boundaries

```text
L=0, v=0, U-2v=0, U+2v=0.
```

The original output was not a valid smooth-branch count for two reasons.

1. Every chart was placed in `F_5[q,t,z,L,U,v]` even though one of
   `L,U,v` had been replaced by `q` or `q*z`.  That unused coordinate
   multiplied every count by 5.
2. Solvability of one inhomogeneous first-order equation was called
   "transverse" even when the fiber Jacobian had rank only 1 or 2.  Such a
   point is not Hensel-smooth.

The corrected code uses a genuine five-variable ring

```text
F_5[q,t,z,X,Y]
```

for each chart and deduplicates the reciprocal overlap

```text
[e:D]=[z:1]=[1:z^(-1)].
```

Only `z != 0` is counted here.  These are the balanced first-blowup cones,
where the base normal `e` and the deeper contact parameter `D` have the same
leading order.

## Strict-transform point that matters

The proper transform of the core ideal is

```text
Saturation(<F1,F2,F3>,<q>).
```

Dividing a power of `q` from each displayed generator separately is not a
replacement for this saturation: combinations of the generators can add
strict-transform equations.

A full symbolic saturation was attempted.  Even the irrelevant
`(b+3,L)`, chart-1 computation produced a Groebner basis with 535 elements;
nontrivial charts did not finish in the allotted run and were stopped.
Accordingly, the code now uses a smaller exact local certificate and does
not label rank-deficient points as strict-transform branches.

At a raw exceptional point let `J_fib` be the `3 x 4` Jacobian in
`(t,z,X,Y)`, omitting the `q` column.  If

```text
rank(J_fib)=3,
```

then the raw pullback is smooth and the implicit-function theorem makes `q`
a free local parameter.  Its component is not contained in `q=0`, so
saturation by `q` changes nothing in that local ring.  Thus this rank-3
condition certifies a genuine Hensel-smooth point of the strict transform.

When `rank(J_fib)<3`, the script says `UNRESOLVED`.  Some such points may be
removed by saturation, and saturation could make a remaining component
smooth after an excess component is removed.  None is claimed here.

## Exact balanced-cone table

The following counts are after removing the dummy coordinate and
deduplicating the two blowup charts.

| base support | deep support | raw cones | certified strict-smooth | unresolved |
|:--|:--|--:|--:|--:|
| `b+3` | `L` | 0 | 0 | 0 |
| `b+3` | `v` | 0 | 0 | 0 |
| `b+3` | `U-2v` | 8 | 8 | 0 |
| `b+3` | `U+2v` | 48 | 0 | 48 |
| `a+3` | `L` | 0 | 0 | 0 |
| `a+3` | `v` | 24 | 0 | 24 |
| `a+3` | `U-2v` | 8 | 8 | 0 |
| `a+3` | `U+2v` | 48 | 0 | 48 |
| `a+b+2` | `L` | 0 | 0 | 0 |
| `a+b+2` | `v` | 16 | 0 | 16 |
| `a+b+2` | `U-2v` | 0 | 0 | 0 |
| `a+b+2` | `U+2v` | 40 | 0 | 40 |
| `DB` | `L` | 0 | 0 | 0 |
| `DB` | `v` | 0 | 0 | 0 |
| `DB` | `U-2v` | 8 | 0 | 8 |
| `DB` | `U+2v` | 48 | 8 | 40 |
| `DC` | `L` | 0 | 0 | 0 |
| `DC` | `v` | 8 | 0 | 8 |
| `DC` | `U-2v` | 8 | 0 | 8 |
| `DC` | `U+2v` | 48 | 8 | 40 |
| **total** |  | **312** | **32** | **280** |

Among the 280 unresolved raw cones, 40 have fiber rank 1 and 240 have
fiber rank 2.  There are no rank-0 cones.

Write a cone key as

```text
<t, e_direction, D_direction, X, Y>,
```

with the projective direction normalized to `[r:1]`.  In the `U-2v` and
`U+2v` rows, `(X,Y)=(L,v)`.  The 32 certified cones are four families:

```text
b+3, U-2v:  <4,r,1,L,3>,  r=1,2,3,4,  L=1,4;
a+3, U-2v:  <4,r,1,L,2>,  r=1,2,3,4,  L=2,3;
DB,  U+2v:  <0,r,1,L,4>,  r=1,2,3,4,  L=2,3;
DC,  U+2v:  <0,r,1,L,4>,  r=1,2,3,4,  L=2,3.
```

The first two lie over `(a,b)=(4,2)` and `(2,4)` modulo 5.  The last two
lie over the intersection `(a,b)=(0,0)` of `DB=0` and `DC=0`.

## Pullback to the exact dual square covers

On the chosen transverse parametrizations,

```text
DB=-8e,             DC=-8e.
```

For the dual classes `R3` and `R2`, respectively, impose the exact necessary
square cover

```text
W^2=-8e.
```

On a balanced chart this eliminates `q` as a unit times `W^2`:

```text
chart 1: q=-W^2/(8z),       chart 2: q=-W^2/8.
```

At `W=0`, differentiation again gives precisely `J_fib`; the `W`
derivative vanishes.  Hence rank 3 is the exact first weighted-layer
smoothness test.  It incorporates the even valuation of the normal
parameter and its constrained unit, rather than merely asking whether the
raw core cone exists.

All eight certified `DB,U+2v` cones persist smoothly on `DB=W^2`, and all
eight certified `DC,U+2v` cones persist smoothly on `DC=W^2`.  No other
balanced `DB` or `DC` cone is certified smooth on this square-cover layer.
In particular, the rank-2 `DB/DC,U-2v` cones must not be promoted to dual
branches.

This is a local existence statement on the necessary square cover, not the
final exact dual-halving test.  The direct `R2/R3` divisibility test and exact
torsion computation are still required on any rational specialization.

## Reproduction and next refinement

Run

```text
magma -b code/contact6_m612_boundary_first_blowups_mod5.m
```

The concise result is also recorded in

```text
data/contact6_m612_affine_blowups_mod5_2026_07_10.txt.
```

If the 280 unresolved cones are revisited, use local saturation or
deflation at their maximal ideals.  A second blind first-order solvability
test would repeat the original overclaim.
