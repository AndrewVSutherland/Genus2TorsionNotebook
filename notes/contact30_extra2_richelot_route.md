# Richelot lift of the order-30 extra-2 locus

## Target produced by this route

This route naturally targets `[2,60]`, not exact cyclic `[60]`.

In the simultaneous contact-5/contact-6 family,

```text
f = Q2*C3,   deg(Q2)=2, deg(C3)=3.
```

If `C3` has a rational root `rho`, then

```text
f = Q2*Q2b*(x-rho).
```

The three pairs, with `x-rho` paired with infinity, give a pointwise-rational
maximal-isotropic Richelot kernel.  The source generically has torsion
`[2,30]`.  Its odd order-15 subgroup survives the degree-4 isogeny.  If a
dual 2-class is divisible by 2, the Richelot codomain has an order-60 point.

However, the dual kernel is a constant `(Z/2)^2` contained in the rational
2-torsion of the codomain.  Thus an exponent-60 codomain has torsion containing

```text
Z/2 x Z/60,
```

and cannot have exact cyclic torsion `[60]`.

A generic second Richelot step merely quotients by the dual kernel and returns
to the `[2,30]` source.  Reaching cyclic 60 from this lane would require a
special second, non-pointwise-rational Galois-stable kernel (equivalently,
correlations among the three dual quadratic discriminants), or a special
additional principal polarization.  Odd isogenies preserve the rational
2-torsion representation, while a quotient by only one rational order-2
point is naturally `(1,2)`-polarized rather than principally polarized.

## Finite viability

`code/contact30_extra2_richelot_finite.m` constructs the distinguished
Richelot neighbor over finite fields and checks the formula by equality of
source and codomain Jacobian orders.  There were no formula failures in the
test range.  The number of codomains with exponent divisible by 60 was:

```text
p       branch -1       branch +1
11         1/1             2/2
13         3/3             3/3
17         2/6             2/5
19         4/8             4/7
23         3/9             3/9
29        10/23           13/25
31         8/8             8/8
```

The denominator is the number of distinguished rational-root kernels.
Observed finite codomain groups include `[2,60]`, `[2,120]`, `[4,4,60]`, and
`[2,600]`.  Hence the `[2,60]` construction is locally healthy; the bottleneck
is finding rational points on the root cover.

## Compact C3-root geometry

The old equation in `(R,rho)` is irreducible of bidegree `(20,3)`.  The
compact core coordinates `(u,s)` are better.  Let `H(u,s)=0` be the
irreducible degree-10 genus-zero order-30 core and put

```text
D = u^6+6*u^4*s-2*u^4+15*u^3*s-u*s^3+u^2,
N = 15*u^5+90*u^4+20*u^3*s-6*u^2*s^2+231*u^3
    +2*u^2*s-15*u*s^2+90*u^2-20*u*s+15*u-2*s,
q = N/D.
```

The equation `C3(rho)=0` is linear in `q`.  After clearing denominators it is

```text
E = u*rho*(rho-s)*N
    + D*(u*(4*rho^3+(s-6)*rho^2+21*rho)+(u-1)^2) = 0.
```

The symbolic recovery check in `code/contact30_c3root_compact_geometry.m`
gives exactly

```text
num(F2 after q=N/D) = s*(u-1)^2*(u+1)^2*H,
num(F3 after q=N/D) =   (u-1)^2*(u+1)^2*H.
```

Moreover, `E=u*E9` with `E9` irreducible of total degree 9.  Eliminating `s`
factors as

```text
Res_s(H,E)
 = u^8*(u+1)^4*(u-1/8)^2*P8(u)^2*P32(u,rho).
```

Here `P8` is irreducible over `Q` and has no rational affine `u`; its
projective component has genus 0 and only the point at infinity in the tested
point search.  The factors `u=0,-1,1/8` are boundary/base-locus pieces.  The
genuine projected root cover is the irreducible total-degree-32 factor `P32`,
with 265 terms.  Its normalization was deliberately not left running after a
short four-minute attempt.

The compact order-30 search through `u`-height 500 found only the four known
open core points

```text
(125,5415), (125,2715),
(1/125,831/3125), (1/125,-69/3125),
```

plus the `u=1/8, D=0` boundary.  Exact arithmetic shows that none of the four
open points has a rational `C3` root.

## Next step

Study the genuine trigonal curve in its original degree-3-over-`R` model,
rather than normalizing the degree-32 plane projection.  Compute its branch
discriminant and singularity corrections, then its normalized genus and
rational points.  A rational root point would be fed to the distinguished
Richelot construction and the three dual 2-classes would be tested for
halving.  This remains a promising search for the new group `[2,60]`, but it
is secondary to the full `M(12)+5` route for exact cyclic `[60]`.

## Geometry-first resolution (2026-07-11)

The promised go/no-go computation is now complete, and this route should be
parked.

- The normalized `C3`-root cover has genus `12`.
- Its sole rational involution gives a genus-`6` trigonal quotient with
  trivial rational automorphism group.
- The necessary descended dual-halving covers have genera `17` and `18`.
  The simultaneous third gate is a genus-`47` `V4` cover.
- A projective 61-prime sieve exhaustively tested `759,920,359` rational
  parameters of height at most `25000` on both signs. Its only survivors
  were `R=1,2,3,7/3`, all exact `c=0` degeneracies; there were no open root
  points.

The exact geometry, bounded search, and reusable Richelot verifier are in

```text
notes/contact30_c3root_trigonal_geometry_2026_07_11.md
notes/contact30_c3root_projective_search_2026_07_11.md
notes/contact30_c3root_richelot_verification_2026_07_11.md
```

The finite-field viability remains correct, but reaching it over `Q` now
requires high-genus rational-point machinery of the same kind this pivot was
intended to avoid.
