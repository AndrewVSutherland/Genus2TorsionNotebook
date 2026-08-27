# Kuru-Sadek [23] plus 2: quotient and boundary analysis

This continues `notes/order23_extra2.md`.  The question is whether the genus-16
condition curve for a rational finite branch point has maps or low-genus
quotients that could prove absence of rational points with `q=t^2` away from
boundary.

## Projective Boundary

For the original plane model `G(q,r)=0`, the leading homogeneous part is

```text
-3240*q^18*r^5.
```

Thus the projective closure has only the two coordinate points at infinity in
this plane model, namely `[q:r:w]=[1:0:0]` and `[0:1:0]`, with high
multiplicity.  This is consistent with the genus loss from the naive degree-23
plane curve down to geometric genus `16`.

## Square-root/sign model

The root condition is

```text
f_t(r) = A(r)^2 - lambda*r^4*(r-1) = 0.
```

For the Kuru-Sadek family,

```text
lambda*alpha = ((alpha-1)^2/(alpha-beta))^2
```

is a square in `Q(t)`.  Therefore any rational branch point forces
`alpha*(r-1)` to be a square.  Introduce

```text
z^2 = alpha*(r-1),   so   r = 1 + z^2/alpha.
```

The sign equation

```text
A(r)*(alpha-beta)*alpha = (alpha-1)^2*r^2*z
```

clears to a polynomial in `q=t^2` and `z`.  After removing the boundary factor
`(q-1)^3`, the main sign model has

```text
deg_q = 16,
deg_z = 5,
geometric genus = 16,
affine singular points = (0,1), (1,0).
```

So this is a cleaner birational model, not a lower-genus quotient.

The birational map back to the root-condition model is

```text
(q,z) -> (q, r = 1 + z^2/alpha(q)).
```

The inverse is given by the sign expression for `z`, away from the boundary
where `r=0`, `alpha=0`, or `alpha=1`.

## Symmetry checks

The following obvious transformations do not preserve the main sign curve and
therefore do not give immediate quotients:

```text
z -> -z
q -> 1/q
q -> -q
q -> 1-q
q -> q/(q-1)
q -> 1/(1-q)
q -> (q-1)/q
```

For each tested `q`-Möbius transform, the transformed polynomial has gcd `1`
with the original model.  A Magma `AutomorphismGroup` probe returned the genus
`16` but did not finish the automorphism computation within about a minute.
Likewise, built-in hyperelliptic/trigonal/gonality probes did not return
quickly after the genus computation.

## Boundary p-adic lifts

The good affine chart has no square-`q` residues at `p=3` or `p=5`, but this
does not prove absence of rational points because every rational point could
reduce to a bad boundary class.  Lifting the square-pullback equation

```text
H(t,r) = G(t^2,r) = 0
```

confirms that the boundary classes have abundant p-adic lifts.

At `p=3`, the mod-3 classes are

```text
(t,r) = (0,0), (0,1), (1,0), (1,1), (2,0), (2,1).
```

They lift through `3^6`; at level `3^6` there are `39366` affine lifts, with
`6561` above each mod-3 class.

At `p=5`, the mod-5 classes are

```text
(t,r) = (0,0), (0,1), (1,0), (1,1), (2,0), (3,0), (4,0), (4,1).
```

They lift through `5^4`; at level `5^4` there are `50000` affine lifts.  The
classes `(2,0)` and `(3,0)` correspond to the bad `q=-1` chart at `p=5`.

Thus there is no simple local contradiction at `p=3` or `p=5`.  The obstruction
is global: rational points, if any, must live on the boundary-lifting branches
of a genus-16 curve.

## Current conclusion

I did not find a low-genus quotient.  The useful structural result is that the
extra-2 condition has two genus-16 birational models:

```text
G(q,r)=0,  deg=(18,5), singular at q=-1,0,1 boundary points;
S(q,z)=0,  deg=(16,5), singular at (0,1),(1,0).
```

The sign model is better for future work, but it still has genus 16 and no
obvious rational automorphisms.  A proof of no `[46]` in this family likely
requires a serious rational-points computation on this genus-16 curve, not a
small-prime or elementary quotient argument.

Main files:

```text
code/order23_extra2_sign_model.sage
data/order23_extra2_sign_model.txt
code/order23_extra2_padic_boundary.sage
data/order23_extra2_padic_boundary.txt
```
