# Exact obstruction on the weighted R3 `P16` component (2026-07-10)

## Result

The irreducible projection factor `P16(e,mu)` carrying the certified
weighted-endpoint `E9` branches over `Q_5` has no rational point on its
normalization.  More precisely, its normalization has no point over either
`Q_2` or `Q_3`.

This closes the `P16` component only.  It does **not** rule out the other
unresolved affine or weighted-boundary strata in the `[6,12]` search.

## Exact genus-one model

The plane factor has bidegree `(16,8)`, is irreducible over `Q`, and has
normalized genus one.  Its four rational plane singularities resolve into
places of degrees

```text
(-3/5:0:1)   [8]
(-6/17:0:1)  [2]
(0:0:1)      [2,2]
(0:1:0)      [8].
```

Thus none of the visible rational singular plane points is a rational point
of the normalization.  Let `D` be the degree-two place above
`(-6/17:0:1)`.  Exact Riemann--Roch computations give

```text
deg(D)=2,  l(D)=2,  l(2D)=4,  l(4D)=8.
```

The products of a basis of `L(2D)` have exactly two quadratic relations,
giving a degree-four genus-one normal model.  Independently, choose a basis
`x,z` of `L(D)` and a complement `y` to
`<x^2,xz,z^2>` in `L(2D)`.  The nine products

```text
y^2, y*x^2, y*x*z, y*z^2, x^4, x^3*z, x^2*z^2, x*z^3, z^4
```

have a one-dimensional relation space in `L(4D)`.  This produces a degree-two
genus-one model in the same function field.  Magma's exact `Minimise` and
`Reduce` transformations send it to

```text
Y^2 = 2 X^4 - 6 Z^4.
```

The degree-four and degree-two constructions both have Jacobian

```text
E: y^2 = x^3 + 3*x.
```

The diagnostic prints the exact minimization and reduction transformations;
in this run the final change-of-variables matrix was
`[[0,1],[-1,0]]`.  The Riemann--Roch construction is birational: `X/Z` has
degree two on the original normalization (its pole divisor is `D`), and the
binary-quartic equation is quadratic over `Q(X/Z)`.

## Elementary local obstructions

For a local point on the weighted projective quartic, scale so that `X,Z`
are integral and `min(v_p(X),v_p(Z))=0`; then `Y` is integral as well.

At `p=2`, at least one of `X,Z` is odd.  Modulo 16, a fourth power is `1`
for an odd input and `0` for an even input.  The right side
`2 X^4-6 Z^4` is therefore respectively `2`, `10`, or `12` modulo 16,
according as only `X`, only `Z`, or both are odd.  None is a square modulo
16.  Hence there is no `Q_2` point.

At `p=3`, if `X` is a unit then reduction modulo 3 gives
`Y^2 = 2 X^4 = 2`, impossible.  Otherwise `Z` is a unit and `3 | X`.
Reduction modulo 3 forces `3 | Y`, while modulo 9 the right side is
`-6 Z^4 = 3`; this contradicts `Y^2 = 0 (mod 9)`.  Hence there is no
`Q_3` point.

For comparison, the model is soluble over `Q_5`: modulo 5,
`(X,Z,Y)=(0,1,2)` is smooth and lifts by Hensel.  This is consistent with the
previously certified `E9` `Q_5` branches, but the `2`- and `3`-adic
obstructions prevent them from globalizing to a rational point.

## Reproduction

From `torsion_jac` run

```bash
magma -b code/contact6_m612_weighted_R3_p16_independent.m
```

The run takes about two seconds on the current machine and asserts the two
Riemann--Roch relation dimensions, the reduced model, its Jacobian, Magma's
local-solubility results, and zero primitive congruence solutions modulo 16
and modulo 9.
