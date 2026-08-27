# Nonautomorphic `[2,6,6]` surface diagnostic

Date: 2026-07-10

This follows `contact6_m36_266_targeted_lift.md`.  It replaces deeper lifting
on the fixed `(r,b)=(4,5)` fiber by saturation and component diagnostics on
the extra-root/cubic-contact equations.

The reproducible files are:

```text
code/contact6_m36_266_surface_decompose.sage
data/contact6_m36_266_surface_decompose_p13_r4_full_open.txt
data/contact6_m36_266_surface_decompose_p19_r4_full_open.txt
data/contact6_m36_266_surface_points_p13.txt
data/contact6_m36_266_surface_points_p19.txt
```

## Scheme and open conditions

On the `eps=+1` branch, put

```text
a = 3 - (b+3)*r - 2/r,     M=L^2.
```

The three contact equations in `(r,b,M,U,v)` have total degrees
`8,10,7`.  The structural saturation used

```text
r*(r-1)*(b+3)*M*v*(U^2-4*v^2)*numerator(a+b+2).
```

The full open diagnostic then also saturated by `disc(f)` and
`Res_x(q,f)`.  Finally, it saturated sequentially by the two ideals generated
by the explicit Mobius-pairing equations `AutoA` and `AutoB`.  This last
operation removes components contained in an automorphism locus; it does not
erase the lower-dimensional automorphism divisor from the closure of a
nonautomorphic component.

## Fixed `r=4` slice

This is a curve slice of the full surface, not a decomposition of the full
two-dimensional surface.  It contains the genuinely nonautomorphic local
point used previously:

```text
p=13: (r,b,M,U,v) = (4,5,1,4,6),
p=19: (r,b,M,U,v) = (4,5,6,9,3).
```

Over both `F_13` and `F_19` the computations give exactly the same Hilbert
data:

```text
raw fixed-r ideal: dimension 2, affine degree 14
after structural saturation: dimension 1
Hilbert polynomial: 95*t - 355
total affine degree: 95
```

The raw dimension-2 part is boundary.  Saturation by `disc(f)` and
`Res_x(q,f)` has exponent zero and does not change the degree-95 curve.
Saturation by `AutoA` and `AutoB` also has exponent zero and does not change
its ideal.  Thus the closure of the smooth nonautomorphic open is a
one-dimensional scheme of total degree 95 in this slice in both tested
characteristics.  In particular, the nonautomorphic point is not an isolated
fixed-fiber accident and the slice is not supported on the two known
automorphism loci.

The degree 95 is the degree of the whole saturated scheme, not the degree of
an individually certified irreducible component.  Primary decomposition at
`p=19` did not finish within 240 seconds.  A separate primality test at
`p=13` was subjected to the same bound.  No component count or individual
component degrees are claimed.

The successful commands were:

```text
sage code/contact6_m36_266_surface_decompose.sage \
  --prime 13 --fixed-r 4 --full-open \
  --output data/contact6_m36_266_surface_decompose_p13_r4_full_open.txt

sage code/contact6_m36_266_surface_decompose.sage \
  --prime 19 --fixed-r 4 --full-open \
  --output data/contact6_m36_266_surface_decompose_p19_r4_full_open.txt
```

## Full-surface finite evidence

Enumeration was done in the `M`-quotient, without assuming that `M` is a
square.  Every nonautomorphic good point in the following table has contact
Jacobian rank 3, hence tangent dimension 2 as expected for the full surface.

```text
       good M-points  nonauto  smooth nonauto  square-M nonauto  bases (r,b)
p=13       42            34          34                28             20
p=19      192           152         152                74             76
```

Each nonzero square `M` has two choices of `L`.  Hence the square counts give
`56` and `148` nonautomorphic `L`-points, exactly matching the earlier direct
enumeration.  This independently checks the `M=L^2` quotient and the
automorphism filtering.

These counts were produced by:

```text
sage code/contact6_m36_266_surface_decompose.sage \
  --prime 13 --count --count-only \
  --output data/contact6_m36_266_surface_points_p13.txt

sage code/contact6_m36_266_surface_decompose.sage \
  --prime 19 --count --count-only \
  --output data/contact6_m36_266_surface_points_p19.txt
```

## Bounded stops and conclusion

A full five-variable saturation over `F_13` reached the raw Groebner basis
but not the structural saturation within 240 seconds.  The exact `Q`, `r=4`
structural saturation also exceeded 240 seconds.  No process was left
running, and no exact rational point was found.

The calculation rules out the proposed strategy of searching for a separate
"nonautomorphic component" obtained merely by saturating away `AutoA` and
`AutoB`: on the tractable `r=4` slice those loci are proper closed
subschemes, not whole components, so the saturation leaves the curve
unchanged.  The arithmetic problem is rational points on a high-degree curve,
not selection of an already separated low-degree component.

The next justified move is to compute the minimal associated primes and
normalizations of the degree-95 `r=4` curve with a system capable of finishing
that decomposition, then retain only the component through
`(b,M,U,v)=(5,1,4,6)` modulo `13` and `(5,6,9,3)` modulo `19`.  Only if that
component has a low-genus quotient should Hensel/CRT reconstruction resume.
Without that component match, more fixed-`b` lifting or a larger rational
height box has no geometric justification.
