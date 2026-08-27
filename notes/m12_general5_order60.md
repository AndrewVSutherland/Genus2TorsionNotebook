# General Mumford order 5 on the compact M(12) surface

## Why this is broader than the point-contact cover

The point-contact calculation in `notes/m12_contact5_exact_cover.md` asks for
an order-5 class of the special form `P-infinity`.  That condition cuts the
two-parameter split `M(12)` surface down to a curve, whose large projected
component is `P55`.

A general rational order-5 class has Mumford support

\[
 q=x^2+Ux+V.
\]

It is characterized instead by the norm identity

\[
 A(x)^2-B(x)^2F(x)=q(x)^5,                 \tag{1}
\]

where `A` is monic of degree 5 and `deg(B)<=2`.  Indeed, reducing (1) modulo
`q` and assuming `gcd(B,q)=1` recovers

\[
 v=-A/B\pmod q,\qquad v^2=F\pmod q.
\]

For the divisor `D=[q,v]`, the function `A+B*y` has divisor `5D`, up to the
choice of sign.  Conversely a generic order-5 Mumford divisor gives (1).
The usual open conditions are

```text
disc(F) * disc(q) * Res(q,F) * Res(B,q) != 0.
```

The coefficient of `x^10` makes the right side of (1) exactly `q^5` after
scaling `A` to be monic.  The coefficients of degrees 9 down to 5 then force
the five remaining coefficients of `A` recursively.  The coefficients of
degrees 4 down to 0 leave five equations in

```text
(b,w,U,V,b2,b1,b0).
```

Thus the general order-5 cover is expected to remain a surface over the
two-dimensional `M(12)` base.  It is materially less restrictive than `P55`.

## Compact M(12) chart

Use

\[
 L=b+(2b-1)x,\qquad H=x+w(1+bx),
\]

and

\[
 F_{b,w}=L\left(LH^2+4b(1+x)^2(wL-x^2)\right).
\]

Every smooth member on

```text
b*w*(b-1)*(2*b-1) != 0
```

carries the marked order-12 class.  The leading coefficient is
`-4*b*(2*b-1)`.  The visible point at `x=-1` has ordinate

\[
 (1-b)(w(1-b)-1),
\]

and its difference from infinity has exact order 12 in every exact survivor
tested below.

## Independent finite masks

The C program

[`code/m12_general5_bw_height_sieve.c`](../code/m12_general5_bw_height_sieve.c)

computes the masks from first principles.  For every residue `(b,w)` modulo
an odd prime `p != 5`, it builds `F`, checks degree and squarefreeness, counts
the curve over `F_p` and `F_(p^2)`, and uses

\[
 \#J(\mathbf F_p)=\frac{\#C(\mathbf F_p)^2+\#C(\mathbf F_{p^2})}{2}-p.
\]

It keeps a good residue only if `5 | #J(F_p)`.  Denominator residues and bad
reductions pass conservatively.  As a positive control, every good residue at
`p != 2,3` is asserted to satisfy `12 | #J(F_p)`.

The first mask counts are

```text
p     good    allowed5    bad
3        0           0      9
7       18           4     31
11      70          16     51
13     110          16     59
17     198          62     91
19     260          72    101
```

These six rows agree exactly with the independent Magma implementation
`code/m12_full_surface_plus5_order60_search.m`; the preserved run is
`m12_general5_bw_mask_crosscheck.log`.

An initial height-3 run produced eight apparent survivors.  Exact inspection
shows that all eight lie on the multiplicity-six discriminant component

\[
 w(b-1)+1=0.
\]

This is a useful positive test that the conservative bad-reduction policy does
not create false exclusions.  The component is now removed exactly before the
rational scan.

## Height-100 and height-200 results

The height-100 run uses every prime through `79` except `2,5`.  It checks
148,172,800 nontrivial rational parameter pairs and leaves no finite-mask
survivor:

```text
SCAN H=100 parameters=12175 pairs=148230625 primes=20
DONE H=100 checked=148172800 survivors=0
```

The log is `data/m12_general5_bw_h100_p79_all.log`.

At height 200 there are 48,927 rational values in each coordinate.  The same
twenty masks check 2,393,618,932 nontrivial pairs and leave ten survivors:

```text
(-61/23,   -87/187)    (-154/43, -52/97)
(183/67,   -51/86)     (141/70,  187/150)
(79/74,    117/145)    (91/114,  -59/142)
(115/134, -160/143)    (124/177,-143/36)
(-32/183,  -76/141)    (-193/198,-26/125).
```

Fresh exact reconstruction with

[`code/m12_general5_bw_verify.m`](../code/m12_general5_bw_verify.m)

gives torsion `[12]` for every one.  Each also has a D4 Frobenius certificate
for geometric simplicity.  Hence none carries rational 5-torsion, and there
is no cyclic-60 hit.  The sieve and exact logs are

```text
data/m12_general5_bw_h200_p79_all.log
m12_general5_bw_verify_h200_p79.log.
```

This proves a bounded-height exclusion in the compact chart, not a global
nonexistence theorem.

## Degree-one B subcover

The first symbolic subcover sets `b2=0` but retains both `b1` and `b0`.  It is
implemented, with exact coefficient recursion and reconstruction assertions,
in

[`code/m12_general5_b2zero_geometry.m`](../code/m12_general5_b2zero_geometry.m).

This slice is expected to be a curve: the base contributes two parameters and
`(U,V,b1,b0)` contribute four more, against five residual equations.  Its
open saturation and projection to `(b,w)` are the next global geometry test.
The constant-`B` locus is lower-dimensional and should be audited separately
if it is removed by leading-coefficient saturation.

## Status

The general-Mumford computation closes an important gap in the earlier search:
the negative result is no longer confined to point-contact 5-torsion.
Nevertheless, no order-60 curve has yet been found.  The two live global tasks
are:

1. factor and normalize the degree-one-`B` subcover;
2. find a rationally parametrized component of the full quadratic-`B` surface,
   then perform exact torsion and D4/root-power certification on its fibers.

