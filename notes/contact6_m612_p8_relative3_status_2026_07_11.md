# `P8` relative-3 status: exact cover, bounded sieve, and local disks

Date: 2026-07-11.

## Current conclusion

The characteristic-zero cover is explicit.  The relevant orthogonal support
component has degree `12`, and recovering the sign gives a connected
degree-`24` cover.  No rational point on that cover is known, and no global
no-point theorem has been proved.

The main results now fit together as follows.

| calculation | certified result | limitation |
|---|---|---|
| exact support algebra | irreducible `Phi12(v)` over `Q(e)`, with exact maps `M(v),U(v)` | describes supports `{Q,-Q}` |
| exact sign lift | `M` is nonsquare; `Phi24(L)=Res_v(Phi12,L^2-M)` is irreducible of degree `24`; `L` is primitive | no rational-point classification yet |
| P8 base change | exact coefficientwise pullback by `e=e(u)`; connected irreducible signed cover of degree `24` | no rational-point classification yet |
| genus | exact pre-P8 support/signed genera `5/10`; exact P8 support/signed genera `73/145` | the values `73/145` are independently reproduced at `7` and `13` |
| relative height-`200` census | all `48,927` reduced `a=1/e` values tested; four mask survivors, no open signed hit after exact checks | bounded result only |
| height-`1000` sieve | `1,216,767` reduced parameters tested; only `u=-3,2`, both `e=0`, survive | bounded result only |
| `p=7` disks | zero-disk E3 layers excluded; two pole disks have signed orthogonal branches | endpoint signatures beyond the bounded fan remain unresolved |
| `p=17` disks | endpoint E3 and common-root affine charts partly excluded; four `DC=0` disks have signed orthogonal branches | other endpoint weights and nonintegral node blowups remain unresolved |

Thus neither `p=7` nor `p=17` obstructs the route: each has genuine local
points on the target signed cover.  The height sieve is strong evidence but
does not replace a global argument.


## Exact cover

The marked rational `3`-line splits the `40` nonzero supports modulo sign as

```text
40 = 1 + 12 + 27.
```

The degree-`1` term is the repeated marked support, the degree-`12` orbit is
the orthogonal complement relevant over `Q`, and the degree-`27` orbit is
nonorthogonal.  The exact reconstruction proves that the degree-`12` factor
is irreducible and contact-open.  It also provides degree-`11` recovery maps
for all contact coordinates.

In the degree-`12` field, `M` is nonsquare.  Hence the signed class requires

```text
L^2=M,
```

and has degree `24`.  The exact resultant `Phi24(L)` is irreducible, while a
linear gcd recovers `v` from `L`; this proves that `L` itself is primitive.
The formulas and verification are in
`notes/contact6_m612_p8_relative3_exact_2026_07_10.md`.


## Exact genera and modular checks

The exact fields before P8 base change satisfy

```text
support genus = 5,
signed genus  = 10.
```

The divisor of `M` has exactly two odd places, both of degree `1` above
`e=0`.  Thus the signed quadratic cover has ramification degree `2`, and
Riemann--Hurwitz gives

```text
2*5-1+2/2 = 10.
```

For the P8 pullback, the relative different has degree `32`.  Exact pullback
of its places through `e=e(u)` gives total different degree `168`, so
Riemann--Hurwitz gives

```text
support genus = 73,
signed genus  = 145.
```

The exact P8 map has ramification index `2` over `e=0`, so the odd valuations
of `M` become even after normalized base change.  The P8 sign lift is
therefore unramified; exact irreducibility of its degree-`24` polynomial
certifies connectedness, giving `2*73-1=145`.  Independent computations over
both `F_7(u)` and `F_13(u)` reproduce `73/145`.


## Bounded rational-parameter sieve

The residue sieve used every prime

```text
7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71
```

on every reduced parameter `u=n/d` with `|n|,d<=1000`.  Of `1,216,767`
parameters, only

```text
u=-3,2
```

survived all masks, and both are the boundary `e=0`.  Therefore there is no
nonboundary specialization in this height box.  The projective parameter at
infinity duplicates the finite `e`-fiber at `u=12` and is also covered.


## Bounded census on the genus-10 relative curve

An independent scan worked directly with the exact support polynomial and
the square condition `L^2=M` over every reduced

```text
a=n/d,  |n|<=200,  1<=d<=200.
```

Masks at the primes `7` through `43` reduced the `48,927` parameters to

```text
a=-3, 0, -15/8, 43/48.
```

Exact specialization leaves no contact-open signed point.  The values
`a=-3` and `a=-15/8` are singular/boundary fibers, with the latter forcing
`M=L=0`; `a=43/48` has no rational support root.  The `a=0` infinity chart
was recomputed directly: its only rational support candidate `v=-1` forces
an irreducible quartic in `U`, so it has no rational lift.  This proves the
bounded height-`200` statement, not global emptiness of the genus-`10`
curve.


## Forced disks at `p=7`

The retained residues are `u=1,2,4,6`.  In rational representatives they
split into

```text
u=2,-3: double-zero endpoint disks,
u=1,-1: parameter-pole disks.
```

The uniformly scaled E3 system is empty at every depth in the zero disks,
and the previously enumerated endpoint signatures through `v(e)=8` add no
point.  This does not exhaust all higher-weight signatures.

Every punctured pole disk has two smooth signed branches with exact Weil
pairing `1`.  They are therefore points of the target degree-`24` cover, so
`p=7` cannot close the global route.


## Forced disks at `p=17`

The ten retained disks have three types.

- At the two `e=0` disks, every uniformly scaled E3 layer is empty, but other
  endpoint weights are not exhausted.
- At the four `a=-2` common-root disks, the raw affine signed contact chart is
  empty, but nonintegral weighted blowups have not been classified.
- At the four `DC=0` disks, each disk has two smooth signed orthogonal
  branches on the degree-`24` cover.

Thus `p=17` also supplies local points rather than a global obstruction.


## Best next work

1. Attack rational points on the exact genus-`10` signed curve over the
   `e`-line first.  Compute its Jacobian and exploit the natural `S3`
   quotients before moving to the much larger P8 pullback.
2. Develop global rational-point methods for the exact genus-`145` P8 signed
   cover after using every lower-genus quotient available over the `e`-line.
3. Analyze the forced disks at `p=19,23,41`.  A single prime can close the
   route only if every retained disk is obstructed.
4. Retain the unresolved high-weight endpoint and nonintegral collision
   charts explicitly; do not promote bounded fan calculations to complete
   local emptiness.
5. Do not prioritize a larger undirected height scan: the height-`1000`
   sieve already has no nonboundary survivor.


## Detailed notes and reproduction

```text
notes/contact6_m612_p8_relative3_exact_2026_07_10.md
notes/contact6_m612_p8_relative3_sieve_2026_07_10.md
notes/contact6_m612_relative3_rational_a_scan_2026_07_11.md
notes/contact6_m612_p8_p7_bad_disks_2026_07_11.md
notes/contact6_m612_p8_p17_bad_disks_2026_07_11.md
```

```bash
magma -b code/contact6_m612_relative3_exact_reconstruct.m
magma -b relative_only:=true do_signed:=false \
  code/contact6_m612_p8_relative3_exact_genus.m
magma -b p8_rh_only:=true \
  code/contact6_m612_p8_relative3_exact_genus.m
magma -b p:=7  do_genus:=true code/contact6_m612_p8_relative3_modular.m
magma -b p:=13 do_genus:=true code/contact6_m612_p8_relative3_modular.m
bash code/contact6_m612_relative3_rational_a_scan.sh 200 43
magma -b code/contact6_m612_relative3_a0_fiber.m
magma -b mode:=search height:=1000 prime_bound:=71 max_exact:=500 \
  code/contact6_m612_p8_extra3_residue_sieve.m
magma -b code/contact6_m612_p8_p7_bad_disks.m
magma -b code/contact6_m612_p8_p17_bad_disks.m
```

The exact relative signed genus is also exposed by
`code/contact6_m612_relative3_exact_m_divisor_continuation.m` and
the P8 base ramification by
`code/contact6_m612_relative3_exact_p8_base_ramification_continuation.m`.
