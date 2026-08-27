# Rational `a=1/e` scan on the exact orthogonal `[3,3]` cover (2026-07-11)

## Result

The exact degree-12 support polynomial `f12_a(v)` and exact recovery
`M_a(v)` were scanned over every reduced rational

```text
a=n/d,  |n| <= 200,  1 <= d <= 200.
```

There are `48,927` such parameters.  Finite-field masks at

```text
p = 7,11,13,17,19,23,29,31,37,41,43
```

retained a residue only if it had an `F_p` root of `f12_a` on which `M_a`
was a square.  Projective-infinity and coefficient-pole residues were marked
bad and retained.  Thus a rejected rational parameter cannot carry a
rational signed point on the exact affine cover.

Only four parameters survived all masks:

```text
a = -3, 0, -15/8, 43/48.
```

Exact checks give **no open signed hit**:

| `a` | exact outcome |
|---|---|
| `-3` | `e=-1/3`; recovery-pole fiber and singular source.  Indeed the source has an `x^2` factor. |
| `0` | `e=infinity`; resolved separately below.  The direct fiber has no rational lift on the limiting degree-12 support. |
| `-15/8` | `e=-8/15`, `v=-3/4`, `M=L=0`, `U=-3/2`; the source has `(x-3/4)^2`, so this is boundary and fails the contact-open condition. |
| `43/48` | `f12_a` has no rational root. |

The identical height-100 scan checks `12,175` parameters and has the same
four survivors and no open hit.

This is a rigorous bounded result, not a global determination of rational
points on the signed cover.

## The `a=0` chart

The rational recovery formula was reconstructed in the coordinate `e=1/a`,
so its coefficients are not regular at `a=0`.  The source curve itself is
smooth there, so the fiber was recomputed directly from the three contact
equations rather than discarded.

After saturation by `M`, the direct `a=0` algebra has length `40` and its
reduced `v`-resolvent factors in degrees

```text
1, 1, 8, 24.
```

The degree-12 polynomial specializes as

```text
(v+1)^4 *
(v^8 - 4v^7 + 16v^6 - 34v^5 + 40v^4
     - 34v^3 + 16v^2 - 4v + 1).
```

The degree-8 factor has no rational root.  Resolving the only rational
support value `v=-1` in the full contact ideal gives

```text
M + U^3/4 + 3U^2/2 + 3U + 2 = 0,
U^4 + 8U^3 + 20U^2 + 20U + 20/3 = 0.
```

The quartic has no rational root.  Hence `a=0` contributes no rational
support lift, signed or otherwise.

## Reproduction

The mask and exact-survivor driver is

```text
code/contact6_m612_relative3_rational_a_scan_continuation.m
```

The standalone driver composes it with the exact reconstruction, so the
complete height-200 calculation is one command with no source edits:

```bash
bash code/contact6_m612_relative3_rational_a_scan.sh 200 43
```

The direct infinity-chart check is:

```bash
timeout 180s magma -b code/contact6_m612_relative3_a0_fiber.m
```

The height-200 run took about 27 seconds on the development machine,
including exact reconstruction.
