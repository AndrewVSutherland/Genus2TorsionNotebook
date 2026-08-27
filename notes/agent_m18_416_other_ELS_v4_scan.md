# Other ELS Fibers: Corrected S_B V_4 Scan

Date: 2026-07-02

This applies the corrected `S_B` genus-3 `V_4` method to the other ELS
fibers from `agent_m18_416_fiber_local_scan_h30.log`.

Reproducible scripts:

```sh
magma code/agent_m18_416_SB_v4_scan_els.m
magma -b Case:="-8" SearchHeight:=1000 code/agent_m18_416_SB_v4_scan_els.m
magma code/agent_m18_416_check_stage_AB_points.m
```

## R = -29/8

This fiber is killed by the same style of rank-zero quotient certificate.

The quotient `E_Y` has minimal model

```text
y^2 = x^3 + x^2 - 1269543240180*x - 20883341354391072
rank bounds = 0
torsion = [2, 2]
```

The forced lambda values give only:

```text
lambda = -7424/21: X=m^2=12789/128, not a rational square
lambda = -5376/37: X=m^2=39701/512, not a rational square
lambda = -172288/777: X=0 boundary
lambda = -256: denominator boundary
```

So `R=-29/8` has no nondegenerate rational points on the corrected
`S_B` cover.  A direct lambda search to height `200` found no
nondegenerate `S_B` hits.

## R = -8

This fiber does **not** get a rank-zero killer quotient from the corrected
`S_B` split.  The three elliptic quotients all have rank `1` and torsion
`[2,2]`:

```text
E_m:   y^2 = x^3 - 1263*x + 17138
E_Y:   y^2 = x^3 - 2847*x + 29986
E_mY:  y^2 + x*y + y = x^3 + x^2 - 326*x + 1874
```

The lambda search to height `1000` found exactly two nondegenerate
`S_B` hits:

```text
lambda = -5/7:
  m = +/-28, S = 448, Y = +/-84, w = -8

lambda = -11/9:
  m = +/-36, S = -576, Y = +/-108, w = 8
```

Exact second-stage diagnostics show that the `B` component passes for
both, but the `A` component degenerates:

```text
A polynomial = x^2 + 256*x + 16384
A discriminant = 0
B polynomial = x^2 + 49*x + 256
B discriminant = 1377
B quadratic-field square = true
```

Thus the visible `R=-8` `S_B` points do not produce valid full
second-stage points.  Unlike `R=-25/4` and `R=-29/8`, however, `R=-8`
is not certified empty by this rank-zero quotient method; it needs a
rank-1 intersection/Mordell-Weil sieve or a direct analysis of the
genus-3 `V_4` curve.

## Status

Among the three ELS fibers from the first local scan:

```text
R = -25/4  killed by rank-zero E_mY quotient; only d_B=0 branch points
R = -29/8  killed by rank-zero E_Y quotient; forced X=m^2 nonsquare/boundary
R = -8     not killed by rank-zero quotient; only small S_B hits found fail A by discriminant 0
```
