# R = -25/4 Corrected S_B Certificate

Date: 2026-07-02

This supersedes the earlier genus-5 interpretation.  The reduced quartic
fiber coordinate differs from the unreduced coordinate by the rational
factor `609/256`.  Including that factor changes the `S_B` pullback from
the spurious genus-5 model to a genus-3 `V_4` cover.

## Corrected Model

Use the reduced fiber coordinate

```text
S^2 = 1024*m^4 - 865600*m^2 + 231800625.
```

The corrected `S_B` condition is

```text
Y^2 =
  29/50*m^4 - 16907/32*m^2 + 268888725/2048
  + (29/1600*m^2 - 17661/2048)*S.
```

Equivalently,

```text
51200*Y^2 - 29696*m^4 - 928*m^2*S
+ 27051200*m^2 + 441525*S - 6722218125 = 0.
```

Magma reports the projective closure has normalized genus `3`.

## V_4 Quotient

Set

```text
lambda = (S - 15225)/(32*m^2).
```

Away from boundary points,

```text
m^2 = 25/16 * (609*lambda + 541)/(1 - lambda^2),

Y^2 =
725/1024 *
(29*lambda + 21)*(609*lambda + 541)*(21*lambda + 25)
/((lambda - 1)^2*(lambda + 1)).
```

The quotient controlling simultaneous square-lifting is

```text
C^2 = -29*(lambda - 1)*(29*lambda + 21)*(21*lambda + 25).
```

Its minimal model is

```text
y^2 = x^3 - x^2 - 415273465*x + 2669108887225,
```

with

```text
rank bounds: 0
torsion: [2, 2]
```

Thus the only rational `lambda` values on this quotient are

```text
lambda = infinity, 1, -21/29, -25/21.
```

There is also the cleared-model boundary `lambda = -541/609`,
corresponding to `m=0`, `S=15225`, `Y=0`, which is outside the original
affine `S_B` condition because the formulas had denominators in `m`.

## Remaining Affine Lifts

The finite affine lifts are exactly

```text
lambda = -21/29:
  m = +/-145/8, S = 15225/2, Y = 0, w = -71/13

lambda = -25/21:
  m = +/-105/4, S = -11025, Y = 0, w = 71/13
```

For all four values,

```text
d_B numerator = 0.
```

So the exact `S_B` cover has rational branch/boundary points, but no
nondegenerate rational points on the `R=-25/4` fiber.  These branch
points cannot pass the original second-stage test, where a zero
discriminant is rejected.

Reproducible script:

```sh
magma code/agent_m18_416_R25_4_SB_v4_certificate.m
```
