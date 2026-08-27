# Saturated global halving cover over an A(2,12) slice

This continues the `A(12) -> A(2,12) -> A(2,24)` lane after the
local-survivor check in `agent_A2_24_halving_cover.md`.  The main rule here is
to work with the square-quartic halving equations and remove the chart boundary
`s4=0`; the raw degree-32 resultants are only useful as a way to identify the
boundary contribution.

## Chosen A(2,12) branch/slice

I used the concrete split slice already recorded in `a2_12_resolvent.tex` and
`a2_12_resolvent.m`:

```text
p = -5/3,  z = 1,
```

and then took the fully split quadratic-factor cover of the residual quartic
`F`.  If `F_m = x^4 + B*x^3 + C*x^2 + D*x + E`, the A(2,12) branch is
described by

```text
Phi_r(eta) = eta^3 - C*eta^2 + (B*D - 4*E)*eta
             + (4*C*E - B^2*E - D^2) = 0,
delta^2 = eta^2 - 4*E,     delta != 0.
```

The recovered quadratic factor is

```text
u0 = (eta - delta)/2,
v0 = (eta + delta)/2,
a  = (D - B*u0)/delta,
U  = x^2 + a*x + u0.
```

For the translated order-12 class I used

```text
D12 = P12 + T_extra,       T_extra = [U, 0].
```

At the known rational split point `(p,z,r)=(-5/3,1,2/3)`, this is

```text
U = x^2 - 136/25*x + 16/3,
D12 = [ x^2 - 308/225*x - 4/9,
        1339/3375*x + 32/135 ].
```

The script checks `D12` has exact order 12 and that `T_extra` is an independent
extra 2-torsion class.

## Saturated halving equations

For `D12=[u,v]`, introduce `M,N` and set

```text
ell_D = v + u*(M*x + N),
S     = (ell_D^2 - f)/u
      = s4*x^4 + s3*x^3 + s2*x^2 + s1*x + s0.
```

The square-quartic equations are

```text
E1 = 8*s4^2*s1 - s3*(4*s4*s2 - s3^2),
E0 = 64*s4^3*s0 - (4*s4*s2 - s3^2)^2.
```

The saturated affine chart is

```text
E1 = E0 = 0,      s4 != 0.
```

This is implemented in
`code/agent_A2_24_saturated_global_cover.m`.

## Exact fiber invariant

On the rational split fiber `(p,z,r)=(-5/3,1,2/3)`, the script prints:

```text
E1: total degree 6, degree_M 6, degree_N 3
E0: total degree 8, degree_M 8, degree_N 4
s4 = M^2 + 225/1331
resultant_N_degree = 32
boundary_factor_degrees = [ <2, 8> ]
saturated_affine_degree_after_s4 = 16
affine_factor_degrees = [ <16, 1> ]
boundary_s4_rational_points = 0
```

So the apparent degree 32 again contains exactly the `s4=0` boundary with
multiplicity 8.  After saturation the honest affine halving fiber has degree
16, as expected for halving a genus-2 Jacobian class.  Unlike the previous
fixed local survivor, this split fiber gives one irreducible degree-16 orbit,
not rational or degree-8 components.

## Finite-field counts on the slice

The same script enumerates the oriented fully split A(2,12) cover over small
good finite fields and counts solutions of `E1=E0=0` with `s4 != 0`.

```text
q   split delta!=0   smooth split   order-12 translates   saturated points   s4-boundary removed
7   2                2              2                     0                  0
11  4                2              2                     0                  0
13  10               8              8                     0                  8
17  14               4              4                     0                  0
19  12               6              6                     16                 4
23  18               8              8                     14                 0
29  16               10             10                    0                  6
37  30               14             14                    16                 10
41  34               24             24                    0                  16
```

These are point counts on the selected slice, not a proof about every
component of the full A(2,12) cover.  Still, they are consistent with the exact
fiber calculation: after removing `s4=0`, this branch does not expose a
rational or otherwise low-degree component.

## Decision

This branch is not a good Chabauty/descent target as it stands.  The saturated
fiber has degree 16 and is irreducible on the rational split point, while the
finite-field counts on the slice are sparse after removing the boundary.  I
would only return to this slice if a later computation finds a quotient or a
special subbranch of the A(2,12) resolvent with smaller factor degrees.
