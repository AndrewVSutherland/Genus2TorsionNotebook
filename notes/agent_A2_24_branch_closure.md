# A(2,12) -> A(2,24) branch closure

Date: 2026-07-02.

This closes the non-quartic pieces left after
`agent_A2_24_quartic_extract.m` for the four best fibers

```text
(-1/3,-1, 4/3),  (-1/3,1, 4/3),
( 1/3,-1,-4/3),  ( 1/3,1,-4/3).
```

Code and logs:

```text
code/agent_A2_24_branch_closure.m
results/A2_24_branch_closure_compact.log
results/A2_24_branch_closure_default.log
```

Commands run:

```text
magma -b PrintPolynomials:=false \
  code/agent_A2_24_branch_closure.m \
  > results/A2_24_branch_closure_compact.log 2>&1

magma \
  code/agent_A2_24_branch_closure.m \
  > results/A2_24_branch_closure_default.log 2>&1
```

Both runs ended with

```text
A2_24_BRANCH_CLOSURE_DONE fibers=4 translated_classes=16
affine_classes=16 affine_factors=32 affine_factors_with_Q_root=0
rational_point_classes=0 boundary_rational_classes=0 exact_divisible_classes=0
```

The common curve is the one recorded in `agent_A2_24_quartic_extract.md`:

```text
R = 1/5*x^2 + 1/5*x + 7/5
F = 9/25*x^4 + 33/25*x^3 + 3*x^2 + 56/25*x + 16/25
f = R*F.
```

The four translated exact order-12 classes which occur over the four sign
fibers are

```text
D_A     = [x^2 - 9*x + 2,              -29*x + 8]
D_B     = [x^2 + 6/11*x + 52/11,       67/121*x - 160/121]
D_minus = [x^2 - 2/3*x - 4/3,          -7/3*x - 8/3]
D_plus  = [x^2 + 13/3*x + 26/3,        5/3*x + 16/3].
```

`D_A` and `D_B` are the O/TR classes; depending on the sign fiber they swap
between `P12` and `P12+TR`.  The two extra 2-torsion translations give
`D_minus` and `D_plus`, as in the quartic extraction note.

## Affine saturated factors

For every one of the 16 translated classes,

```text
raw resultant degree in M = 32
s4 = M^2 - 9/125
boundary contribution = [ <2, 8> ]
saturated affine degree = 16
gcd(E1,E0) = 1
```

After removing `s4=0`, the saturated affine factorization is:

| class | saturated M-factor degrees | Q-factorization degrees | rational M-roots |
|---|---:|---:|---:|
| `D_A` | `[16]` | `[16]` | none |
| `D_B` | `[16]` | `[16]` | none |
| `D_minus` | `[4,4,8]` | `[4],[4],[8]` | none |
| `D_plus` | `[4,4,8]` | `[4],[4],[8]` | none |

Thus there is a projection obstruction on the saturated affine chart: any
rational halving point would project to a rational root of one of these
univariate M-factors, but all 32 factors checked in the four fibers have
`rational_M_roots=[]`.

The quartic factors are exactly the three quartics `qA,qB,qC` from
`agent_A2_24_quartic_extract.md`.  The two residual degree-8 factors are:

```text
r_minus = M^8 + 60/13*M^7 + 51066/4225*M^6 + 51408/4225*M^5
          + 1070739/105625*M^4 + 2395008/528125*M^3
          + 3873906/2640625*M^2 + 3577932/13203125*M
          + 1266273/66015625

r_plus  = M^8 + 192/143*M^7 + 110412/511225*M^6
          - 8892288/28117375*M^5 - 174495546/1546455625*M^4
          + 225364032/7732278125*M^3
          + 741675852/38661390625*M^2
          - 1160054784/193306953125*M
          + 440118441/966534765625.
```

Magma factors both as irreducible degree-8 polynomials over `Q`, with no
rational roots.  The two O/TR degree-16 polynomials are printed in the default
log; each is irreducible over `Q` and has no rational roots.

## Boundary check

The removed boundary is always

```text
s4 = M^2 - 9/125.
```

Since `Roots(s4)` over `Q` is empty, there are no rational points on the
`s4=0` boundary.  As a sanity check, over `Q(sqrt(5))` the specialized system
has `gcd_N_degree=3`, but the specialized `S` has degree `3`, not `4`; these
are the expected discarded degenerate square-quartic boundary solutions, not
legitimate rational halvings.

## Exact check

For each translated class the script also maps to the integral square model

```text
y_I^2 = 125^2*f,  [u,v] -> [u,125*v],
```

and calls Magma's exact `IsDivisibleBy(D,2)`.  All 16 checks return `false`.

## Verdict

No rational halving occurs in any translated order-12 class for any of the
four best fibers.  The affine saturated cover has no rational M-projection,
the `s4=0` boundary has no rational M-coordinate, and Magma's exact
integral-model divisibility check agrees.  Consequently no candidate reaches
the A(2,24) torsion-certification stage.
