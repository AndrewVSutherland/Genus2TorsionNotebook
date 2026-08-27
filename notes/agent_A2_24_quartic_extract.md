# A(2,12) -> A(2,24) quartic branch extraction

Date: 2026-07-02.

This note extracts the explicit quartic saturated halving branches for the four
best third-pass fibers

```text
(-1/3,-1, 4/3),  (-1/3,1, 4/3),
( 1/3,-1,-4/3),  ( 1/3,1,-4/3).
```

Code and log:

```text
code/agent_A2_24_quartic_extract.m
results/A2_24_quartic_extract_default.log
```

Run:

```text
magma code/agent_A2_24_quartic_extract.m
```

The default run used `SearchHeight=8` and ended with

```text
A2_24_QUARTIC_EXTRACT_DONE fibers=4 extra_classes=8
quartic_branches=16 quartics_with_Q_root=0
direct_Q_points_height_8=0
```

## Common fiber data

All four sign fibers give the same genus-2 curve data:

```text
R = 1/5*x^2 + 1/5*x + 7/5
F = 9/25*x^4 + 33/25*x^3 + 3*x^2 + 56/25*x + 16/25
  = (9/25)*(x^2 + x + 1/3)*(x^2 + 8/3*x + 16/3)

f = R*F
  = 9/125*x^6 + 42/125*x^5 + 171/125*x^4
    + 362/125*x^3 + 597/125*x^2 + 408/125*x + 112/125.
```

The cyclic rational 2-torsion class is

```text
TR = [x^2 + x + 7, 0].
```

The two extra rational 2-torsion translations are

```text
T1 = [x^2 + x + 1/3, 0],
T2 = [x^2 + 8/3*x + 16/3, 0].
```

For both `T1` and `T2`, in every one of the four fibers, the saturated affine
resultant after removing `s4=0` has

```text
raw_resultant_degree_M = 32
boundary_factor_degrees = [ <2, 8> ]
saturated_degree = 16
affine_factor_degrees = [ <4, 1>, <4, 1>, <8, 1> ]
s4 = M^2 - 9/125.
```

The translated order-12 classes occurring are only the two Mumford pairs

```text
D_minus = [x^2 - 2/3*x - 4/3,  -7/3*x - 8/3],
D_plus  = [x^2 + 13/3*x + 26/3,  5/3*x + 16/3].
```

## Quartic branches

Let `a` denote the halving parameter `M` on a quartic branch.  Each quartic
component is zero-dimensional: over `Q(a)`, `gcd(E1,E0)` is linear in `N`, and
`s4(a)=(125*a^2-9)/125` is nonzero.

The three distinct quartic `M`-polynomials are:

```text
qA = M^4 + 12/13*M^3 + 1728/4225*M^2
     + 1512/21125*M + 324/105625

qB = M^4 + 12/13*M^3 + 2214/4225*M^2
     - 3348/21125*M - 1863/105625

qC = M^4 - 84/143*M^3 + 56646/511225*M^2
     - 14364/2556125*M - 3807/12780625.
```

The branch equations are:

```text
D_minus:
  qA(a)=0,  N = -(4225*a^3 + 2730*a^2 + 846*a + 198)/126
  qB(a)=0,  N =  (325*a^2 + 390*a - 51)/120

D_plus:
  qB(a)=0,  N = (-325*a^2 + 90*a - 69)/120
  qC(a)=0,  N = (2556125*a^3 - 1083225*a^2 + 105975*a - 9171)/18720
```

Attachment to the four sign fibers:

```text
(-1/3,-1, 4/3): T1 -> D_minus, T2 -> D_plus
(-1/3, 1, 4/3): T1 -> D_plus,  T2 -> D_minus
( 1/3,-1,-4/3): T1 -> D_minus, T2 -> D_plus
( 1/3, 1,-4/3): T1 -> D_plus,  T2 -> D_minus
```

Each extra translated class also has one saturated degree-8 residual factor;
the run printed it and found no rational roots for it.

## Arithmetic checks

All three quartics are irreducible over `Q` and have no rational roots.

For `qA`:

```text
disc(qA) = -184312513953792/33656858361572265625
         = -2^18*3^15*7^2 / (5^12*13^10)
resolvent = t^3 - 1728/4225*t^2 + 73872/1373125*t
            - 93312/34328125
          = (t - 72/325)*(quadratic irreducible).
```

For `qB`:

```text
disc(qB) = -3944197523094110208/33656858361572265625
         = -2^38*3^15 / (5^12*13^10)
resolvent = t^3 - 2214/4225*t^2 - 104004/1373125*t
            - 1615464/34328125
          = (t - 18/25)*(quadratic irreducible).
```

For `qC`:

```text
disc(qC) = -3944197523094110208/5165516130531461748291015625
         = -2^38*3^15 / (5^12*11^10*13^8)
resolvent = t^3 - 56646/511225*t^2 + 8210484/1827629375*t
            - 397374984/6533775015625
          = (t - 1098/39325)*(quadratic irreducible).
```

Since the components are zero-dimensional with linear `N` over `Q(a)`, there
is no plane curve or elliptic/low-genus model to extract here.  The rational
point question on the quartic components is exactly the rational-root question
for `qA`, `qB`, and `qC`.

## Verdict

No rational halving candidate was found in these quartic saturated branches.
In fact, the quartic components themselves have no rational points: their
`M`-coordinates satisfy irreducible quartics over `Q`, and `N` is then forced
linearly over the corresponding quartic field.  The direct affine search
through height `8` found no saturated rational `(M,N)` points, agreeing with
the exact obstruction.

No exact genus-2 torsion certification was triggered because no rational
candidate reached that stage.  The next proof target, if this fiber is to be
closed completely in a standalone lemma, is to package the same saturation
argument for the two `[16]` O/TR components and the degree-8 residual factors.
