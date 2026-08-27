# Exact contact-30 C3-root Richelot verifier

## Construction

Write

```text
A0 = h6-(x-1)^3,
C3 = h6+(x-1)^3.
```

At a rational root `rho` of `C3`, put

```text
B0 = C3/(x-rho),   L = x-rho.
```

Then the contact-30 curve is

```text
y^2 = A0*B0*L
```

and the three classes `[A0,0]`, `[B0,0]`, and `[L,0]` form a
pointwise-rational maximal isotropic Richelot kernel.  Define

```text
U = [B0,L],   V = [L,A0],   W = [A0,B0],
Delta = det(coeff(A0),coeff(B0),coeff(L)),
```

where `[F,G]=F'G-FG'`.  The distinguished codomain is

```text
y^2 = Delta*U*V*W.
```

The new verifier reconstructs either an `(R,branch,rho)` point or a compact
`(u,s,rho)` point, checks both contact identities, verifies the source
orders `5,6,30`, constructs the source and dual kernels, tests all three
dual classes with `IsDivisibleBy`, computes the exact source and codomain
torsion subgroups, extracts an exact order-60 point when one exists, and
applies the D4/root-power geometric-simplicity certificate.

## Discriminant reduction

For general quadratic `A0,B0` and linear `L`, direct symbolic calculation
gives

```text
Res(U,V) = Delta^2,
Res(U,W) = Delta^2*disc(B0),
Res(V,W) = Delta^2*disc(A0).
```

Consequently the norm conditions for halving `[U,0]`, `[V,0]`, and
`[W,0]` are respectively

```text
disc(B0) square,
disc(A0) square,
disc(A0)*disc(B0) square.
```

These are necessary, not sufficient.

There is a stronger condition for `W`.  Let `phi:J -> J'` be the Richelot
isogeny and `phihat` its dual.  If `2H=T` for a nonzero dual-kernel class,
then

```text
S = phihat(H) in J[2](Q),    phi(S)=2H=T.
```

Thus `J[2](Q)` must strictly contain the source Richelot kernel.  For the
factorization `(2,2,1,1)`, this requires at least one of `A0,B0` to split.
Combining this with the product norm condition shows that `W` can halve only
when **both** quadratics split.  The useful necessary filters are therefore

```text
U: B0 splits,
V: A0 splits,
W: A0 and B0 both split.
```

In particular, the hoped-for `W` shortcut from two correlated nonsquare
discriminants does not exist.

For a definitive test, the script also constructs the exact affine Mumford
halving equations.  For a dual factor `R`, set

```text
ell = R*(M*x+N),    S=(ell^2-g)/R.
```

Writing `S=s4*x^4+...+s0`, the ordinary residual chart is

```text
8*s4^2*s1 - s3*(4*s4*s2-s3^2) = 0,
64*s4^3*s0 - (4*s4*s2-s3^2)^2 = 0.
```

The chart-independent certification remains Magma's exact
`IsDivisibleBy([R,0],2)`.

## Odd torsion transport

If `D=D5+D6` has order 30, then `2D` has order 15.  The Richelot kernel is
2-primary, so the degree-4 isogeny is injective on `<2D>`.  The codomain
therefore retains a rational order-15 subgroup.  Its dual kernel is
`(Z/2)^2`; if any nonzero dual class has a rational half, the codomain
contains

```text
Z/2 x Z/60.
```

The exact torsion computation decides whether this is exactly `[2,60]` or
a larger group.

## Positive controls

The self-test proves the symbolic identities and exercises all three exact
halving branches over finite fields:

```text
p=11, R=9, branch=-1, rho=2:
    halves [U,V,W] = [false,true,false],  J'=[2,60]

p=17, R=6, branch=-1, rho=8:
    halves [U,V,W] = [true,false,false],  J'=[2,120]

p=47, R=8, branch=-1, rho=18:
    source factors = [1,1,1,1,1], J=[2,2,2,300]
    split gates = [true,true,true]
    halves [U,V,W] = [false,false,true],  J'=[2,2,600]
```

It also uses a rational synthetic Richelot example with torsion `[2,4]` to
verify four reconstructed Mumford halves against exact Jacobian doubling.
The complete self-test finishes with `SELFTEST_OK`.

Run it with

```text
magma -b mode:="selftest" \
  code/contact30_c3root_richelot_verify.m
```

The finite controls establish that the construction and all three halving
branches are locally viable; they do not assert a rational contact-30 root
point.  The separate root search through height 25,000 found none, so no
global candidate currently reaches the expensive torsion stage.

## Files

```text
code/contact30_c3root_richelot_verify.m
notes/contact30_c3root_richelot_verification_2026_07_11.md
```
