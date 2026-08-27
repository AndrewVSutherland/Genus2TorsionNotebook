# A(12) to A(2,24): halving translated order-12 classes

This note turns the direct `A(12)` / `A(2,12)` negative search into a
structural halving-cover statement.  The relevant existing files are
`torsion_goal_log.md`, `search_A2_24_from_A12.m`,
`certify_A2_24_A12_survivor_exact_halves.m`, `a12_parameterization.tex`, and
`a2_12_resolvent.tex`.

## Existing search state

The rational `A(12)` chart has

```text
C: Y^2 = R F,        F = R x^2 + 4(R+x-1)(R-1) = Q^2 + R ell^2,
```

with `R,Q` quadratic and `ell` linear.  The visible order-12 class is

```text
P12 = P4 + P6,
P4 = [Q, R ell],       P6 = [R+x-1, xR].
```

The `A(2,12)` step asks the residual quartic `F` to acquire an extra rational
2-torsion class.  The recorded cubic resolvent over `Q(p,z,r)` is irreducible;
after fixing small `z`, the numerator degrees in `r` remain large, so this did
not give a low-degree rational parametrization.

The corrected search `search_A2_24_from_A12.m` tests translated classes
`P12 + T`, where `T` is a rational 2-torsion class available to the script.  The
height-10 run checked `7134` translated order-12 classes:

```text
local_half_reject=7104, local_half_pass=30,
half_candidates=0, target_hits=0.
```

The only height-10 local passes printed with `LocalPrimeBound:=1000` are
duplicate A12 presentations of one curve,

```text
f = 81/15625*x^6 + 369/5000*x^5 + 2499/10000*x^4
    - 1/18*x^3 - 151/648*x^2 - 625/648*x + 6875/11664.
```

The four distinct printed translated order-12 classes are:

```text
A: u = x^2 - 25/6*x + 625/162,
   v = 32/9*x - 400/81

B: u = x^2 + 1150/693*x + 3125/6237,
   v = -2949/5929*x + 13775/17787

C: u = x^2 + 50/9*x - 625/27,
   v = 4/9*x + 100/27

D: u = x^2 + 1100/189*x - 625/1701,
   v = 1601/7938*x - 52825/71442
```

The exact certifier finds torsion invariants `[2,2,12]` on an integral square
model and no rational halves of these four classes.

## Generic halving-cover equations

Let `D_T = P12 + T = [u_T,v_T]` be any translated order-12 class on a smooth
member of the `A(2,12)` cover.  On the standard affine Mumford halving chart,
introduce two parameters `M,N` and put

```text
ell_T = v_T + u_T*(M*x + N),
S_T   = (ell_T^2 - f)/u_T
      = s4*x^4 + s3*x^3 + s2*x^2 + s1*x + s0.
```

The residual divisor is double precisely when `S_T` is a scalar multiple of a
square quadratic.  On the open chart `s4 != 0`, this is equivalent to

```text
E1 = 8*s4^2*s1 - s3*(4*s4*s2 - s3^2) = 0,
E0 = 64*s4^3*s0 - (4*s4*s2 - s3^2)^2 = 0.
```

Then the half is

```text
H = [G, -ell_T mod G],
G = x^2 + (s3/(2*s4))*x + (4*s4*s2 - s3^2)/(8*s4^2).
```

Thus the structural `A(2,24)` cover is obtained by imposing `E1=E0=0` for the
translated order-12 classes and saturating by the chart boundary `s4=0`, along
with the usual smoothness and exact-order open conditions.  Since multiplication
by 2 on a genus-2 Jacobian has geometric degree 16, the saturated fiber over a
fixed order-12 class should have degree 16 over an algebraic closure.

## Survivor component behavior

For the known local survivor, the exact equations for each of the four classes
have:

```text
deg(E1) = 6,     deg_M(E1)=6, deg_N(E1)=3,
deg(E0) = 8,     deg_M(E0)=8, deg_N(E0)=4,
gcd(E1,E0)=1,    dim <E1,E0> = 0.
```

The certifier resultants have degree 32 and factor as

```text
two rational linear factors, each with multiplicity 8,
two irreducible degree-8 factors.
```

The new probe `code/agent_A2_24_halving_cover_probe.m` explains the linear
factors: for all four classes,

```text
s4 = M^2 - 81/15625 = (M - 9/125)(M + 9/125).
```

So the repeated rational linear factors are exactly the `s4=0` boundary, not
honest affine square-quartic halves.  The boundary rational points are:

```text
A: (M,N)=(-9/125,-13/16),       ( 9/125, 13/16)
B: (M,N)=(-9/125,-2421/6160),   ( 9/125, 2421/6160)
C: (M,N)=(-9/125,-9/80),        ( 9/125, 9/80)
D: (M,N)=(-9/125,-157/1680),    ( 9/125, 157/1680)
```

At every one of these boundary points, `deg(S)=2` and the square-quartic test is
false.  After removing the `s4=0` boundary, each class has affine degree 16,
with factor degrees

```text
[ <8,1>, <8,1> ].
```

This matches the expected 16 geometric halves, organized here as two degree-8
Galois orbits.  There is no rational component and no rational half.

## Probe and next decision

Implemented symbolic probe:

```text
code/agent_A2_24_halving_cover_probe.m
```

Run result:

```text
torsion_invariants_integral_square_model=[ 2, 2, 12 ]
contains_Z2_Z24=false
for A,B,C,D:
  boundary_s4_rational_points=2
  boundary points have deg(S)=2 and square_quartic=false
  resultant_N_affine_degree_after_s4=16
  affine_factor_degrees=[ <8,1>, <8,1> ]
```

Outcome: the known local survivor is not evidence for a rational halving
component; the only rational-looking pieces of the raw resultant are boundary
artifacts.  A Chabauty/descent object should not be derived from the unsaturated
degree-32 resultants.  If this route continues, the next structural target is
the saturated global cover over a chosen `A(2,12)` resolvent branch, with the
`s4=0` boundary removed first; only if that produces a low-dimensional quotient
or a stable degree-8 component should a Chabauty/descent calculation be built.
