# Full extra-Weierstrass `M(12)` halving route

This leaves the one-parameter line `a=(1-r)/4` and works on the full two-dimensional surface where `T+1` is split and `Q4=W/(T+1)` has an additional rational root.

Start with

```text
T = a x^2 - x + r,
W = ((x-r)(T+1))^2 + 4a x^2T(T+1).
```

Split `T+1` by

```text
a = (1-z^2)/(4(r+1)).
```

Then a root `u` of

```text
Q4 = W/(T+1)
```

gives an independent rational 2-torsion candidate.  For a root `w` of `T+1`, move `w` to infinity via

```text
X = 1/(x-w),   Y_old = Y_new/X^3.
```

The extra 2-torsion class is represented in the odd model by

```text
beta = 1/(u-w).
```

## Algebraic `[4,12]` conditions

The point of this route is not just to search the surface, but to write the
condition that the independent 2-torsion class is divisible by `2`.

After moving `w` to infinity, write the odd quintic as

```text
f5(X) = (X - beta) g(X),    g(X) = g4 X^4 + g3 X^3 + g2 X^2 + g1 X + g0.
```

A half of the class `[beta - infinity]` exists exactly when there are
parameters `m,n,A,B` such that

```text
g(X) - (X - beta)(mX+n)^2 = g4 (X^2 + AX + B)^2.
```

Thus the algebraic `[4,12]` locus on this chart is cut out by

```text
Q4(u) = 0
H0 = H1 = H2 = H3 = 0,
```

where `Hi` are the four non-leading coefficient equations in the identity
above.  The raw equations live in `r,z,u,A,B,m,n`.

Equivalently, the top two coefficient equations solve for

```text
A = (g3 - m^2)/(2*g4),
B = (g2 - 2mn + beta*m^2 - g4*A^2)/(2*g4),
```

and the remaining conditions are

```text
R1 = g1 - n^2 + 2 beta mn - 2 g4 A B = 0,
R0 = g0 + beta n^2 - g4 B^2 = 0.
```

So a compact construction problem is:

```text
Q4(u) = 0,    R1(r,z,u,m,n) = 0,    R0(r,z,u,m,n) = 0,
```

with `A,B` recovered from the displayed formulas.  The chart exclusions are
`r != -1`, `z != -1` for this chosen root of `T+1`, `z^2 != 1`, `z != 0`,
`u != w`, and the usual smoothness/nondegeneracy conditions.

The script

```text
code/m12_full_surface_z12x4_equations.m
```

generates these equations.  Its summary is in

```text
data/m12_full_surface_z12x4_equations_summary.txt
```

and the explicit polynomial equations are in

```text
data/m12_full_surface_z12x4_equations.txt
```

## Reduced `s=m^2` condition

The two remaining halving equations can be reduced further without a full
Groebner elimination.  Put

```text
s = m^2,   t = mn,   N = n^2.
```

After solving for `A,B`, the equations are linear in `t,N`; the only
nonlinear relation left is

```text
t^2 = sN.
```

Writing

```text
A = (g3 - s)/(2g4),
D = g2 + beta*s - g4*A^2,
P = 2(beta + A),
Q0 = g1 - A*D,
L = 8(g4*beta - s)(beta + A) + 4D,
C = 4(g4*beta - s)Q0 + 4g4*g0 - D^2,
```

the generic branch `L != 0` has

```text
t = -C/L,   N = P*t + Q0.
```

The single compatibility equation is

```text
F_s = C^2 + s*P*C*L - s*Q0*L^2 = 0.
```

Thus, away from the special branch `L=C=0`, the `[4,12]` construction problem
is:

```text
Q4(u) = 0,   F_s(r,z,u,s) = 0,
```

with `s` required to be a rational square for an actual rational `m`.

Reducing `F_s` modulo `Q4(u)` gives

```text
F_s mod Q4 = (r+1)^3 * G_s.
```

Since `r=-1` is excluded on this chart, the main nonboundary condition is

```text
Q4(u) = 0,   G_s(r,z,u,s) = 0.
```

The script

```text
code/m12_full_surface_z12x4_s_condition.m
```

generates this reduction.  The degree summary is:

```text
F_s total_degree 64, degree_s 8, terms 27785
F_s_mod_Q4 total_degree 76, degree_u 3, degree_s 8, terms 29349
G_s_main total_degree 73, degree_u 3, degree_s 8, terms 24865
```

The summary and explicit equations are in

```text
data/m12_full_surface_z12x4_s_condition_summary.txt
data/m12_full_surface_z12x4_s_condition.txt
```

## Constructive `s`-condition search

The script

```text
code/m12_full_surface_z12x4_s_construct_search.m
```

uses the reduced equations constructively.  For each rational `(r,z)` it:

```text
1. finds rational roots u of Q4,
2. specializes the generic condition to a univariate polynomial F_s(s),
3. keeps rational square roots s=m^2,
4. recovers t=mn and n,
5. constructs the Mumford half
   [X^2 + AX + B, (X-beta)(mX+n) mod (X^2 + AX + B)],
6. verifies 2H = [beta - infinity].
```

It also handles the special branch `L=C=0` separately by solving

```text
t^2 - s P t - s Q0 = 0.
```

A complete height-20 run gave:

```text
checked_split 258548
extra_roots 1640
order12_tests 3280
s_polys 3280
s_roots 0
square_s 0
zero_s_polys 0
special_s_roots 0
special_square_s 0
constructed 0
verified_halves 0
hits 0
```

So, up to height `20`, every independent extra 2-torsion candidate produces a
nonzero degree-8 `s` polynomial with no rational root; the exceptional
`L=C=0` branch also contributes no rational `s` values.

The run outputs are in

```text
data/m12_full_surface_z12x4_s_construct_h10.txt
data/m12_full_surface_z12x4_s_construct_h15.txt
data/m12_full_surface_z12x4_s_construct_h20.txt
```

## Exact search

The script

```text
code/m12_full_surface_z12x4_search.m
```

searches rational `r,z` and exact rational roots `u` of `Q4`, then uses Magma's exact `IsDivisibleBy` on an integral model.

A height-15 run gave

```text
checked_split 80924
extra_roots 864
exact_tests 1728
hits 0
```

## Finite-field sieve

The script

```text
code/m12_full_surface_z12x4_finite_field_sieve.m
```

checks the full surface over finite fields.  It found:

```text
p 5 good_rz_with_qroot 2 qroots 2 classes 4 divisible 4 bad 8
p 7 good_rz_with_qroot 4 qroots 6 classes 12 divisible 0 bad 14
NO GOOD FULL-SURFACE HALVES MOD 7
```

Thus, on the full two-dimensional extra-Weierstrass surface, there are good mod-`7` points with extra roots, but none has the independent 2-torsion class divisible by `2` in the finite-field Jacobian.

Conclusion: any rational `Z/12Z x Z/4Z` example on this full surface must reduce to the boundary modulo `7`.  The nonboundary part of the full surface is obstructed modulo `7`.


## Boundary residues

The script

```text
code/m12_full_surface_z12x4_boundary_sieve.m
```

records the mod-`7` boundary classification and then applies the same residue test at several auxiliary primes before doing exact rational halving.

Modulo `7`, the affine `(r,z)` residues split as follows:

```text
r=-1                  7
z=+-1                12
z=0                   6
disc(W)=0             8
good, no Q4 root     12
good Q4 root, no half 4
good Q4 root, half    0
```

The four good residues with a `Q4` root are

```text
(3,3), (3,4), (4,3), (4,4),
```

and every resulting independent 2-torsion class is nondivisible by `2` in the finite-field Jacobian.

A height-30 run with primes

```text
7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43
```

gave

```text
DONE height 30 total 1234321 survivors 5587 exact 0 hits 0
param_boundary 4441
W_bad 1142
T_bad 0
Q4_no_root 4
odd_bad 0
order12_bad 0
no_independent 0
exact_tested 0
```

So, after the local filters, every small-height survivor either lies on the explicit parameter boundary (`r=-1`, `z^2=1`, or `z=0`), is singular, or has no rational extra root of `Q4`.  No nonsingular survivor reaches an exact rational halving test.

This is not yet a complete p-adic proof along every boundary branch.  It does, however, make the computational picture strongly negative for this `Z/12Z x Z/4Z` route on the full extra-Weierstrass surface.


## Height-50 boundary extension

For the simple-target `[4,12]` question, I extended the same boundary sieve to
height `50`:

```text
magma -b height:=50 progress_interval:=1000000 \
    code/m12_full_surface_z12x4_boundary_sieve.m \
    > data/m12_full_surface_z12x4_boundary_h50.txt
```

The result was still negative:

```text
DONE height 50 total 9579025 survivors 15591 exact 0 hits 0
param_boundary 12377
W_bad          3166
T_bad          0
Q4_no_root     48
odd_bad        0
order12_bad    0
no_independent 0
exact_tested   0
```

Thus the local boundary filters leave more residue survivors at height `50`,
but none reaches the actual rational halving test: after exact rational
inspection every survivor is still on the parameter boundary, singular, or has
no rational extra root of `Q4`.  No simple `[4,12]` example was found on this
full extra-Weierstrass `M(12)` surface.
