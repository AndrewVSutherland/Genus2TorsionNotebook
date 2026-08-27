# Agent notes: full degree-2 Mumford norm on one-contact-5 family

Date: 2026-07-02.

This is the third pass for the `Z/5 x Z/5` lane.  I did **not** repeat the
old `B=1` degree-2 contact slice.  The target was the full norm/Cantor system

```text
f = h^2 - K*x^5,
A(x)^2 - B(x)^2*f = Lambda*U(x)^5,
U = x^2 + s*x + t,   deg(A)<=5, deg(B)<=2.
```

Code:

```text
code/agent_z5x5_full_mumford_norm.m
```

## Normal form

On the open contact chart `h(0) != 0`, scale `y` so

```text
h = 1 + h1*x + h2*x^2.
```

On the generic norm chart, scale the function `A+B*y` so that `A` is monic of
degree `5`; then `Lambda=1`.  Write

```text
U = x^2 + s*x + t,
B = b0 + b1*x + b2*x^2,
A = x^5 + a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0.
```

The coefficients `a4,...,a0` are forced recursively by canceling degrees
`9,8,7,6,5` in

```text
E = A^2 - B^2*(h^2 - K*x^5) - U^5.
```

The remaining full norm equations are

```text
E_i = coeff_x^i(E) = 0,    i=0,...,4.
```

In the eight variables `(h1,h2,K,s,t,b0,b1,b2)`, their sizes are:

```text
E0: degree 30, terms 852
E1: degree 27, terms 560
E2: degree 24, terms 347
E3: degree 21, terms 208
E4: degree 18, terms 123
```

Immediate boundaries:

```text
B = 0       => A^2 = U^5, so U is a double-root boundary.
B constant => contact-type slice after rescaling; skipped here.
```

The lower linear branch `b2=0` is much smaller but still not trivial:

```text
L0: degree 10, terms 51
L1: degree 9,  terms 41
L2: degree 8,  terms 32
L3: degree 7,  terms 24
L4: degree 6,  terms 18
```

## Bounded rational probe

Command:

```sh
magma -b do_symbolic:=false rational_height:=2 prime_bound:=31 finite_trials:=100000 cantor_trials:=10000 code/agent_z5x5_full_mumford_norm.m
```

The normalized integer box

```text
h1,h2,K,s,t,b0,b1,b2 in [-2,2],
K != 0,
(b1,b2) != (0,0)
```

gave:

```text
rational_tested=300000
independent_hits=0
bad_U_boundary=128
bad_curve=0
bad_B_mod_U=0
dependent_norm=0
other_boundary=0
```

All rational norm solutions in this box were the bad `U=x^2` boundary
(`s=t=0`).  There was no smooth degree-2 class and no rational
`Z/5 x Z/5` candidate.

## Finite-field full norm hits

The same run found a genuine full-norm point over `F_11` with quadratic `B`:

```text
p=11
h = 2*x^2 + 7*x + 1
K = 6
f = 5*x^5 + 4*x^4 + 6*x^3 + 9*x^2 + 3*x + 1
U = x^2 + 9*x + 2
V = 10*x
A = x^5 + 2*x^4 + 7*x^3 + 4*x^2 + x + 2
B = 4*x^2 + 2*x + 7
#J(F_11) = 100
```

Magma verified the recovered Mumford class `[U,V]`, the contact class at
`x=0`, `5*D0=5*D2=0`, and no nontrivial `F_5`-linear relation between them.

An earlier short run also found:

```text
p=7
h = 4*x + 1
K = 4
f = 3*x^5 + 2*x^2 + x + 1
U = x^2 + 3
V = 5*x
A = x^5 + 6*x^4 + 3*x^3 + 4*x^2 + 6*x + 4
B = 2*x^2 + 4*x + 2
#J(F_7) = 100
```

Again `deg(B)=2`, so this is outside the discarded `B=1` slice.

## Cantor sanity check

The script also samples the equivalent condition `[5][U,V]=0` directly on the
one-contact family.  In the wider run it found:

```text
p=7
h = 3*x + 1
K = 3
f = 4*x^5 + 2*x^2 + 6*x + 1
U = x^2 + x + 4
V = x + 2
#J(F_7) = 100
```

This confirms the norm formulation and the Cantor formulation are seeing the
same live finite-field geometry.

## Verdict

The full degree-2 Mumford norm system **does not immediately collapse**.  Unlike
the `B=1` contact slice, it has genuine finite-field points with independent
`5`-torsion and `deg(B)=2`.

The tiny rational normalized box is cold: the only norm solutions there are
the double-root `U=x^2` boundary.  The next useful attack is therefore not to
discard the full norm system, but to exploit the five residual equations above:
either eliminate/saturate on the open chart

```text
K != 0, disc(f) != 0, disc(U) != 0, gcd(B,U)=1, b1 or b2 != 0,
```

or use the finite-field hits over `F_7`/`F_11` as local charts for a focused
lifting or modular obstruction.
