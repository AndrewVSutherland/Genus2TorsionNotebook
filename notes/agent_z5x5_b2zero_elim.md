# Agent notes: `Z/5 x Z/5` full Mumford norm, `b2=0` branch

Date: 2026-07-02.

Code:

```text
code/agent_z5x5_b2zero_elim.m
```

Logs:

```text
results/z5x5_b2zero_count_p7.log
results/z5x5_b2zero_slice_sat_p7.log
results/z5x5_b2zero_count_p11.log
```

## Setup

I reconstructed the full normalized norm system

```text
f = h^2 - K*x^5,       h = 1 + h1*x + h2*x^2,
A^2 - B^2*f = U^5,
U = x^2 + s*x + t,
A = x^5 + a4*x^4 + ... + a0,
B = b0 + b1*x + b2*x^2.
```

As before, `a4,...,a0` are forced by the coefficients of degrees
`9,8,7,6,5`.  The full residual sizes are unchanged:

```text
E0..E4 degrees 30,27,24,21,18
terms       852,560,347,208,123
```

On the linear branch `b2=0`, with `B=b0+b1*x`, the residual sizes drop to:

```text
E0..E4 degrees 10,9,8,7,6
terms       51,41,32,24,18
```

The open product used for counts and sliced saturation is

```text
K * b1 * b0 * disc(U) * Res(B,U) * disc(f).
```

Here `b1=0` removes the constant-`B` contact slice, and `b0=0` was also
excluded as the linear-`B` contact-root boundary.  The explicit nontrivial
factors recorded by Magma are:

```text
disc(U) = s^2 - 4*t,
Res(B,U) = -s*b0*b1 + t*b1^2 + b0^2,
disc(f) =
  16*h1^4*h2^3*K^2 - 128*h1^2*h2^4*K^2 + 108*h1^5*K^3
  + 256*h2^5*K^2 - 900*h1^3*h2*K^3 + 2000*h1*h2^2*K^3
  + 3125*K^4.
```

## Commands

```sh
magma -b prime_bound:=7 count_prime_bound:=7 sample_limit:=3 \
  code/agent_z5x5_b2zero_elim.m \
  > results/z5x5_b2zero_count_p7.log 2>&1

magma -b prime_bound:=7 count_prime_bound:=7 sample_limit:=2 \
  do_slice_saturation:=true do_dimension:=true do_primary:=true \
  code/agent_z5x5_b2zero_elim.m \
  > results/z5x5_b2zero_slice_sat_p7.log 2>&1

magma -b prime_bound:=11 count_prime_bound:=11 sample_limit:=3 \
  code/agent_z5x5_b2zero_elim.m \
  > results/z5x5_b2zero_count_p11.log 2>&1
```

I also tried unsliced exact saturation over `F_3`; even with dimension skipped,
it was too slow for this pass and was interrupted.  The sliced saturation below
is the bounded exact algebra computation.

## Modular Evidence

Over `F_3`, the branch has no open rational points:

```text
raw_residual_points=103
open_points=0
```

This is only an `F_3` point-count failure, not a saturation obstruction over
the algebraic closure.

Over `F_7`, the branch is live:

```text
raw_residual_points=2521
open_points=12
```

A representative open point is

```text
h1=1, h2=0, K=2, s=0, t=6, b0=3, b1=1
f = 5*x^5 + x^2 + 2*x + 1
U = x^2 + 6
V = 3
A = x^5 + x^3 + 6*x^2 + 2*x + 6
B = x + 3
#J(F_7)=100
jacobian_rank=5
```

Magma verifies that the contact class and `[U,V]` are both killed by `5` and
are independent over `F_5`.

The slice through this point, `h1=1,h2=0`, has:

```text
raw_residual_points=63
open_points=2
```

Exact sliced saturation by the same open product gives:

```text
slice_raw dimension=2 component_degrees=[1,3]
slice_sat_product dimension=0
slice_sat_product primary_components=11
```

The two open points on this slice are the pair with `B=x+3` and `B=6*x+4`,
both independent.

Over `F_11`, the branch is also live:

```text
raw_residual_points=15281
open_points=240
```

One smooth independent lift target is

```text
h1=1, h2=0, K=7, s=5, t=10, b0=2, b1=4
f = 4*x^5 + x^2 + 2*x + 1
U = x^2 + 5*x + 10
V = 8*x + 6
A = x^5 + 7*x^4 + 10*x^3 + 9*x^2 + 6*x + 5
B = 4*x + 2
#J(F_11)=150
jacobian_rank=5
```

## Verdict

The `b2=0` full Mumford norm branch does **not** die after removing the obvious
open-chart boundaries.  It has smooth open finite-field points over `F_7` and
`F_11`, and the sampled points give independent `5`-torsion classes.

The strongest bounded algebra result is the `F_7`, `h1=1,h2=0` slice: before
saturation it has dimension `2`, and after saturating by
`K*b1*b0*disc(U)*Res(B,U)*disc(f)` it becomes a zero-dimensional chart with
`11` primary components and two open rational points.  This is a good local
chart for a finite-field lift attempt, not a modular obstruction.
