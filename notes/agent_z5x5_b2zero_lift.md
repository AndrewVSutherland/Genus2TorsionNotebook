# Agent notes: `Z/5 x Z/5`, `b2=0`, lift of live finite-field charts

Date: 2026-07-02.

Worker: `Z5_B2ZERO_LIFT`.

Code:

```text
code/agent_z5x5_b2zero_lift.m
```

Logs:

```text
results/z5x5_b2zero_lift_default.log
results/z5x5_b2zero_lift_wide.log
```

## Setup

I used the live open slice found by `agent_z5x5_b2zero_elim`:

```text
h1 = 1, h2 = 0,
f = (1+x)^2 - K*x^5,
U = x^2 + s*x + t,
B = b0 + b1*x.
```

The five slice variables are ordered as

```text
K, s, t, b0, b1.
```

The script reconstructs the same forced-`A` norm equations

```text
A^2 - B^2*f = U^5
```

and clears only powers of `2` from the residuals, so the congruence tests at
`p=7,11` are equivalent to the original equations.  On this slice the cleared
residual sizes are:

```text
E0 degree=10 terms=26
E1 degree=9  terms=23
E2 degree=8  terms=20
E3 degree=7  terms=15
E4 degree=6  terms=11
```

## Commands

Default lift through `7^3` and `11^3`, with a smaller rational box:

```sh
magma \
  code/agent_z5x5_b2zero_lift.m \
  > results/z5x5_b2zero_lift_default.log 2>&1
```

Wider bounded pass used for the verdict:

```sh
magma -b max_7_power:=4 max_11_power:=3 \
  rational_num_bound:=200 rational_den_bound:=100 max_search_tuples:=5000000 \
  code/agent_z5x5_b2zero_lift.m \
  > results/z5x5_b2zero_lift_wide.log 2>&1
```

## Lift Tables

The two `F_7` open points on the slice both have full slice Jacobian rank.
The determinant is `4 mod 7`, so Hensel lifting is nonsingular and unique.

```text
F7_slice_A base = [2,0,6,3,1]
  mod 7:    1 lift, [2,0,6,3,1]
  mod 49:   1 lift, [30,35,13,3,36]
  mod 343:  1 lift, [30,84,258,52,134]
  mod 2401: 1 lift, [2088,2142,258,1424,1849]

F7_slice_B base = [2,0,6,4,6]
  mod 7:    1 lift, [2,0,6,4,6]
  mod 49:   1 lift, [30,35,13,46,13]
  mod 343:  1 lift, [30,84,258,291,209]
  mod 2401: 1 lift, [2088,2142,258,977,552]
```

Balanced representatives modulo `2401` are:

```text
F7_slice_A: [-313,-259,258,-977,-552]
F7_slice_B: [-313,-259,258, 977, 552]
```

The comparison `F_11` point is also nonsingular.  The determinant is
`4 mod 11`.

```text
F11_slice_A base = [7,5,10,2,4]
  mod 11:   1 lift, [7,5,10,2,4]
  mod 121:  1 lift, [84,5,87,101,81]
  mod 1331: 1 lift, [689,489,571,1190,202]
```

Balanced modulo `1331`:

```text
[-642,489,571,-141,202]
```

## Rational Search

I searched rationals `a/b` with

```text
|a| <= 200, 1 <= b <= 100, gcd(b,p)=1,
```

constrained to the final p-adic residue class for each coordinate.  Each
candidate tuple was checked exactly against the cleared rational residuals.

```text
chart          modulus   candidate counts (K,s,t,b0,b1)   tested    exact hits
F7_slice_A     2401      (9,10,10,10,10)                  90000     0
F7_slice_B     2401      (9,10,10,10,10)                  90000     0
F11_slice_A    1331      (17,18,17,17,18)                 1591812   0
```

I also tried CRT-constrained comparison boxes joining each `F_7` lift to the
sample `F_11` lift:

```text
CRT_F7A_F11A modulus 3195731: no coordinate candidates in the same box
CRT_F7B_F11A modulus 3195731: no coordinate candidates in the same box
```

No exact rational tuple appeared, so no rational `f,U,V` certification over
`Q` was triggered.

## Verdict

The live `F_7` slice points are genuinely p-adically smooth: both lift
uniquely through `7^4`, and already through the requested `7^2` and `7^3`.
The sampled `F_11` comparison point likewise lifts uniquely through `11^3`.

Thus these charts are not locally obstructed; they look like isolated smooth
p-adic points of the saturated open chart.  The bounded rational reconstruction
search found no rational candidate of numerator height `<=200` and denominator
`<=100` in the tested congruence classes.
