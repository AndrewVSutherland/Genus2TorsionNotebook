# Agent notes: `Z/5 x Z/5`, `b2=0`, global/local CRT sieve

Date: 2026-07-02.

Worker: `Z5_GLOBAL_SIEVE`.

Code:

```text
code/agent_z5x5_b2zero_global_sieve.m
```

Log:

```text
results/z5x5_b2zero_global_default.log
results/z5x5_b2zero_global_frontier1M.log
```

## Setup

I kept the previous normalized slice

```text
h1 = 1, h2 = 0,
f = (1+x)^2 - K*x^5,
U = x^2 + s*x + t,
B = b0 + b1*x,
variables = (K,s,t,b0,b1).
```

The script reconstructs the forced-`A` equations

```text
A^2 - B^2*f = U^5
```

and clears the same powers of `2` as the lift script:

```text
integral_residual_denominators=[65536, 16384, 16384, 256, 512]
E0..E4 terms=[26,23,20,15,11]
```

It uses the two `F_7` slice points and both sign choices of the `F_11`
slice point:

```text
F7A  = [2,0,6,3,1]
F7B  = [2,0,6,4,6]
F11A = [7,5,10,2,4]
F11B = [7,5,10,9,7]
```

## Command

Main bounded search:

```sh
magma \
  code/agent_z5x5_b2zero_global_sieve.m \
  > results/z5x5_b2zero_global_default.log 2>&1
```

Extra frontier-only checkpoint:

```sh
magma -b do_low_crt_search:=false do_adjacent_probe:=false \
  frontier_height:=1000000 \
  code/agent_z5x5_b2zero_global_sieve.m \
  > results/z5x5_b2zero_global_frontier1M.log 2>&1
```

Default parameters:

```text
low CRT search:   7^4 * 11^3 = 3195731
high CRT frontier: 7^7 * 11^6 = 1458956660623
frontier height: 500000
search boxes:
  |num|<=5000,  den<=5000
  |num|<=20000, den<=1000
  |num|<=1000,  den<=20000
```

The main run took about `38s`; the frontier-only checkpoint took about `37s`.

## Lift Data

All four selected slice charts are nonsingular and lift uniquely to the high
moduli.

```text
F7A  mod 7^7:  [-379671, -86695, -350288,  162291,  225142]
F7B  mod 7^7:  [-379671, -86695, -350288, -162291, -225142]
F11A mod 11^6: [ -72516, 218773,   48487,  445744,  784161]
F11B mod 11^6: [ -72516, 218773,   48487, -445744, -784161]
```

The low CRT balanced residues are:

```text
F7A/F11A: [-1258437, 1082592, 21867, -106621,   347593]
F7A/F11B: [-1258437, 1082592, 21867,  515238, -1313899]
F7B/F11A: [-1258437, 1082592, 21867, -515238,  1313899]
F7B/F11B: [-1258437, 1082592, 21867,  106621,  -347593]
```

## Rational CRT Search

For each of the four sign combinations and each of the three boxes above, the
script generated all coordinate rationals satisfying the low CRT class, then
sieved the resulting tuple product at two additional primes larger than the
denominator bound.

Candidate tuple counts:

```text
box |num|<=5000, den<=5000:
  4536, 4032, 4032, 4536

box |num|<=20000, den<=1000:
  8064, 12096, 12096, 8064

box |num|<=1000, den<=20000:
  14112, 18432, 18432, 14112
```

Total tuple products considered: `122544`.

In all twelve searches the first added good prime already killed every tuple:

```text
modular_survivors_after_prefix_primes=[0,0]
exact_tested=0
exact_hits=0
```

So no rational tuple was found, and no `f,U,V,A,B` or independent
`Z/5 x Z/5` certificate over `Q` was triggered.

## High CRT Height Frontier

Using the unique high lifts and CRT modulus

```text
M = 7^7 * 11^6 = 1458956660623,
```

the main script scanned denominators `1..500000` for each coordinate and each
sign combination.  For all four sign combinations and all five coordinates:

```text
min_height_le_500000=none
tuple_height_lower_bound_from_coordinates=>500000
```

The extra frontier-only checkpoint scanned denominators `1..1000000`.  In
every sign combination, the coordinates `K`, `s`, and `b1` still had no
rational representative of naive height at most `1000000`; `t` first appeared
at height `861780`, and the same-sign `b0` coordinates first appeared at
height `699357`.

```text
K:  none <= 1000000
s:  none <= 1000000
t:  -682447/861780, height 861780
b0: +/-699357/658586, height 699357 in the same-sign cases;
    none <= 1000000 in the mixed-sign cases
b1: none <= 1000000
tuple_height_lower_bound_from_coordinates=>1000000
```

Thus, any rational point integral at `7` and `11` and lying simultaneously in
one of these selected `F_7/F_11` smooth charts must have at least one of
`K,s,t,b0,b1` of naive height greater than `1000000`.  This is already a
strong height obstruction from CRT alone, before using the equations over `Q`.

## Adjacent Local Probe

I also lifted one adjacent open point from each count log on fixed
`(h1,h2)` slices:

```text
F7 adjacent h1=2,h2=0:
  base [1,0,5,1,3]
  unique lift to 7^3 = [274,42,236,127,206]
  jacobian determinant = 4 mod 7

F11 adjacent h1=1,h2=1:
  base [4,3,2,2,10]
  unique lift to 11^3 = [1280,850,1135,134,131]
  jacobian determinant = 5 mod 11
```

These adjacent samples are also smooth local points, not cheap local
obstructions.

## Verdict

The selected `h1=1,h2=0` `F_7/F_11` smooth charts remain locally alive, but the
combined CRT constraints are very restrictive globally.  The bounded rational
CRT search found no tuple in the three tested denominator/height partitions,
and the high-power CRT classes alone exclude tuple coordinate height
`<=1000000` for all four sign combinations.

This does not eliminate the whole `b2=0` branch over `Q`, because it only uses
the selected smooth local charts.  It does, however, push any rational point
passing through these charts beyond a very large naive-height threshold.
