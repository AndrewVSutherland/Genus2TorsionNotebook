# The known split fiber at `s=59/49`

The first-radicand elliptic quotient has a particularly natural
specialization at

```text
s = 59/49.
```

It contains the known labelled CK point

```text
(r1,r2,r3,r4,r5) = (1,8/7,-9/7,-5/7,-1/7),
```

which is a useful positive control: the complementary cubic splits and
the first literal radicand is a square, but the other three are not.

The reproducible calculation is

```text
magma -b code/elkies22210_orbit12_seed_fiber_search.m
```

and its compact output is recorded in
`data/elkies22210_orbit12_seed_fiber_b250_summary.txt`.

## Exact Mordell--Weil group

The specialized quotient is

```text
y^2 = x^3-(7103/2401)x^2-(1503792/117649)x
        +162409536/5764801.
```

Its minimal model is

```text
y^2+xy = x^3-5656459x+2951658593.
```

Magma proves rank bounds `[2,2]` and returns Mordell--Weil invariants
`[4,0,0]`.  Thus

```text
E(Q) = Z/4 + Z^2.
```

For a reproducible search box, the script pins the basis

```text
T  = (0,12744/2401),                 order 4,
P1 = (-288/343,101160/16807),
P2 = (-944/441,367688/64827).
```

Their images on the minimal model are

```text
P1 = (-1096,89063),
P2 = (-16892/9,2277427/27).
```

The script certifies that the free-coordinate change from the full
Mordell--Weil group has determinant `+-1`.  The known seed point is

```text
(1593/343,38232/16807) = 2T+P1-P2.
```

## Search and result

We exhaustively searched

```text
n1*P1+n2*P2+k*T,
-250 <= n1,n2 <= 250,  0 <= k < 4.
```

This is `1,004,004` distinct Mordell--Weil elements, or `502,002`
representatives after identifying `P` and `-P` and deleting the
identity.  To avoid constructing enormous rational coordinates for all
of them, the script first tests complete splitting of the complementary
cubic modulo 18 good primes.  This is a rigorous necessary condition:
at primes where `q` has finite reduction, a monic cubic split over `Q`
must remain completely split modulo `p`; indeterminate reductions are
kept rather than rejected.

Only five representatives survive all local tests.  Exact arithmetic on
all five finds only two completely split cubics:

```text
q = 59/49: roots 0,-1,-59/49  (boundary),
q = 8/7:   roots -1/7,-5/7,-9/7 (the known smooth seed).
```

For the smooth seed the exact mask is again

```text
[true,false,false,false].
```

There is no full four-radicand hit in the box, and no second smooth
split lift at all.  This is a strong negative result for this very
natural fiber, though not a proof about all of its infinitely many
rational points.

## Larger independent basis box

The supplementary script

```text
code/elkies22210_orbit12_seed_fiber_mw_sieve.m
```

uses a second certified full Mordell--Weil basis and 40 good primes.  The
completed command

```text
magma -b coefficient_bound:=500 \
  code/elkies22210_orbit12_seed_fiber_mw_sieve.m
```

exhausts `4,008,004` elements with `|n1|,|n2| <= 500`.  The local
split-cubic sieve leaves exactly six elements: the four torsion points
and the two signs of the known seed.  Exact reconstruction therefore
again gives only the boundary value `q=59/49` and the known smooth value
`q=8/7`, with no full hit.  The compact run record is
`data/elkies22210_orbit12_seed_fiber_mw_sieve_b500_summary.txt`.
