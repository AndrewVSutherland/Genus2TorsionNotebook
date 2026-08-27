# New-base search on the degree-2 `[5,5]` contact family

This is the bounded follow-up to `notes/z5x5_degree2_hensel_reconstruction.md`.
It changes the base residue cells instead of lifting the old pair more deeply.
No rational `[5,5]` example was found.

## Files and runs

```text
code/z5x5_degree2_new_base_screen.m
code/z5x5_degree2_new_base_exact_scan.m
data/z5x5_degree2_new_base_screen_p7_p11.txt
data/z5x5_degree2_new_base_exact_scan_adaptive_i30_h6.txt
```

The recorded outputs were produced by:

```bash
magma -b p1:=7 p2:=11 code/z5x5_degree2_new_base_screen.m \
  > data/z5x5_degree2_new_base_screen_p7_p11.txt

timeout 180s magma -b integer_height:=30 rational_height:=6 \
  code/z5x5_degree2_new_base_exact_scan.m \
  > data/z5x5_degree2_new_base_exact_scan_adaptive_i30_h6.txt
```

Both commands completed with exit status zero.  Magma was V2.29-4.

## Exact search box

The family is

```text
f = (1+a*x+b*x^2)^2-k*x^5.
```

The integer phase covers every triple with

```text
|a|, |b|, |k| <= 30.
```

The rational phase covers every reduced coordinate `n/d` with

```text
|n| <= 6, 1 <= d <= 6,
```

and omits only triples already covered by the integer phase.  Thus 226981
integer triples and 101626 additional rational triples were enumerated,
328607 unique triples in total.  The cases `k=0` and `disc(f)=0` were rejected
exactly over `Q`.

The discriminant factorization used to detect singular reductions of the
displayed equation is

```text
k^2 * (16*a^4*b^3 - 128*a^2*b^4 + 108*a^5*k
       + 256*b^5 - 900*a^3*b*k + 2000*a*b^2*k + 3125*k^2).
```

## Actual finite subgroup screen

The screen computes `AbelianGroup(Jacobian(C))` over finite fields and counts
invariant factors divisible by 5.  Hence it tests containment of
`(Z/5 Z)^2`, not merely divisibility of the Jacobian order by 25.

Cached complete residue tables were built at `p=3,7,11,13,17`.  Their smooth
rank-two cell counts were respectively

```text
0, 6, 262, 60, 176.
```

The corresponding numbers of smooth cells were

```text
12, 252, 1120, 1872, 4352.
```

Singular or undefined reductions of the displayed equation are skipped rather
than killed.  If fewer than two smooth displayed reductions remain, the script
tests `p=19,23,29,31,37,41,43` directly.  Every final survivor passed the
actual rank-two test at at least two primes where the displayed equation is
smooth; the insufficient-good-prime count was zero.

The fact that no smooth displayed mod-3 cell has 5-rank two is a useful new
chart constraint: an integral rational `[5,5]` hit represented by this
equation must have singular displayed reduction modulo 3, equivalently the
displayed polynomial discriminant must vanish modulo 3.  This does **not**
prove intrinsic bad reduction of the curve or Jacobian at 3; a different or
minimal integral model could be smooth.  When the displayed model is smooth
modulo 7, a hit must lie in one of only six cells; all six have `b=0 mod 7`
and finite invariants `[10,10]`.

## Simplicity and exact torsion

The adaptive sieve left 84 bases:

```text
integer phase: 62
rational phase: 22
```

For every one, Magma found a good prime at which the resultant of the
Frobenius polynomial with `W-T^12` is irreducible.  This is the standard
geometric-simplicity certificate.  An independent Sage check on the first
survivor `(a,b,k)=(-30,-4,-15)` reproduced the mod-11 Frobenius polynomial

```text
T^4 + 5*T^3 + 18*T^2 + 55*T + 121
```

and verified that its twelfth-power transform is irreducible.

Magma `TorsionSubgroup` was then run on an integral square-scaled model for
all 84 survivors.  The exact distribution was

```text
80 with torsion [5],
 4 with torsion [10],
 0 with torsion [5,5].
```

The marked class `(x,1)` was independently checked in Sage on the first
survivor and has exact order 5.

## Why the Hensel stage is empty

The original workflow proposed enumerating finite contact branches, pairing
them across primes, and Hensel lifting.  Exact `TorsionSubgroup` is a stronger
and much cheaper gate once the rational base is fixed: `[5]` or `[10]`
proves that no independent rational order-5 class, open or boundary, exists
above that base.  Consequently all 84 candidates were eliminated exactly
before branch enumeration, and the number sent to the Hensel stage is zero.

This is a bounded negative result, not a nonexistence theorem for `[5,5]`.
The next justified enlargement is to exploit the forced singular reduction of
the displayed equation modulo 3 and the six smooth displayed mod-7 cells,
rather than deepen lifts over the old mod-11 and mod-19 base pair.  A larger
height scan should parameterize or stratify the displayed discriminant-zero
congruence modulo 3 and continue using adaptive smooth-prime rank tests plus
exact torsion as the final gate.
