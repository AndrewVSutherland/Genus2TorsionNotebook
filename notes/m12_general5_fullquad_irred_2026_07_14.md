# The irreducible quadratic-`B` chart for `M(12)+5`

## Outcome

The genuinely irreducible quadratic-`B` locus is covered by a root-free
quadratic-algebra parametrization.  It is locally alive, with exact cyclic
order-`60` controls, but no rational example was found.  An integral search
through `|r|,|n| <= 100` leaves three fixed quadratics after primes through
`19`; every one has zero open point already at both `23` and `29`.

## Complete root-free parametrization

Normalize

```text
B=x^2+r*x+n.
```

In `Q[x]/B`, the identity `A^2=q^5` and `gcd(B,q)=1` show that
`u=A/q^2` satisfies `u^2=q`.  Write `u=s*x+t`.  Since `q` is monic,

```text
q=u^2+(1-s^2)B.
```

Put

```text
D=t^2-r*s*t+n*s^2 != 0,
h=(1-s^2)*(-s*x+t-r*s)/(2D),
v=u+B*h.
```

Then `v^2=q (mod B^2)`.  Let `A0=q^2*v mod B^2`; every monic quintic
with the required congruence is uniquely

```text
A=A0+B^2*(x+a).
```

Consequently

```text
S=(A^2-q^5)/B^2
```

is a quintic, and the exact norm identity over the compact `M(12)` family is
equivalent to

```text
S=tau*F_(b,w),       tau in Q^{*2}.
```

This proves completeness on the `deg(B)=2`, `gcd(B,q)=1` chart after
normalizing the leading coefficient of `B`.  No root of `B` is selected, so
irreducible quadratic `B` over `Q` is genuinely included.

## Finite controls

The nonsplit open-match counts are

| `p` | 7 | 11 | 13 | 17 | 19 |
|---:|---:|---:|---:|---:|---:|
| matches | 2 | 10 | 14 | 46 | 42 |
| fixed-`B` masks | 2 | 9 | 14 | 40 | 38 |

At `p=7`, the two nonsplit controls are

```text
(b,w)=(3,5), B=x^2+x+4;
(b,w)=(5,6), B=x^2+6*x+3.
```

For both, direct Magma arithmetic gives `#J(F_7)=60` and

```text
(ord(D5),ord(D12),ord(D5+D12))=(5,12,60).
```

The split count at `p=7` is one, exactly half the ordered-root count two in
the separate root-pair chart, providing an independent normalization check.

## Integral fixed-`B` sieve

Using all-discriminant masks is essential: a quadratic irreducible over `Q`
may split modulo a prime.  Among the `39,228` pairs with `|r|,|n| <= 100` and
nonzero nonsquare rational discriminant, the mask sizes at `7,11,13,17,19`
are `3,17,20,77,89`.  Only

```text
(r,n,disc)=(-83,75,6589), (-15,73,-67), (78,19,6008)
```

survive, with exactly one open point in each tested fiber.  All three fixed
quadratics have zero open points at `p=23`, and again at `p=29`.

This is a bounded, good-open exclusion rather than a global theorem.  A
global point on one of these three quadratics would have to be nonintegral,
bad, or on a removed boundary at `p=23`; points outside the coefficient box
are not covered.

## Reproduction

```text
python3 code/m12_general5_fullquad_irred.py --mode full \
  --primes 7,11,13,17,19 --discriminant nonsplit
python3 code/m12_general5_fullquad_irred.py --mode box \
  --primes 7,11,13,17,19 --bound 100
magma -b MemGB:=2 code/m12_general5_fullquad_irred_local_verify.m
```

The Python runs were capped at 2.4 GB.  The `p=7` census used about 14 MB;
the Magma order checks used about 32 MB.
