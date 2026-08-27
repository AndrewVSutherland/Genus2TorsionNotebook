# Agent Z/48: A(16) plus rational 3-torsion

Date: 2026-07-02.

Task lane: test the route "simple Z/16 from the A(8)->A(16) cover plus
independent rational 3-torsion" for the missing `Z/48`.

## Context read

I read:

- `a8_order16_cover_status.md`
- `a8_a16_dimension_analysis.md`
- `torsion_goal_log.md`, especially the 2026-07-01 A(16) simple `[16]`
  discovery and follow-up
- relevant `torsion_jac` plus-3/contact notes and scripts:
  `m2228_three_torsion_summary.txt`, `m3222_plus3.md`,
  `m3222_plus3_search.m`, `contact6_m36.md`,
  `contact6_m36_symbolic.m`, `m212_three_torsion.md`,
  and `m212_extra3_search.m`.

The A(8)->A(16) material gives a certified simple `[16]` curve on
`(r,t)=(3,1/3)` and additional small certified `[16]` slices
`(-1,1/3)`, `(-1,1/2)`, `(3,1/2)`.  All recorded exact torsion invariants are
`[16]`; the known first slice has factor type `[1,4]`, and the newer small
slices have factor type `[2,4]`.

## Exact criterion

Let `C/Q` be a smooth genus-2 curve whose Jacobian `J` has a verified rational
point `P` of exact order `16`.

Then `J(Q)` contains a cyclic subgroup `Z/48` if and only if `J(Q)[3]` is
nonzero.

Reason: if `Q` is a rational point of order `3`, then `P+Q` has order
`lcm(16,3)=48`.  Conversely, a rational point of order `48` has nonzero
3-primary part.  So the word "independent" is automatic here because the
2-primary and 3-primary parts have coprime orders.

For an odd-degree model, an exact algebraic 3-torsion condition can be imposed
by the usual cubic-contact test: find a reduced Mumford divisor with
`q` of degree at most `2` and a cubic `h` such that

```text
h(x)^2 - f(x) = m^2 q(x)^3,
```

with the nonboundary conditions `q` squarefree enough and `gcd(q,f)=1`.
This is the same mechanism used in the existing `m2228_three_torsion` and
`contact6_m36` scripts.

## Fast necessary filter

For any good prime `p != 3`, reduction injects rational `3`-torsion:

```text
J(Q)[3] -> J(F_p).
```

Thus a necessary condition for this route is

```text
3 divides #J(F_p)
```

for every good prime `p != 3`.  For an already certified A(16) curve this is
the relevant small-prime obstruction.  Equivalently, the gcd of `#J(F_p)` over
good primes must remain divisible by `3`.

## Probe

I added:

```text
code/agent_Z48_A16_plus3_probe.m
```

It reconstructs the A(8) model from certified A(16) tuples, verifies the
square-root halving relation and exact order-16 relation, computes exact
torsion when requested, and then computes `#J(F_p)` for good primes up to
`97`.

Run:

```text
magma code/agent_Z48_A16_plus3_probe.m
```

Result summary:

| A16 tuple | exact torsion | point-count gcd | first good prime killing 3 |
|---|---:|---:|---|
| `(3,1/3,p=2,mu=9)` | `[16]` | `16` | `p=5`, `#J(F_5)=16` |
| `(3,1/3,p=2,mu=27/11)` | `[16]` | `16` | same curve, `p=5` |
| `(3,1/3,p=34/9)` | `[16]` | `16` | `p=13`, `#J(F_13)=128` |
| `(-1,1/3,p=1/3)` | `[16]` | `32` | `p=5`, `#J(F_5)=32` |
| `(-1,1/2,p=-35/6)` | `[16]` | `16` | `p=5`, `#J(F_5)=16` |
| `(3,1/2,p=17/6)` | `[16]` | `16` | `p=5`, `#J(F_5)=16` |

The first two rows are two A16 presentations of the same underlying curve.

## Outcome

The concrete known/simple A16 examples have immediate small-prime
obstructions to adding rational 3-torsion.  In every distinct curve tested,
the gcd of the good-prime point counts is not divisible by `3`, so
`J(Q)[3]=0` and this route does not produce `Z/48` from the current A16
examples.

This is not a proof that the full A(16)+3 locus is empty.  It says the present
small rational A16 points/slices are dead on arrival for `Z/48`.  Any further
A(16) search for this lane should include the cheap filter
`3 | #J(F_p)` for several good primes before exact torsion or cubic-contact
work.

