# Route B2' for generic `[2,22]`: the first search of the `[1,1,2,2]` ∩ CF-order-11 locus

Date: 2026-07-25.  Lane 3 of the overnight campaign (aws-spot-5), completed and
exact-tested on claude-box after the campaign session died.

## The construction

Both known geometrically simple `[2,22]` curves have factor type `[1,1,2,2]`
(2-rank 2) and carry their order-11 class at infinity.  The BLP2009 ansatz
`y^2 = R(x)^2 - 4c^2 S(x)^2` makes `f` a product of two cubics; the
`[1,1,2,2]` shape is exactly "each cubic acquires a rational root", and that
condition is **linear** — given `(c,d,r1,r2)` it solves for `(a,b)`.  So the
2-structure is free, the chart becomes the monic

```text
y^2 = x (x - 1) (x^2 + p1 x + p2) (x^2 + p3 x + p4),
```

and the only remaining condition is `ord(D_inf) = 11`, which the polynomial
continued fraction decides exactly in microseconds with no Jacobian arithmetic.
This is the first time that locus has been searched at all.

## What was run

- C sieve over `(p1,p2,p3,p4)` at height 60 (4,407 rationals per coordinate,
  5-prime divisibility gate, `recbound=500`), 192-way on aws-spot-5, ~6.7 h.
  `code/claude_ov_b2p_scan.c`, `results/claude_ov_b2p_sweep_H60.log`.
- Exact stage on all **288** distinct candidate tuples: exact `CFOrder(D_inf)`,
  factor type, exact `TorsionSubgroup`, the real-subfield-disc RM census over
  the good primes `p <= 200`, a two-prime strict `End = Z` certificate attempt,
  and absolute (Igusa) invariants as an isomorphism fingerprint.
  `code/claude_ov_l3_exact.m`, `results/claude_ov_l3_exact.log`.
  (Three bugs were fixed to make the exact stage run: the candidate parser
  accepted only integers, `TorsionSubgroup` needs an integral model — the sweep
  emits rational `p_i`, so scale `x -> X/m`, `y -> Y/m^3` — and the Frobenius
  loop had to move to the same integral model.)

## Result

```text
288 candidates ->  270 exact torsion [2,2]        (false positives of the 5-prime gate)
                     1 exact torsion [2,2,2]
                    17 exact torsion [2,22]
```

All 17 have `CFOrder(D_inf) = 11` and factor type `[1,1,2,2]`, as designed.
But grouping by absolute invariants, the 17 tuples are only **TWO distinct
curves**, in 12 and 5 coordinate presentations respectively — and they are
**exactly the two already-known witnesses**:

- curve 1 (12 presentations, incl. `(1,6,-6,32)` and `(4,27,-3,8)`) = `19044.h.2`;
- curve 2 (5 presentations, incl. `(-34/15,32/45,31/15,-2/9)`) = the corrected
  BLP row C4.

Both positive controls therefore reproduce, which validates the funnel end to
end.  **No new `[2,22]` curve appears at height 60.**

## The structural point

Every one of the 17 tuples returns

```text
ndisc = 1,  discs = [5],  RMflag = true
```

i.e. the squarefree core of the real-Weil-subfield discriminant is *constantly*
5 across every good prime `p <= 200`.  Every `[2,22]` point on this locus at
accessible height is `RM(sqrt 5)`, and no tuple anywhere in the 288 produced a
scattering signature with `[2,22]`.

This is worth stating carefully, because it is stronger than a null result and
weaker than a theorem.  Addendum 2 of the plan note established that the
*ambient* BLP stratum is NOT inside Humbert-5 (15 of the 18 BLP curves have an
`End = Z` signature).  What this search suggests is that the **2-rank-2
sublocus** of it may be Humbert-forced even though the ambient is not — which
would explain, rather than merely record, why `[2,22]` is RM-only.

**Next computation (the one that would settle it):** write the Humbert-5
equation on the `(p1,p2,p3,p4)` chart and test whether the CF-order-11 locus is
contained in it, symbolically rather than by sampling.  If it is, `[2,22]`
generic is impossible on this chart and the target moves to a genuinely
different 11-source; if it is not, the generic points are simply of larger
height and the sweep should go to height 200+ (the C sieve is already written
and shards 192 ways).

## Caveats

Height 60 in all four coordinates is not a large box, and the sieve used only a
5-prime gate before the exact stage (270 of 288 survivors were spurious), so
the false-positive budget was loose rather than the search being deep.  This is
a calibrated negative on a locus that had never been searched, not evidence of
impossibility.
