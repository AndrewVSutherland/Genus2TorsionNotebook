# The (2,2,2,12) moduli surface: factorization, configuration structure, and an exhaustive search

Date: 2026-07-19 (Opus session, from Jen's (A,B,C) equations)

Jen's surface, in `P^2_{A,B,C}`:

```text
S:  y^2 = F1 = A^4B^2 - 2A^2B^4 + B^6 - A^4C^2 + 3A^2B^2C^2 - 2B^4C^2 - A^2C^4 + B^2C^4
    z^2 = F2 = A^6 - 2A^4B^2 + A^2B^4 + 2A^4C^2 - 3A^2B^2C^2 + B^4C^2 + 2A^2C^4 - 2B^2C^4 + C^6
```

with `(A,B,C) -> (s,m,n) -> ` the `M(2,2,2,6)` quintic, giving `[2,2,2,12]`.

## 1. Both sextics factor, and they share a conic

Verified symbolically (`code/opus_abc_verify.m`):

```text
F1 = (A-B)(A+B)(B-C)(B+C) * W          [four lines times a conic]
F2 = W * Q
W  = A^2 - B^2 + C^2                   [the SHARED conic]
Q  = A^4 - A^2B^2 + A^2C^2 - B^2C^2 + C^4     [irreducible]
```

Consequences:

- `S` is **not** a generic fibre product of two K3s.  Because the two
  branch sextics share the component `W`, `S` is the `(Z/2)^2`-cover of
  `P^2` with the three branch divisors `W` (deg 2), the four lines
  (deg 4), and `Q` (deg 4); the third quotient is `w^2 = P1*Q` (deg 8),
  `P1 = (A^2-B^2)(B^2-C^2)`.
- The two conditions are equivalent to the single statement:
  **`W`, `P1`, `Q` all lie in the same square class.**
- Kodaira dimension: for a `(Z/2)^2`-cover with reduced branch data of
  degrees `2+4+4 = 10`,
  `K_S = pi^*(K_{P^2} + (1/2)(D1+D2+D3)) = pi^*(O(-3) + O(5)) = pi^*O(2)`,
  which is ample: **general type (`kappa = 2`)**, hence Bombieri--Lang
  predicts the rational points are NOT Zariski dense.  Caveat: the branch
  configuration is very degenerate (4 lines + conic + quartic with many
  intersections), so a resolution could drop `kappa`; the computation
  above is the naive one.
- Every "obvious" rational curve on `S` is a degenerate locus: on
  `A = +-B`, `B = +-C`, `A = 0`, `C = 0` one has `F1 = 0`, and on `W = 0`
  both vanish.  These are exactly factors of Jen's discriminant, so the
  trivial points are completely accounted for and the interesting points
  are sporadic.

## 2. The configuration structure (the useful discovery)

The 12 representations of each known curve organize as: **three** values
`{x1,x2,x3}` for the `A,C` slots and **two** values `{b1,b2}` for the `B`
slot, with all `3 x 2` ordered combinations `(x_i, b_k, x_j)` on `S`
(`S_3 x Z/2`, order 12 — this matches the independent "12 CURVE1 + 12
CURVE2" fingerprint census in `data/claude_sib_t5surf_hitclassify.txt`).

```text
curve #2:  x = {120, 143, 266},   b = {218, 241}     N = 105605
curve #1:  x = {143, 408, 1015},  b = {437, 1013}    N = 1217138
```

and in both cases, exactly:

```text
    x1^2 + x2^2 + x3^2  =  b1^2 + b2^2  =  N .
```

Writing `alpha_l = x_l^2-b1^2`, `beta_l = x_l^2-b2^2`, the surface
conditions at `(x_i,b_k,x_j)` become `W = -beta_l`, `F1 = alpha_i alpha_j beta_l`,
so the whole system is equivalent to:

```text
(C1)  x1^2+x2^2+x3^2 = b1^2+b2^2 = N
(C2)  gamma_l := (x_l^2-b1^2)(x_l^2-b2^2) = square,   l = 1,2,3
(C3)  P := prod_l (x_l^2 - b1^2) = square   (equivalently R = prod beta_l)
```

**(C2) says each `x_l` is a rational point of the SINGLE elliptic curve**

```text
E_{b1,b2}:  w^2 = (x^2-b1^2)(x^2-b2^2)
        <=> Y^2 = X(X-b1^2)(X-b2^2),     X = x^2,  Y = x*w
```

i.e. a point of `E` whose `X`-coordinate is a **square**.  So:

> A `(2,2,2,12)` curve in this family is exactly a pair `(b1,b2)` together
> with three points of `E_{b1,b2}` with square `X`-coordinate whose
> `X`-coordinates sum to `N = b1^2+b2^2`, subject to `(C3)`.

Verified exactly for both curves: the three `x` values are precisely the
nontrivial integer points of `E_{b1,b2}` (checked to `x < 4000`), and
`P`, `R` are squares.  In every case tested, `F2 = square` came along for
free once `(C1)-(C3)` held, so `(C3)` appears to be the only extra
condition; that implication was not proved symbolically.

### Boundedness (why a box search is the right tool)

Since `X_l = x_l^2 > 0` and `sum_l X_l = N`, every `x_l < sqrt(N)`, and
`b_k < sqrt(N)`.  Hence **all five coordinates of a configuration are
within a factor `sqrt(3)` of one another** — the rational points cannot
have long thin tails (large `C`, small `A,B`).  A box search of side `H`
is therefore complete for *every configuration whose largest coordinate
is at most `H`*, and no unbounded Mordell--Weil enumeration is needed.
(The elliptic-fibration search over the `(A,B)`-plane,
`code/opus_abc_fibration.m`, was built and then abandoned for exactly this
reason: it is not needed, and `MordellWeilGroup` on those fibres is slow.)

## 3. Searches run

All code in `code/opus_abc_*`, logs/data in `results/` and the session
scratchpad.

**(a) Direct box search** — `code/opus_abc_search.c`.  Uses the
factorization, an 8-small-prime residue sieve indexed by `C mod p`
(rebuilt per `(A,B)`), 3 medium primes, then exact `__int128` square
tests; applies Jen's full 16-factor discriminant condition; exploits the
genuine automorphism `(A,B,C) -> (C,B,A)` of `F1,F2` to halve the domain.
Throughput ~7 ns/triple (height 600 in 1.6 s).

```text
COMPLETE to height 10000  (8000x the volume of the height-500 scan)
  exact 128-bit tests : 61,274,954
  degenerate hits     : 28,768
  non-degenerate hits : 339
  distinct primitive classes : 12   (6 for each known curve)
  NEW curves          : NONE
```

**(b) Configuration search** — `code/opus_config_search.c`.  Scans
`b1<b2`, and for each pair finds all `x < sqrt(N)` on `E_{b1,b2}`, then
looks for triples with `sum x^2 = N`; post-filtered by `(C3)`.

```text
COMPLETE for all configurations with b1 < b2 <= 6000
  pairs scanned                    : 17,997,000
  triples passing (C1)+(C2)        : 37
  distinct primitive configurations:  4
  GENUINE (also passing (C3))      :  2   = exactly the two known curves
```

The two rejected primitive configurations are

```text
b=(1000,1261) x=(275, 936,1280)   P not square -> F1,F2 not squares
b=(1254,1449) x=(608,1078,1463)   P not square -> F1,F2 not squares
```

so `(C3)` is a real and necessary cut, not a formality.

Census of `#{x < sqrt(N) : x on E_{b1,b2}}` over all 18M pairs:

```text
 0 pts : 17,063,073      4 pts :  15,613      8 pts :  505
 1 pt  :    697,445      5 pts :   5,935      9 pts :  258
 2 pts :    166,450      6 pts :   2,589     10 pts :  133
 3 pts :     43,778      7 pts :   1,093    >=11    :  144
```

So ~95% of pairs give no usable point at all, and although `43,778`
pairs do carry three or more points, only `37` triples anywhere in the
range have `X`-coordinates summing exactly to `N` — and only `2` of
those survive `(C3)`.  The bottleneck is the exact sum, not the supply
of points.

**(c) Pipeline validation** — `code/opus_abc_classify.m`.  The seven
primitive `(A,B,C)` classes of height `<= 600` map through
`(s,m,n)` to exactly **two** curves by `G2Invariants`; both have exact
torsion `[2,2,2,12]` and two-prime root-power geometric-simplicity
certificates (`p=71,103` for curve #2; `p=37,73` for curve #1).

## 4. Are there infinitely many interesting points?

Evidence points to **very few, plausibly finitely many**:

- `K_S = pi^*O(2)` (naively ample) puts `S` in general type, so
  Bombieri--Lang predicts non-density.
- The configuration census is brutal: over the scanned range almost every
  pair `(b1,b2)` yields **0 or 1** point of `E_{b1,b2}` with square `X`
  below `sqrt(N)`; three such points with an exact sum `N` is a triple
  coincidence.
- Heuristic count: there are `~B^2/2` pairs with `b <= B`; for each, `N`
  is of size `~B^2` and the chance a given triple of `X`-values sums
  exactly to `N` is `~1/N ~ 1/B^2`.  With `O(1)` points per curve on
  average, the expected number of solutions up to height `B` is
  **`O(1)`**, growing only as fast as the average point count itself
  (very slowly, with rank).  This matches the data: two curves up to
  10000, and it predicts the next one — if any — needs a substantially
  larger search, not a slightly larger one.

## 5. Recommended next steps

1. Prove `(C1)+(C2)+(C3) => F2 = square` symbolically; that would reduce
   the whole family to a clean statement about `E_{b1,b2}` and make the
   census above rigorous rather than empirical.
2. Prove that every non-degenerate rational point of `S` extends to a
   configuration (observed for both curves and consistent with the 12-fold
   representation census, but not proved).  With (1) it would turn the box
   search into a genuine classification theorem in the scanned range.
3. If more curves are wanted, run the configuration search — not the box
   search — over `b <= 10^5` on a cluster; it is the cheaper of the two
   per unit of coverage, and shard by `b2` (work grows as `b2^3`).
4. Target high-rank `E_{b1,b2}`: the whole game is pairs whose elliptic
   curve carries three square-`X` points below `sqrt(N)`, so prescreening
   `(b1,b2)` by rank (or by the count of small square-`X` points) is far
   more efficient than uniform scanning.
