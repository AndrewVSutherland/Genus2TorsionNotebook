---
name: simplicity-certificates
description: Certify a genus-2 Jacobian J is GEOMETRICALLY SIMPLE via the good-prime Frobenius L-polynomial irreducibility + n-th-power-transform ("D4"/Leprevost/root-power) criterion, and the independent Sage/Lombardo endomorphism test. WHEN — before claiming ANY result: a non-simple curve (split, isogeny-decomposable, CM, bielliptic) is worthless for this project, so simplicity must be certified, never assumed. Trigger words — simple, geometrically simple, End(J), Frobenius, L-polynomial, chi(T), power transform, D4, Lombardo, split, isogeny factor, quadratic factorization.
---

# Certifying geometric simplicity of a genus-2 Jacobian

## When to use this

Load this skill the moment you have a candidate curve `C: y^2 = f(x)` (deg
`f in {5,6}`) whose torsion you want to publish. **Geometric simplicity is
the whole game** (`g2-torsion-lab`): torsion of a *simple* abelian surface is
the object of interest; a `[6,6]` or `[2,2,20]` curve whose Jacobian splits
into (or is isogenous to a product of) elliptic curves is worthless here. You
must certify simplicity *before* reporting a hit, and record the witness. This
skill gives:

1. the **primary certificate** — a single good-prime Frobenius test that is
   fast, self-contained, and provable (the `SimplicityCertificate` function in
   `code/agent_a2_24_composite8x3.m`);
2. the **independent historical confirmation** — the Sage/Lombardo
   endomorphism test (`geometric_endomorphism_*`);
3. the **anti-signature** — how a *non-simple* Jacobian betrays itself, so you
   can reject split curves early.

"Simple" here always means **geometrically simple**: `End(J_Qbar) ⊗ Q` is a
field (no nontrivial idempotent, hence no isogeny splitting over `Qbar`).

## The primary certificate (good-prime Frobenius, "D4"/power-transform)

Source of truth: the `SimplicityCertificate` function in
`code/agent_a2_24_composite8x3.m` (lines 75–97). Copy its logic exactly; do
not paraphrase the arithmetic.

### Step 1 — pick a good prime

Given the **integer model** `fInt` (clear denominators first — see
`validate-and-record-a-hit`; simplicity is tested on the integral model, never
the rational one), iterate `pp` over `[3,5,7,11,13,17,19,23,29,31,37,41,43,47]`.
A prime `pp` is admissible iff **all** hold:

- `pp` does not divide the leading coefficient: `(Z!LeadingCoefficient(fInt)) mod pp ne 0`;
- `pp` does not divide the numerator of the discriminant: `(Z!Numerator(dsc)) mod pp ne 0`, where `dsc := Discriminant(fInt)`;
- `f mod pp` has degree `>= 5` and is squarefree: `Degree(fp) ge 5 and IsSquarefree(fp)`.

These conditions guarantee good reduction of the genus-2 curve at `pp`, so
Frobenius acts on a genuine abelian surface `J mod pp`.

### Step 2 — build the Frobenius characteristic polynomial chi(T)

Count points over `F_p` and `F_{p^2}` with the standard hyperelliptic count
(the `CountCurve` function, lines 67–73): for each `xx in F_q`, add 1 if
`f(xx)=0`, add 2 if `f(xx)` is a nonzero square; then add 2 more iff the
**leading coefficient is a square** (this is the two-points-at-infinity
contribution for an even/sextic model — see `two-rank-and-factor-types`).

Then, with `a1`, `a2` the first two Frobenius traces:

```magma
a1  := pp + 1 - CountCurve(fp);                              // = p+1-#C(F_p)
a2  := (CountCurve(fp2) - pp^2 - 1 + a1^2) div 2;            // from #C(F_{p^2})
chi := T^4 - a1*T^3 + a2*T^2 - a1*pp*T + pp^2;              // Frobenius char. poly
```

`chi(T) = T^4 - a1 T^3 + a2 T^2 - a1 p T + p^2` is the degree-4 characteristic
polynomial of Frobenius on `J mod pp` (functional-equation symmetric: the `T`
and `T^3` coefficients are tied by the factor `p`, and the constant term is
`p^2`).

### Step 3 — the certificate test

```magma
if not IsIrreducible(chi) then continue; end if;   // reducible => not a witness at pp
K := NumberField(chi); pi := K.1; drop := false;
for nn in [2..12] do
    if Degree(MinimalPolynomial(pi^nn)) lt 4 then drop := true; break; end if;
end for;
if not drop then return true, pp, chi; end if;      // CERTIFIED simple, witness (pp, chi)
```

**Certificate.** If, at an admissible prime `pp`,

1. `chi` is **irreducible over Q**, AND
2. for a root `pi` of `chi`, `MinimalPolynomial(pi^n)` has **degree 4 (no
   degree drop)** for every `n in [2..12]`,

then `J` is **geometrically simple**, and `(pp, chi)` is the witness to record.

### Why one good prime suffices

If `J` split geometrically (isogenous over `Qbar` to a product of elliptic
curves, or had extra endomorphisms / CM), then at *every* good prime the
Frobenius eigenvalues would generate a proper subfield after passing to some
power `n` (the abelian surface becomes isogenous to a product over the field
where the split is defined, and a suitable power of Frobenius lands in the
smaller endomorphism-eigenvalue field). Concretely: a **degree drop** of
`MinimalPolynomial(pi^n)` below 4 means `pi^n` lies in a proper subfield of
`Q(pi)`, exposing exactly such a subfield / CM-by-power structure. Conversely,
if `chi` is irreducible and *no* power `pi^n` (`n <= 12`) drops degree, that
single reduction is a simple abelian surface with no power giving a subfield —
which rules out geometric splitting and CM. So **a single witnessing good
prime certifies geometric simplicity of `J/Q`.**

The `[2..12]` power range is the standard root-power / Leprevost "D4" window
used throughout the lab; it covers the small cyclotomic obstructions that a
CM/split surface would exhibit.

### What a degree drop means (retry, do not reject)

A degree drop at some `n` means `pi^n` lies in a proper subfield, so `pp` is
**not a witness at that prime** — but this does *not* prove `J` non-simple. Try
the next prime. Only after exhausting the prime list (function returns
`false, 0, 0`) is the certificate inconclusive; then widen the prime list or
fall back to the Lombardo test below. Contrast this with the reducible-`chi`
case, which is also just "not a witness here, try another prime" — reducibility
at one good prime does NOT by itself prove non-simplicity either (see Pitfalls
for the case that *does*).

## Worked examples (grounded in the notes)

**`[2,2,20]` at p=71** (`notes/how_we_found_2220_examples.md`, lines 91–97).
The documented geometrically-simple `[2,2,20]` curve has good-prime
certificate at `p=71` with Frobenius polynomial

```text
chi = X^4 + 2*X^3 + 14*X^2 + 142*X + 5041
```

(`5041 = 71^2`, matching the `p^2` constant term), whose **12th-power
transform is irreducible** — i.e. no degree drop at `n=12` (hence at every
divisor of 12). This is the gold-standard recorded witness.

**`[6,6]` at p=23** (`notes/contact6_m36.md`, line 576). The
geometrically-simple `[6,6]` curve found on the `[1,2,2]` core
(`a=133/39, b=-7/13`, `y^2 = 11389248*x^5 - 18252000*x^4 + 42399396*x^3 -
10288044*x^2 + 29659500*x`) has

```text
simple certificate: p=23, L_p(T) = 529*T^4 - 26*T^2 + 1 irreducible
```

(`529 = 23^2`; here `L_p(T)` is the reciprocal `T^4 chi(1/T)` normalization,
with leading coefficient `p^2` — the irreducible degree-4 local factor is the
same witness).
This curve has factor type `[1,2,2]`, *not* the extra-root type `[1,1,1,2]`
(see `two-rank-and-factor-types`).

## Independent confirmation: the Sage/Lombardo endomorphism test

The historical geometric confirmation is Lombardo's endomorphism-algebra
computation in Sage, over `J = C.jacobian()`
(`code/m212_extra3_geom_simple_check.py`; results in
`notes/how_we_found_2220_examples.md` lines 84–89 and
`notes/contact5_order40_family.md` lines 341–342):

```python
J.geometric_endomorphism_algebra_is_field(B=100)  # True  => End^0(J_Qbar) is a field
J.geometric_endomorphism_ring_is_ZZ(B=100)        # True  => End(J_Qbar) = Z (generic, RM-free)
```

Both `True` means `End(J_Qbar) ⊗ Q` is a field with `End(J_Qbar) = Z`: `J` is
geometrically simple with no extra endomorphisms. This is the confirmation
recorded for the `[2,2,20]` example. Note the house rule (`AUTHORING.md`):
**the project runs in Magma**; the Lombardo/Sage test is the *one* historical
exception, used as an independent cross-check, not the day-to-day tool. Prefer
the Magma Frobenius certificate above for new hits and cross-check with
Lombardo when you can.

(For a *negative* example the same Sage call returns `False`; see
`notes/m212_three_torsion.md` lines 162–163, where
`geometric_endomorphism_algebra_is_field(B=100) = False` flags a non-simple
Jacobian.)

## The anti-signature: how a NON-simple Jacobian betrays itself

Learn to reject split curves fast — several long lab searches produced
`[6,6]`/`[3,12]` curves that were **all non-simple**, wasting exact-torsion
time.

- **L_p factoring into quadratics at good primes is THE signature of a
  non-simple (split / isogeny-decomposable) Jacobian.** In the non-simple
  `[2,6,6]` extra-root examples of `notes/contact6_m36.md` (lines 459–468),
  *every* checked good-prime local `L`-polynomial factored into quadratics,
  e.g.

  ```text
  p=11: L_p = (11*T^2 + 1)^2
  p=13: L_p = (13*T^2 - 2*T + 1)^2
  ```

  When `L_p` (equivalently `chi`) is a product of two quadratics at good
  primes, `J mod p` is isogenous to a product of two elliptic reductions — the
  fingerprint of a split Jacobian. **Reject it.**

- If the SAME quadratic factorization pattern persists at several good primes,
  and Sage returns `geometric_endomorphism_algebra_is_field = False`, the
  Jacobian is genuinely split. You can then *exhibit* the split with the
  `Degree2Subcovers` / `AutomorphismGroup` route
  (`code/m212_extra3_split_certificate.m`): it produces two degree-2 maps to
  elliptic curves (e.g. identified as Cremona `90c3` and `510g1`), which is a
  positive *split* certificate — the opposite of what you want.

## Pitfalls (how a cheaper model gets this wrong)

- **Never claim "simple" from irreducibility of `chi` at a single prime alone.**
  Irreducible `chi` without the `n in [2..12]` power-transform check misses
  CM/subfield surfaces whose `chi` is irreducible but some `pi^n` drops degree.
  The power-transform loop is not optional.
- **Never test simplicity on the rational (unscaled) `f`.** Run the certificate
  on `fInt` (denominators cleared), because `LeadingCoefficient`,
  `Discriminant`, and the `mod p` reductions must be integral for the
  good-prime tests to be valid. See `validate-and-record-a-hit`.
- **A single reducible `chi` does NOT prove non-simplicity** — it just means
  that prime is not a witness; try another. What *does* signal non-simplicity is
  `chi` splitting into **quadratics at many good primes** (the anti-signature
  above), ideally confirmed by Lombardo `= False` and an explicit degree-2
  subcover.
- **Do not confuse this with the covering/boundary "certificates"** in files
  like `code/agent_m18_416_R25_4_SB_v4_certificate.m` — those certify that a
  *search locus* has no non-boundary rational points (a rank-zero elliptic
  quotient argument), a completely different notion from geometric simplicity
  of a found Jacobian.
- **Bielliptic even sextics are never simple.** Searches pre-skip them (in
  `agent_a2_24_composite8x3.m`, the `&and[Coefficient(f,i) eq 0 : i in [1,3,5]]`
  test drops even sextics because they cover elliptic curves and split). Don't
  send an obviously-bielliptic curve to the certificate expecting `true`.

## See also

- `g2-torsion-lab` — the hub; the two-things-must-hold rule for any result.
- `two-rank-and-factor-types` — the point-count leading-coefficient-square
  term, the `CountCurve` idiom, and factor types (the `[1,2,2]` simple locus vs
  the `[1,1,1,2]` split extra-root locus).
- `validate-and-record-a-hit` — where the certificate sits in the full hit
  checklist, and how to record the witness `(pp, chi)`.
