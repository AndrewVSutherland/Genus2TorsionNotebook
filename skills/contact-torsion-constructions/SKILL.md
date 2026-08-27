---
name: contact-torsion-constructions
description: Build a marked rational n-torsion class by contact (f = h^2 - c*(x-r)^n) and an INDEPENDENT rational 3-torsion class (f = h3^2 + kappa*q3^3). WHEN constructing a curve with a prescribed torsion class — any contactN route (quintic-contact-5, contact-6, contact-7, contact-9), or adding an independent 3-torsion class to reach [3,6], [6,6], [2,24], [3,12].
---

# Contact torsion constructions

## When to use this

Load this when you need to **force** a rational torsion class of a chosen
order onto a genus-2 curve `C: y^2 = f(x)`, rather than search for one. Two
independent moves live here:

1. **Marked-class contact** `f = h^2 - c*(x-r)^n` — gives one marked rational
   class of order dividing `n` (generically exactly `n`), coming from the
   function `h - y`. This is the primitive behind quintic-contact-5,
   contact-6, contact-7, contact-9.
2. **Independent cubic-contact 3-torsion** `f = h3^2 + kappa*q3^3` — bolts a
   rational 3-torsion class onto an existing curve, independent of any marked
   class. Use it to compose coprime orders: order-8 + order-3 = `Z/24`,
   marked order-6 + independent order-3 = `[3,6]`/`[6,6]`, etc.

Coprime composition (order `a` class + order `b` class, `gcd(a,b)=1`, gives
order `ab` for free) is the whole reason move 2 exists. See
`target-playbook` and the hub `g2-torsion-lab`.

---

## 1. Marked-class contact: `f = h^2 - c*(x-r)^n`

### The pattern

Take `deg h = ceil(n/2)`. Write

```text
f = h^2 - c*(x-r)^n
```

Then `div(h - y) = nP - n*infinity` where `P = (r, h(r))`, so the class

```text
D = P - infinity,   P = (r, h(r))
```

has order dividing `n`, generically **exactly** `n`, provided `f` is
squarefree of the correct degree and `h(r) != 0` (a nonzero `y`-coordinate at
`P`, otherwise `P` is a Weierstrass point of order 2). The contact is at the
single point `x = r`, where `f - h^2` vanishes to order `n`.

Different pole structures at infinity realize different `n`; the higher-order
contact families (`contact7`, `contact9`) use `x*y - h` and `x^2*y - h` to
raise the pole order — see `notes/contact7_family.md`,
`notes/contact9_family.md`. The two instances below are the load-bearing ones
for this skill; copy them verbatim.

### Worked instance A — quintic-contact-5 (order 5)

From `notes/how_we_found_2220_examples.md` (and the elementary form in
`notes/m10_quintic_contact5.md`). Here `n = 5`, `deg h = ceil(5/2) = 3` but
the used `h` is degree 2 with the contact point at `x = 0`:

```text
h(x) = 1 + t*x + ((t^2 - 1)/2)*x^2,
f(x) = h(x)^2 - ((t + 1)^4/4)*x^5.
```

The contact identity gives a rational 5-torsion class. `P = (0, 1)` here
(the marked point is `x = 0`, `h(0) = 1`); `div(y - h) = 5(P) - 5(infinity)`
when `f` is squarefree of degree 5.

The elementary two-parameter version (from `notes/m10_quintic_contact5.md`),
with the marked point also at `x = 0`:

```text
h(x) = 1 + a*x + b*x^2,
f(x) = h(x)^2 - (1 + a + b)^2*x^5.
```

On `C: y^2 = f(x)` the point `P = (0,1)` satisfies
`div(y - h(x)) = 5(P) - 5(infinity)`, so `[P - infinity]` is rational
5-torsion when `f` is squarefree of degree 5. The normalization forces `x = 1`
to be a root of `f`.

In the one-parameter subfamily of `how_we_found_2220_examples.md` there is an
explicit extra divisor class

```text
H = [x^2 + 2*x/(t+1), (t+2)*x + 1]   with   2H = [x-1, 0],
```

so a smooth specialization has a rational point of order `20` (the marked
5-torsion doubled into a Weierstrass-2 class). That note then factors the
residual quartic `f/(x-1)` to build in extra rational 2-torsion, landing the
geometrically simple `[2,2,20]` example.

### Worked instance B — contact-6 (order 6)

From `notes/contact6_m36.md`. Here `n = 6`, `deg h6 = 3`, monic, contact at
`x = 1`:

```text
h6 = 1 + a*x + b*x^2 + x^3,
f  = h6^2 - (x - 1)^6.
```

The resulting `f`-coefficient expansion (copy VERBATIM from the note):

```text
f = (2*b + 6)*x^5
  + (2*a + b^2 - 15)*x^4
  + (2*a*b + 22)*x^3
  + (a^2 + 2*b - 15)*x^2
  + (2*a + 6)*x.
```

Marked class: for smooth `f`, with side conditions

```text
b != -3          (else lead coeff 2*b+6 = 0, model drops to degree 4),
h6(1) = a+b+2 != 0   (else P is a Weierstrass point, wrong order),
```

the function `h6 - y` has divisor

```text
6P - 6*infinity,    P = (1, a+b+2),
```

so `D = P - infinity` is a marked class of order dividing `6`, generically
exact order `6`. Note `f(0) = 0`, so this chart has **one** rational
Weierstrass point built in (`T0 = (0,0) - infinity`), but not full rational
2-torsion — it is deliberately *not* the split `[2,2,2,6]` family.

Useful factorization of the contact-6 `f` (from the same note), which exposes
the 2-rank structure:

```text
f = x*((b+3)*x^2 + (a-3)*x + 2)*(2*x^2 + (b-3)*x + (a+3)).
```

When the two quadratics are irreducible over `Q`, the curve has factor type
`[1,2,2]` and rational 2-rank 2 — the productive locus for simple `[6,6]`
(see `two-rank-and-factor-types`, `halving-and-doubling`).

---

## 2. Independent cubic-contact 3-torsion: `f = h3^2 + kappa*q3^3`

This is the **VALIDATED** identity from `notes/agent_a2_24_contact.md`
(validated on both known simple `Z/24` curves) and implemented in
`code/agent_a2_24_contact_extract.m`.

### The identity

A rational order-3 class `D = (u, v)` (Mumford; `u` monic deg 2, `v` deg <= 1)
on `J(y^2 = f)`, `deg f = 6`, gives the cubic-contact identity

```text
f  =  h3^2  +  kappa * q3^3
```

with

```text
q3 := u          (the "contact conic", monic degree 2),
h3               (the "contact cubic", degree 3),
kappa            a constant,
```

and the leading-coefficient relation

```text
lead(f) = lead(h3)^2 + kappa.
```

The 3-torsion class is `[q3, h3 mod q3]` in Mumford coordinates.

Geometric origin (from the note): `div(y - h3) = (zeros of f - h3^2) - 3*inf+
- 3*inf-`. For order 3 we need `f - h3^2 = kappa*q3^3`, i.e. `3*(P1+P2) -
3*(inf+ + inf-)` where `P1, P2` lie over the roots of `q3` with `y = h3`.

**Non-degeneracy side conditions** (from `contact6_m36.md`, in that note's
convention): `m*v*(U^2 - 4*v^2) != 0` and `gcd(q, f) = 1` (here `v` is the
constant of `q = x^2 + U*x + v^2`; see below).

**Two flavours of the 3-torsion divisor** (from `agent_a2_24_contact.md`):

- `q3` irreducible (e.g. curve A): the two contact points are conjugate over
  a real quadratic field; the class is still rational.
- `q3` split (e.g. curve B, `q3 = x*(x - 2/3)`): the 3-torsion is a sum of two
  rational points `P1 + P2`.

This split/non-split dichotomy is a natural way to stratify a 3-torsion search.

### The equivalent contact-6 convention (same identity, opposite sign)

`notes/contact6_m36.md` states the same identity as

```text
h3^2 - f = m^2*q^3,   with   q = x^2 + U*x + v^2,   h3 = m*x^3 + N*x^2 + R*x + S.
```

This is the **same identity up to sign**: comparing `f = h3^2 - m^2*q^3` with
`f = h3^2 + kappa*q3^3` gives

```text
kappa = -m^2.
```

(Correspondingly the extraction script header of
`code/agent_a2_24_contact_extract.m` writes `f = h3^2 - kappa*u^3` internally,
then reports `f = h3^2 + kappa*q3^3` after a sign flip; both live in that one
file. Watch the sign convention when you copy formulas between the two notes.)

The contact-6 note gives explicit formulas for the cubic coefficients. Put
`L = 1/m`. If `c_i` are the coefficients of `f = c5*x^5 + c4*x^4 + c3*x^3 +
c2*x^2 + c1*x`, define

```text
B     = c5*L^2 + 3*U,
Delta = 4*c4*L^2 + 12*(U^2+v^2) - B^2.
```

Then

```text
N = B/(2L),    R = Delta/(8L),    S = v^3/L,
```

and the remaining `[3,6]` cover (the three cubic-contact equations that must
vanish) is

```text
B*Delta + 16*v^3 - 8*c3*L^2 - 8*U^3 - 48*U*v^2 = 0,
Delta^2 + 64*B*v^3 - 64*c2*L^2 - 192*(U^2*v^2+v^4) = 0,
Delta*v^3 - 4*c1*L^2 - 12*U*v^4 = 0.
```

Together with exact order of the marked class `D` and independence of this
3-torsion class from `2D`, this gives the `[3,6]` condition on the contact-6
chart.

### When to use it

Use `f = h3^2 + kappa*q3^3` when you want to **add** a rational 3-torsion
class that is independent of a marked contact class, to compose coprime
orders:

- order-8 marked class (A(8) chart) + independent order 3 = `Z/24`
  (`code/agent_a2_24_composite8x3.m`), and with 2-rank >= 2, `[2,24]`;
- marked order-6 (contact-6) + independent order 3 = `[3,6]`, then halve to
  `[6,6]` or add a rational Weierstrass point (`contact6_m36.md`);
- order-12 (A(12)) + independent order 3 = `[3,12]`.

**Warning (thin-set condition).** Rational 3-torsion cut this way is a
*covering* condition, not a hypersurface: `notes/agent_a2_24_d0_derivation.md`
proves that on the tractable `d = 0` slice of the A(8) chart the genuine
`kappa != 0` locus collapses (the Rabinowitsch ideal `<Eq4, Eq3,
w*kappa - 1>` over `Q(r,t)` is the unit ideal), and the genuine cover `X` is a
**genus-1, rank-0** curve carrying essentially only one genuine point (curve
B). So do not expect a positive-dimensional 3-torsion family from a naive
slice; expect a thin cover. This is *why* `[2,24]` is hard.

---

## 3. Extraction: recover `(q3, h3, kappa)` from a known 3-torsion class

Implemented in `code/agent_a2_24_contact_extract.m`, function
`ExtractContact(f, u, v)`. Use it to *recover* the contact data from any curve
already known to have a rational order-3 Mumford class `(u, v)` with `u` monic
degree 2 — e.g. a class returned by `TorsionSubgroup(J)`.

Algorithm (verbatim logic from the script):

1. Since `(root(u), v)` lies on the curve, `u | (f - v^2)`. Compute the
   quotient `g4 = (f - v^2) div u` (degree 4). If the remainder is nonzero,
   `(u, v)` is not a valid 3-torsion contact class — fail.
2. Write the contact cubic as

   ```text
   h3 = v + u*(s*x + w)
   ```

   with unknowns `s, w`, and impose `f - h3^2 = kappa*u^3`. Dividing by `u`
   once gives the **degree-4 coefficient identity**

   ```text
   g4 - 2*v*(s*x + w) - u*(s*x + w)^2 - kappa*u^2 = 0,
   ```

   i.e. 5 coefficient equations (degrees 0..4) in the 3 unknowns `(s, w,
   kappa)`. It is quadratic in `(s, w)`; solve exactly.
3. In the script this is done by building `R<S,W,K> := PolynomialRing(Q, 3)`,
   forming the polynomial identity in a `PolynomialRing(R)`, taking
   `Coefficient(expr, i)` for `i in [0..4]`, and calling `Variety` on the
   ideal. Take `V[1]`; set `h3 = v + u*(sv*x + wv)`, `q3 = u`, `kappa = kv`.
4. Verify: `f - (h3^2 + kappa*q3^3)` must be `0`, and `lead(h3)^2 + kappa`
   must equal `lead(f)`.

### Validated data (the two known simple `Z/24` curves)

From `notes/agent_a2_24_contact.md` and reproduced by the extraction script
(`curves := [ <5, -5/2, -9/2, ...>, <1/3, -1/9, -1, ...> ]`):

- **Curve A** (`r=5, p=-5/2, t=-9/2`; `chi_17` simple): `q3 = x^2 - 435/73*x +
  2529/292` (disc `= 4608/5329 = (48/73)^2 * 2`, IRREDUCIBLE over `Q`, splits
  only over `Q(sqrt2)`); `h3 = 1205280*x^3 - 10769040*x^2 + 31705560*x -
  30713580`; `kappa = -1451998172160`; `lead(h3)^2 + kappa = 701706240 =
  lead(f)`.
- **Curve B** (`r=1/3, p=-1/9, t=-1`; `chi_13` simple; on the A(8) `d=0`
  slice, `Q = x^2`): `q3 = x^2 - 2/3*x = x*(x - 2/3)` (SPLITS RATIONALLY);
  `h3 = 13122*x^3 + 26244*x^2 - 34992*x + 15552`; `kappa = 1377495072`; the
  order-3 divisor is `P1 + P2` with `P1 = (0, 64/243)`, `P2 = (2/3, 32/243)`,
  both rational.

---

## Pitfalls

- **Sign convention on `kappa`.** `agent_a2_24_contact.md` writes
  `f = h3^2 + kappa*q3^3`; `contact6_m36.md` writes `h3^2 - f = m^2*q^3`
  (i.e. `f = h3^2 - m^2*q^3`). They agree with `kappa = -m^2`. The extraction
  script internally uses `f = h3^2 - kappa*u^3` and then flips. Always check
  which sign a formula is in before substituting, or your identity will be off
  by a sign on the cubic term.
- **`deg h` for the marked contact.** It is `ceil(n/2)` in general, but the
  *worked* quintic-contact-5 uses a degree-2 `h` with the contact at `x = 0`
  (not `ceil(5/2) = 3`); contact-6 uses monic degree-3 `h6` with contact at
  `x = 1`. Copy the instance, don't re-derive the degree from `n` alone.
- **Forgetting the side conditions.** Dropping `b != -3` (contact-6) silently
  drops `f` to degree 4; dropping `h6(1) != 0` makes the marked point a
  2-torsion Weierstrass point of the wrong order. In searches these show up as
  "boundary" hits (see `contact6_m36.md`: at height 5 all six cubic-contact
  solutions had `b = -3`).
- **Expecting a family from rational 3-torsion.** It is a thin/covering
  condition. Do NOT model it as "solve one equation for one parameter"; the
  `d=0` slice's genuine cover is genus-1 rank-0
  (`agent_a2_24_d0_derivation.md`). Budget for a covering-space search, or
  compose with a chart that already carries the coprime factor.
- **`q3` split vs irreducible.** Both give rational 3-torsion. If you only
  accept split `q3` (two rational points) you miss curve A; if you only accept
  irreducible you miss curve B. Handle both, or state which you target.
- **No invented Magma.** The extraction idiom (`PolynomialRing(Q,3)`, build
  the identity, `Coefficient(expr,i)`, `Variety(ideal<...>)`) is exactly the
  one in `code/agent_a2_24_contact_extract.m`. Reuse it; do not invent a
  "solve for contact" intrinsic.

## See also

- `g2-torsion-lab` — the hub; the five reusable moves and the search funnel.
- `named-charts-reference` — the A(8), A(12), M(2,12), M_1(8,4) charts that
  carry a marked high-order class for free (compose them with the independent
  3-torsion here).
- `halving-and-doubling` — double the order of the marked class (`n -> 2n`);
  the halving wall.
- `two-rank-and-factor-types` — pick the factor type of `f` for extra rational
  2-torsion (the `[1,2,2]` locus of contact-6).
- `target-playbook` — decompose a target into coprime factors and pick the
  route.
