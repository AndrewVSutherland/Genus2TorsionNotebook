# Lane 1 — the contact-7 `u = -1/2` family: `[2,2,14]` becomes a THEOREM

2026-07-25.  Three sessions: the overnight run (killed before it wrote a note),
a resumed session (committed `d36bb18`, `981cdc7`, `e824c15`, then killed by an
API outage), and this one.  All computation on **claude-box** (32 vCPU /
128 GB), all Magma through `code/claude_magma_slot.sh`.  No spot box was
started.

---

## 0.  One-paragraph summary

The contact-7 three-root surface `R(s,t,u) = 0` — whose rational points are
exactly the genus-2 curves in the contact-7 chart whose quintic has factor type
`[1,1,1,2]`, i.e. torsion containing `[2,2,14]` — carries an explicit **rational
curve**: the fibre `u = -1/2`, whose smooth model is the rank-1 elliptic curve
`E : y^2 = x^3+6x^2-16` of conductor 288.  This session established:

* **every** nondegenerate member has `J(Q)_tors ⊇ [2,2,14]` and rational 2-rank
  **exactly** 3 — unconditional, no per-member computation (§2);
* **57** members (`k = 1..30`, `eps = 0,1`; coefficient heights up to 1177
  decimal digits) have `J(Q)_tors` **exactly** `[2,2,14]` *and* a **strict
  two-prime `End(J_Qbar) = Z` certificate** — every single one (§3);
* a **congruence argument** upgrades this to an **infinite** subfamily, modulus
  `N = 11` (§4);
* a **residue-class covering** upgrades *that* to **thousands of infinite
  classes at once**: **7898 of the 7920** congruence classes `(k mod 3960, eps)`
  are certified outright (**99.72%**), including every genuine member with
  `|k| ≤ 40`; of the 22 left, 8 are the classes of the degenerate points, which
  reduction can never see, and the rest are prime-pool-limited (§5);
* the order-112 rung (`[2,2,2,14]`) on this family is **closed by Chabauty**,
  unconditionally, for every member (§2.3);
* within the `u = const` pencil, `u = -1/2` is the **only** fibre with
  infinitely many rational points, so a second infinite component of the
  surface would have to be a multisection (§6).

As far as this repo's literature intake goes (`notes/internet_infinite_families_sweep.md`,
`notes/claude_torsion_refs_dossier.md`, LPS2004 / BLP2009 / Flynn /
Daowsud–Schmidt), no infinite family with torsion containing `(Z/2)^2 x Z/7`
was previously recorded; the `[2,2,14]` row of `torsion_realizations` was filled
by five sporadic curves (commit `39b71da`) plus Lane 7's five more.

---

## 1.  The construction (verbatim, reproducible)

Contact-7 chart:

```text
h(x) = 1 - (7/2) x + a x^2 + b x^3 ,      f(x) = (h(x)^2 + (x-1)^7)/x^2
```

`f` is a monic quintic; `C : y^2 = f(x)`; the marked class is
`D = [(1, h(1)) - infty]`.  A rational root `r` of `f` is `r = 1 - v^2` with
`v` a rational root of

```text
Q(v) = v^5 + c4 v^4 + (2c4-1) v^3 + (c4+c0-1/2) v^2 + 2 c0 v + c0 ,
c4 = b+2 ,  c0 = 5/2-(a+b) .
```

Three rational roots `(s,t,u)` <=> `R(s,t,u) = 0`, the `S3`-symmetric
tridegree-`(3,3,3)` surface (rebuilt symbolically and re-verified against the
eleven known points in `results/claude_ov_lane1_family.log`); with
`e1,e2,e3` the elementary symmetric functions,

```text
R = 2(e2-2) e3^2 + (4 e2^2 + 4 e1 e2 - 2 e1 - 4 e2 - 1) e3 + e2(2(e1+e2)^2 + e1).
```

**The family.**  Put `u = c`, `p = s+t`, `q = st`.  `R` becomes a plane cubic
in `(p,q)` whose linear part is `c^2(2c+1) p`, so the origin is a singular
point of the fibre exactly for `c in {0, -1/2}`; `c = 0` is chart-degenerate
(`v = 0` is the marked root `x = 1`).  So `c = -1/2` is the unique
non-degenerate nodal fibre; its cubic

```text
-p^3 + 3 p q^2 + 2 p q + 2 q^3 = 0
```

is nodal at the origin, hence rational: `p = m q` gives

```text
q(m) = 2m / (m^3 - 3m - 2) ,   p(m) = 2m^2 / (m^3 - 3m - 2) ,
```

and `s,t` are rational iff `p^2 - 4q` is a square, i.e. iff

```text
w^2 = g(m) := -m^4 + 6 m^2 + 4 m = -m (m+2) (m^2 - 2m - 2) .
```

That quartic is the elliptic curve

```text
E : y^2 = x^3 + 6 x^2 - 16 ,   x = 4/m ,  y = 4w/m^2 ,
conductor 288,  disc 2^12 * 3^3,  E(Q) = <G=(4,-12)> x <T=(-2,0)> = Z x Z/2 .
```

**Member recipe.**  For `P = k G + eps T` in `E(Q)`:

```text
m = 4/x(P),  w = y(P) m^2/4,  d = (m+1)^2 (m-2),
s = (m^2+w)/d,  t = (m^2-w)/d,  u = -1/2,
G7(v) = -(v^5 - v^3 - v^2/2)/(v+1)^2,
c4 = (G7(s)-G7(t))/(s^2-t^2),  c0 = G7(s) - c4 s^2,
b = c4 - 2,   a = 9/2 - c0 - c4,
h = 1 - (7/2)x + a x^2 + b x^3,   f = (h^2+(x-1)^7)/x^2 .
```

Verified identically on `R` over `Q(m)[w]/(w^2-g(m))`
(`results/claude_ov_lane1_family.log`, line `PARAMETRISATION identically on R:
true`).

**The degenerate points.**  The construction fails exactly when
`m in {0, ±1, ±2, ∞}`, i.e. `x(P) in {∞, ±4, ±2}`, i.e. at the **eight** points

```text
D = { O , ±G , ±G+T , ±2G , T }   =   (k,eps) with k in {0,±1,±2} :
  (0,0)=O,  (0,1)=T [x=-2],  (±1,0)=±G [x=4],  (±1,1)=±G+T [x=-4],  (±2,0)=±2G [x=2].
```

Note `(±2,1)` is *not* degenerate.  Everything else in `E(Q)` gives a genuine
genus-2 curve.  (`results/claude_ov_lane1_cover.log`, first block.)

---

## 2.  Theorem A (containment) — unconditional, every member

> For every `P in E(Q) \ D` the curve `C_P` is a genus-2 curve whose Jacobian
> satisfies `J(Q)_tors ⊇ (Z/2)^3 x Z/7 = Z/2 x Z/2 x Z/14`, with rational
> 2-rank **exactly 3**.

Proof, in three pieces, none of which needs a per-member computation.

1. **2-rank ≥ 3.**  `1-s^2, 1-t^2, 1-u^2` are three distinct rational roots of
   the quintic `f`, so `f` has ≥ 4 irreducible factors over `Q`; for
   `deg f = 5` the rational 2-torsion is `(Z/2)^{#factors-1}`.

2. **The marked class has order exactly 7.**  On `C_P`,
   `(xy - h)(xy + h) = x^2 f - h^2 = (x-1)^7`.  At `Pt = (1, h(1))` we have
   `xy - h = 0` and `xy + h = 2h(1) ≠ 0`, so all seven zeros of `(x-1)^7` above
   `x=1` on the `Pt` side belong to `xy-h`; at the (single, since `deg f = 5`)
   point at infinity `ord(x) = -2`, `ord(y) = -5`, `ord(h) = -6`, so
   `ord(xy-h) = -7`.  Hence `div(xy-h) = 7(Pt) - 7(infty)` and
   `7 [Pt - infty] = 0`.  The class is nonzero because a genus-2 curve carries
   no rational function with a single simple pole.  So the order is exactly 7.
   *(An identity of the chart, valid for every member — it replaces the
   per-member `Order(D)` check.)*

3. **2-rank is exactly 3, not 4.**  The residual quadratic of `Q(v)` after
   removing `s,t,u` has discriminant
   `Delta = 4 g(-m) / ((m-1)^4 (m+2)^2)` identically on the family
   (`results/claude_ov_lane1_splitall.log`: verified as an identity in
   `Q(m)[w]/(w^2-g(m))`, plus 14 numerical spot checks).  So the quintic splits
   completely iff `g(-m)` is a square.  Combined with `w^2 = g(m)`, that forces
   a rational point on the genus-2 curve
   `H : y^2 = g(m) g(-m)/m^2 = m^6 - 12 m^4 + 36 m^2 - 16`.
   `rank Jac(H) = 1` (`RankBounds = [1,1]`, confirmed by the even-model split
   `Jac H ~ E_a x E_b` with `rank E_a = 1` (cond 288), `rank E_b = 0`
   (cond 36)), so **Chabauty applies**:
   `H(Q) = {(1:±1:0), (-1,±3), (1,±3), (-2,0), (2,0)}`, i.e. `m in {±1, ±2, ∞}`
   (`results/claude_ov_lane1_order112.log`).  The step `y^2 = g(m)g(-m)` →
   `H` divides by `m^2`, so `m = 0` has to be excluded by hand — and it is
   degenerate anyway (`g(0) = 0` forces `w = 0`, i.e. `s = t`).  Every one of
   `m in {0, ±1, ±2}` is in the degenerate set `D` of §1, and `m = ∞`
   (`x(P) = 0`) is not attained on `E(Q)` at all (`y^2 = -16`).  So **no
   nondegenerate member has 2-rank 4.**

   *(Corollary — a clean unconditional negative: the order-112 rung
   `[2,2,2,14]` is **impossible on this family**, for every member, with no
   height bound.  A direct search along `E(Q)` up to `|k| ≤ 400` and a scan of
   1994 family members up to `|n| ≤ 500` also found 0, as they must.)*

If the residual quadratic in `v` is irreducible its two roots `v, v'` give
`1-v^2, 1-v'^2` conjugate; they can only be rational if `v' = -v`, which makes
`1-v^2` a double root of `f` and `disc f = 0` — excluded.  So the factor type
is exactly `[1,1,1,2]`.

**Not a seed claim.**  The `never-claim-a-family-from-its-seed` rule
(`validate-and-record-a-hit`) is respected here: Theorem A is proved by two
polynomial identities on the parametrised family (`x^2 f - h^2 = (x-1)^7` and
`f(1-v^2) = 0` for `v = s,t,u`, both checked in the function field
`Q(m)[w]/(w^2-g(m))`, not at a point), plus one Chabauty computation on a curve
that does not depend on the member.  The seed `k=3, eps=0` is used only to
*name* a base point for the congruence argument of §4, and §5 removes even that
dependence.

---

## 3.  Theorem B (57 explicit members) — `code/claude_ov_lane1_exact.m`

`results/claude_ov_lane1_exact.log`.  Forked 10 ways, `PMAX = 200`,
`KMAX = 30`, wall time < 3 min.

For each member: containment as above, then `N_p = #J(F_p) = L_p(1)` at every
good odd `p` in `[11,200]` (24–38 primes per member).  `J(Q)_tors` injects into
`J(F_p)`, so it divides `g := gcd_p N_p`.

```text
members generated (k = 1..30, eps = 0,1) : 57      (3 degenerate skips: k=1 both, k=2 eps=0)
factor type [1,1,1,2], 2-rank 3, ord(D)=7 : 57/57
gcd_p #J(F_p) = 56                        : 57/57   => TORSION EXACTLY [2,2,14]
strict two-prime End=Z certificate        : 57/57
RM pre-screen: #distinct real-subfield squarefree cores = 9..17 (never constant)
coefficient height:  log10 ht = 9  (k=2)  ...  1177  (k=30, eps=1)
```

The seven members the overnight session had done by exact
`TorsionSubgroup` (`results/claude_ov_lane1_certify.log`, `k = 2..5`) agree
with the gcd bound, and the certifying pairs reproduce **exactly**
(`k=2: (17,31,32)`, `k=3e0: (17,31,64)`, `k=3e1: (17,29,32)`,
`k=4e0: (17,29,32)`, `k=4e1: (17,79,32)`, `k=5e0: (17,31,32)`,
`k=5e1: (43,53,64)`), i.e. the cheap route is validated against the expensive
one on a 7-member positive control.

The strict certificate used is the project's non-negotiable one: at `p0` the
Frobenius quartic `chi` is irreducible over `Q` **and**
`deg MinPoly(pi^n) = 4` for `n = 2..12`; a second prime `q0` with the same
strictness whose splitting field is **linearly disjoint** from `p0`'s
(`[L_{p0} L_{q0} : Q] = [L_{p0}:Q][L_{q0}:Q]`).  Primes are always *scanned*,
never fixed.

---

## 4.  Theorem C (an INFINITE certified subfamily) — `code/claude_ov_lane1_theorem.m`

> **Theorem C.**  Let `P_0 = 3G` and `N = 11`.  For every `n in Z` the member
> `C_{P_0 + 11 n G}` is a genus-2 curve over `Q` with
> `J(Q)_tors = Z/2 x Z/2 x Z/14` **exactly** and `End(J_Qbar) = Z`.
> The set `{P_0 + 11 n G}` is infinite, and the absolute Igusa invariants take
> infinitely many values on it.

Mechanism.  Every ingredient of §1 is a rational function of `P in E(Q)` with
coefficients in `Z[1/2]`.  If `P ≡ P_0` in `E(F_p)` for an odd `p` of good
reduction for `E` at which the construction is defined mod `p` (all the
denominators `x(P)`, `(m+1)^2(m-2)`, `s^2-t^2`, `(s+1)^2`, `(t+1)^2` reduce to
units), then the coefficients of `f_P` are `p`-integral,
`f_P ≡ f_{P_0} (mod p)`, and `C_P` has good reduction at `p`.  Therefore
`#J_P(F_p)` and `chi_{P,p}` depend only on the image of `P` in `E(F_p)`.  Taking

```text
S = {43, 263, 307},      N = lcm_{p in S} ord( Gbar in E(F_p) ) = lcm(11,11,11) = 11,
```

every `P in P_0 + N<G>` has the *same* `#J(F_p)` and the *same* `chi_p` for
`p in S`.  For the base member `k=3, eps=0` (`m = 1/13`, `(s,t,u) = (-13/49, 13/50, -1/2)`):

```text
p = 43  : ord(Gbar)=11  #J(F_p) = 1568   chi = x^4 -  8x^3 +  70x^2 -  344x + 1849   STRICT
p = 263 : ord(Gbar)=11  #J(F_p) = 71064  chi = x^4 +  6x^3 + 310x^2 + 1578x + 69169  not strict
p = 307 : ord(Gbar)=11  #J(F_p) = 88928  chi = x^4 - 18x^3 + 222x^2 - 5526x + 94249  STRICT
gcd over S = 56 ;  (43,307) is a strict linearly-disjoint pair
```

so `J(Q)_tors | 56` and, with Theorem A, `= [2,2,14]`; and `End(J_Qbar) = Z`.
Since `E(Q)` has rank 1, `P_0 + 11<G>` is infinite.  QED.

Checks in the log:

* "built from `Pbar` alone equals the reduction of `f`" verified at every
  `p in S` (`Pbar = (9:9:1)`, `(52:133:1)`, `(52:89:1)`);
* the prediction re-verified by direct recomputation at `n = ±1, ±2`, i.e.
  `k = -19, -8, 14, 25` (coefficient heights 470, 80, 253, 815 decimal digits),
  each in 0.45 s, all four: three rational roots, `g(-m)` not a square,
  identical `(N_p, chi_p)` on `S`, `gcd = 56`;
* the absolute Igusa invariant `j1` of `k = 2..8` (`eps=1`) is pairwise
  distinct, so `P |-> C_P` is finite-to-one and the coset really does give
  infinitely many `Qbar`-isomorphism classes;
* a scan over base members `k0 = 1..12` and all prime sets from `p ≤ 500`
  found `N = 11` minimal overall; `k0 = 9, eps = 0` achieves `N = 12` with only
  **two** primes `S = {61,313}`.

---

## 5.  Theorem C+ (7898 of 7920 congruence classes) — `code/claude_ov_lane1_cover.m`

**This is the answer to "what would it take to go from 57 members to the
generic member".**  Theorem C uses one residue class; nothing stops us from
doing *all* of them.

Two facts shape the computation.

* For a **fixed** `n`, only **finitely many** primes have `ord(Gbar_p) = n`
  (they divide the `n`-th term of the elliptic divisibility sequence of `G`).
  With `PMAX = 800` the only primes with `ord(Gbar_p) = 11` are `43, 263, 307`,
  and three primes are not enough to pin the gcd down to 56 on most classes
  (they give 112 or 224).  So one must combine **many different orders**, not
  many primes of one order.
* Cost is controlled by tabulating `(#J(F_p), chi_p, strict?)` once per
  `(p, r mod n_p, eps)` — `2 * sum_p n_p` Euler factors — and then doing the
  `2L` classes by table lookup.

Set-up: `L = 3960 = 2^3·3^2·5·11`, prime pool = all good `p ≤ PMAX` with
`n_p := ord(Gbar_p)` dividing `L`.

**The certification criterion.**  Write `T = J_P(Q)_tors` and let `g` be the gcd
of `#J(F_p)` over the usable primes.  By Theorem A, `T ⊇ (Z/2)^3 x Z/7` and `T`
has 2-rank **exactly** 3, so `T_2 = Z/2^a x Z/2^b x Z/2^c` with `1 ≤ a ≤ b ≤ c`.
A class is certified when

1. `odd(g) = 7` — then `T_odd = Z/7`;
2. `v_2(g) = 3` **or** at some usable `p` the group `J(F_p)_2` is **elementary
   abelian**, i.e. `v_2(#J(F_p)) = 2-rank J(F_p) = (#irreducible factors of
   f_bar over F_p) - 1`.  In the second case `T` cannot contain an element of
   order 4, so `c = 1` and `T_2 = (Z/2)^3`;
3. two usable primes carry strict Frobenius quartics with linearly disjoint
   splitting fields — then `End(J_Qbar) = Z`.

Criterion 2 in its second form is **strictly stronger than demanding `g = 56`**
and it is free: the factor count of `f_bar` is a by-product of building it.  It
was added after run 1 showed that most failures were `g = 112`, i.e. exactly the
`[2,2,28]`-versus-`[2,2,14]` ambiguity that this test settles.

> **Theorem C+.**  Let `P = kG + eps T` in `E(Q)`.  If `(k mod 3960, eps)` is
> **not one of the 22 classes** printed as `undecided classes` in
> `results/claude_ov_lane1_cover.log` (they are listed verbatim below), then
> `C_P` is a genus-2 curve with `J(Q)_tors = Z/2 x Z/2 x Z/14` **exactly** and
> `End(J_Qbar) = Z`.  That is 7898 of the 7920 classes; each is an infinite set.

**Run 1** (`PMAX = 800`, 26 primes, `results/claude_ov_lane1_cover_p800.log`):

```text
(p, n_p) = (5,5) (7,8) (11,6) (13,10) (19,10) (23,12) (37,20) (43,11) (47,8)
           (59,30) (61,12) (71,36) (239,30) (263,11) (307,11) (313,12) (349,60)
           (373,60) (397,60) (503,36) (571,22) (577,24) (613,36) (647,36)
           (673,30) (769,60)
1312 Euler factors tabulated in 285 s ; 195 distinct strict chi ; sweep 64 s

CERTIFIED (exact [2,2,14] AND strict disjoint End=Z pair) : 7388 of 7920  (93.3%)
  no usable prime (reduction degenerate at all 26)        :   38
  gcd_p #J(F_p) > 56                                      :  444
  gcd = 56 but no linearly disjoint strict pair           :   50
gcd histogram : 56 -> 7438,  112 -> 396,  224 -> 12,
                2184 -> 6, 3136 -> 6, 334712 -> 10, 337680 -> 10, 351456 -> 4
```

The large gcd values are classes with a *tiny* usable-prime pool (one or two
primes), not evidence of extra torsion.  The 8 genuine members with `|k| ≤ 40`
that run 1 fails to cover are `k = ±10, ±28, ±31, ±32`, all with `eps = 1`; and
all eight are already known to be exactly `[2,2,14]` by Theorem B, which is free
to use *any* prime.  For instance `k=10, eps=1` needs the pair `(11,17)` and
`n_17 = 26 ∤ 3960`, so `p = 17` is not in the pool.  **This is a limitation of
the prime pool, not a mathematical obstruction.**

**Run 2** (`PMAX = 5000`, 42 primes, criterion 2 in its strong form,
`results/claude_ov_lane1_cover.log`) — **this is the run to quote**:

```text
42 primes, orders n_p in {5,6,8,10,11,12,18,20,22,24,30,36,55,60,66,90,110,120,132,180,198}
3874 Euler factors tabulated in 605 s ; 754 distinct strict chi ; sweep 104 s

CERTIFIED (exact [2,2,14] AND strict disjoint End=Z pair) : 7898 of 7920  (99.72%)
  no usable prime at all                                  :   10 classes
  odd(g) != 7 or 2-part not settled                       :   12 classes
  gcd = 56 but no linearly disjoint strict pair           :    0 classes
gcd histogram : 56 -> 7840,  112 -> 68,  168 -> 2

  no usable prime :  (0,0) (0,1) (1,0) (1,1) (2,0) (1980,0) (1980,1)
                     (3958,0) (3959,0) (3959,1)
  not settled     :  (358,0,112) (900,0,112) (1080,0,112) (1439,1,112)
                     (1979,0,112) (1979,1,168) (1981,0,112) (1981,1,168)
                     (2521,1,112) (2880,0,112) (3060,0,112) (3602,0,112)

EVERY genuine member with |k| <= 40 is certified.
```

**58 of the 68 classes with `g = 112` were saved by the elementary-abelian
2-part test** (7898 − 7840 = 58); without it they would all have been left open,
as in run 1.

Three structural checks on that output.

* **The undecided set is stable under `k -> -k`**, as it must be: `x(-P) = x(P)`
  so `m(-P) = m(P)`, and `w` flips sign, which merely swaps `s` and `t` — the
  same curve.  Indeed the 22 classes pair up as
  `358 ↔ 3602`, `900 ↔ 3060`, `1080 ↔ 2880`, `1439 ↔ 2521`, `1979 ↔ 1981`,
  `1 ↔ 3959`, `2 ↔ 3958`, with `0` and `1980` self-paired.
* **Eight of the ten "no usable prime" classes are the degenerate classes**
  `(0,0), (0,1), (1,0), (1,1), (2,0), (3958,0), (3959,0), (3959,1)` — those can
  never be decided by reduction (see below).
* **The other two are `k ≡ 1980 = L/2`, and that is explained, not mysterious.**
  For every pool prime, `3960·Gbar = O`, so `1980·Gbar` has order dividing 2.  At
  the pool primes with `8 | n_p` (`p = 7, 47, 577, 1201, 1487, 2833`) it has
  order exactly 2, and there `x^2+4x-8` is irreducible mod `p`, so
  `E(F_p)[2] = {O, Tbar}` and `1980·Gbar = Tbar` is forced — degenerate for
  `eps = 0`, and `Tbar + Tbar = O` for `eps = 1`.  Unlike the eight degenerate
  classes, this one is an artefact of `L`: a pool with `16 | L` would split it.

The 12 remaining classes are prime-pool-limited exactly as in run 1; ten sit at
`g = 112` (`[2,2,14]` vs `[2,2,28]`, undecided because no usable prime had an
elementary-abelian 2-part) and two at `g = 168` (odd part 21, so `Z/21` is not
yet excluded).  Nothing suggests they are different in kind.

**Validation.**  For `k = 3, 5, 8, 14, 25` the pairs `(#J(F_p), chi_p)` read off
the **class table** were compared with those of the actual rational curve `C_P`
at 30–35 usable primes each: identical in every case, `gcd = 56` each time.
This is the check that the Key Lemma is being applied correctly.

**Why eight classes can never be removed by this method.**  A class is invisible
to every prime in the pool exactly when its reduction is a degenerate point at
all of them, and `P ≡ D (mod L<G>)` does precisely that for every `p` with
`n_p | L`.  Enlarging the pool enlarges `L` and shrinks the excluded set to
`P ≡ D (mod L'<G>)`, but never empties it: the reduction of `D` itself is always
degenerate.  So the residue-class method can reach density `1 - 8/(2L)` — as
close to "all" as one likes — but never literally all.  **§8 says what would.**

---

## 6.  Where the eleven known rational points of `R` sit; the second-component question

`code/claude_ov_lane1_pencil.m`, `results/claude_ov_lane1_pencil.log`.

**5 of the 11 lie on the `S3`-orbit of the `u = -1/2` curve** — they are exactly
the ones with a coordinate equal to `-1/2`, and each is a small member:

```text
(-1/2, -13/49, 13/50)     = k=3, eps=0  (m = 1/13)
(-15/8, -15/19, -1/2)     = k=3, eps=1  (m = -9/5)
(-164/297, -1/2, 164/361) = k=4, eps=0  (m = 16/41)
(-511/61, -511/625, -1/2) = k=4, eps=1  (m = -98/73)
(-10, -10/7, -1/2)        = k=2, eps=1  (m = -4/5)
```

Note the perfect match with §1: the members with `|k| ≤ 4` that are *not*
degenerate are exactly `k=2 eps=1`, `k=3` both, `k=4` both — five of them, and
all five are among the eleven known points.  The eleven were found by a height
search on `R`, so this is what one expects, and it is a consistency check on the
parametrisation rather than new information.

**6 do NOT lie on the family**, including the `RM(sqrt2)` curve:

```text
(-5, -15/8, -15/22)        (-3, -3/4, -3/5)  <- RM(sqrt 2)
(-5/18, -10/49, 4/17)      (-4/9, -4/25, 4/17)
(-165/41, -33/16, -165/289) (-17/50, -34/189, 34/121)
```

**What that says.**  `C_c = R ∩ {u=c}` is a `(3,3)` curve; via `p=s+t, q=st` it
is the double cover of the plane cubic `E_c` branched at the 6 points
`E_c ∩ {p^2=4q}`, so generically `g(C_c) = 2·1 - 1 + 6/2 = 4`.  The genus can
drop only where `E_c` is singular or the branch divisor is non-reduced, and
(`results/claude_ov_lane1x_specialc.log`) those two discriminants — degree 28
and degree 26 in `c` — have **rational** roots only at

```text
disc_k(E_c pencil) :  c in {0, -1, -1/2}      (other factors 4c^2+2c-1, 8c^3+4c^2-4c-1)
disc_p(branch)     :  c in {1, 0, -1, -1/2}   (other factor a degree-6 irreducible)
```

`c = 0` and `c = -1` are chart-degenerate (each fibre splits into 3 components,
`v = 0` is the marked root `x = 1` and `v = -1` blows up `G7`), `c = 1` has
genus 3, `c = -1/2` has genus 1.  Hence **for every rational `c` outside
`{0, ±1, -1/2}` the fibre has geometric genus 4** and, by Faltings, finitely
many rational points.  Two independent confirmations:

* the **17 distinct fibre parameters** occurring as coordinates of the six
  off-family points all give genus-4 fibres (each listed in the log);
* an **exhaustive genus scan** over all 1957 rational `c = n/d` with
  `|n|, d ≤ 40`: 1955 of 1957 have genus exactly 4, and the two exceptions are
  `c = 1` (genus 3) and `c = -1/2` (genus 1).

> So the `u = -1/2` curve (and its two `S3`-images) is the ONLY fibre of the
> `u = const` pencil with infinitely many rational points.  A second infinite
> component of `R = 0` would have to be a **multisection** of the pencil.

The two points sharing `u = 4/17` sit on one genus-4 fibre; that fibre has ≥ 2
rational points but only finitely many.  Nothing here contradicts the family:
the six off-family points are simply sporadic points on genus-≥3 fibres, which
is exactly what Bombieri–Lang predicts for a surface of general type with one
rational curve on it.

**Where the RM points sit.**  The known `RM(sqrt2)` point `(-3,-3/4,-3/5)` is
*off* the family (its three fibres are `c = -3, -3/4, -3/5`, all genus 4).  All
57 tested family members have a scattering real-subfield-disc signature (9–17
distinct squarefree cores over 24–38 primes) and a strict `End=Z` certificate;
the control `(-3,-3/4,-3/5)` returns the constant signature `{2}`.  No RM member
of the family has been found.  *This is not a proof that the family contains no
RM member* — the RM locus in a 1-parameter family of abelian surfaces is a
countable union of points and could meet it — but Theorem C+ gives an `End=Z`
certificate on 7898 of the 7920 congruence classes unconditionally, which is
stronger than a screen.

---

## 7.  Global geometry of the surface — `code/claude_ov_lane1_geometry.m`

The overnight version of this script crashed in `RationalPointsGeneric`.
**Fixed** (not dropped): `RationalPoints` refuses a scheme living in a *product*
of projective spaces (`Argument must lie in affine or ordinary projective
space`); pushing each component through `SegreEmbedding` into `P^7` first makes
point enumeration work, with a `GroebnerBasis` fallback.  The script now runs to
`LANE1_GEOMETRY_DONE` (`results/claude_ov_lane1_geometry.log`).

* `Xbar ⊂ P^1 x P^1 x P^1` is an irreducible surface of tridegree `(3,3,3)`;
  `K_{(P^1)^3} = O(-2,-2,-2)`, so by adjunction `K_Xbar = O(1,1,1)|_Xbar`,
  which is **ample**.  `p_g = h^0(O(1,1,1)) = 8` (the restriction sequence is
  exact on `H^0` because `H^1(O(-2,-2,-2)) = 0` by Künneth) and
  `K^2 = (H1+H2+H3)^2·(3H1+3H2+3H3) = 18`.
* `Sing(Xbar)` is **0-dimensional**: **14** reduced points, listed in the log by
  their Segre coordinates and defining ideals.  Decoding them: **8 affine**
  points — `(0,0,0)`, the `S3`-orbit of `(-1/2,0,0)` (3 points), and the four
  points `(-1,-1,-1)`, `(-1,-1,1)`, `(-1,1,-1)`, `(1,-1,-1)` (the `(±1,±1,±1)`
  with at least two `-1`s) — and **6 at infinity** (one or two coordinates
  `= ∞`).  Projected to the `u`-line the affine singular ideal is
  `u^6 (u-1)(u+1)(u+1/2)^3`.
* So `Xbar` is a surface with ample canonical class and isolated singularities:
  of **general type** provided those 14 singularities are canonical, which was
  **not** verified (the singularity types were not classified).  Under
  Bombieri–Lang its rational points should lie on finitely many curves; we have
  exhibited one (`u = -1/2` and its two `S3`-images) and shown in §6 that no
  other *fibre* of the `u`-pencil qualifies.
* **The `S3` quotient and the precise obstruction.**  In the symmetric
  coordinates `R` is a *quadratic* in `e3`:
  `R = 2(e2-2)e3^2 + (...)e3 + (...)`, so `X_e := X/S3` is the double cover of
  the `(e1,e2)`-plane branched along
  `Disc_{e3}(R) = 16 e1^2 e2 + 4 e1^2 + 8 e1 e2^2 + 24 e1 e2 + 4 e1 + 8 e2^2 + 8 e2 + 1`,
  an **irreducible plane cubic** whose projectivisation is singular at
  `(-1/2:0:1)` and has geometric genus 0.  (The projective branch divisor of the
  double plane is that cubic together with the line at infinity — a quartic with
  several nodes, so `X_e` is a singular/weak degree-2 del Pezzo and hence
  rational over `Qbar`; *that last sentence is a classification assertion, not a
  computed parametrisation*.)  All eleven known points of `X` have square
  `Disc` value, as they must.
  The obstruction is then exactly this: rational points of `X_e` are cheap, but
  a point of `X_e(Q)` gives a `[2,2,14]` curve only if
  `T^3 - e1 T^2 + e2 T - e3` **splits completely over `Q`**, and that splitting
  cover *is* `X` itself, which is of general type.  A systematic attack on that
  covering condition is the natural next move for a second component
  (**untouched**).

---

## 8.  What would be needed for "EVERY member", and what would be needed for the surface

**(a) Every member of the family.**  After Theorem C+ the gap is (i) the
prime-pool-limited classes, which shrink as the pool grows and cost only compute,
and (ii) the eight classes `P ≡ D (mod L<G>)`, which reduction can **never** see
however large the pool.  For (ii) a different argument is needed.  What is
*not* enough:

* a generic-fibre computation over `K = Q(m)[w]/(w^2-g(m))` gives
  `J_eta(K)_tors`, and specialisation of an abelian scheme injects
  `J_eta(K)_tors ↪ J_P(Q)_tors` — the wrong direction.  It cannot bound the
  special fibres, where torsion *jumps up*.  (It does immediately give
  `J_eta(K)_tors = [2,2,14]` exactly, from Theorem A plus any one member with
  `gcd = 56` — but that is a weaker statement than it looks.)

What *is* enough: rule out each way the group can grow, uniformly.  With
Theorem A (`2`-rank exactly 3, marked class of order exactly 7) the possible
proper overgroups of `[2,2,14]` inside `J(Q)_tors` are constrained to:

1. `(Z/2)^4 x Z/7` (2-rank 4) — **already impossible**, Chabauty on `H`, §2.3;
2. `Z/2 x Z/2 x Z/28` and higher 2-power: requires one of the **7** nonzero
   rational 2-torsion classes to be divisible by 2 in `J(Q)`.  This is the one
   real remaining computation.  Concretely: for `f = (x-r1)(x-r2)(x-r3)·q(x)`
   with `r_i = 1-s^2, 1-t^2, 1-u^2`, 2-divisibility of the class
   `[(r_i,0) - (r_j,0)]` is a set of explicit square-class conditions in the
   differences of the roots of `f` (Cassels' `x - a` map to `L*/L*^2`,
   `L = Q[x]/f`); imposing them on the family gives, for each of the 7 classes,
   a cover of `E` whose rational points must be determined (each is a curve in
   `m`, `w`, and one or two new square roots — expect genus ≥ 2, hence a
   Chabauty/Mordell–Weil-sieve job like the one already done for `H`);
3. odd growth `Z/49`, `Z/21`, `Z/35`, `(Z/7)^2`: killed on any single class by
   one prime with the corresponding non-divisibility of `#J(F_p)` — already
   automatic in the `gcd = 56` computation for every class Theorem C+ covers,
   and for the eight residual classes it would follow from the same 2-divisibility
   machinery plus one more explicit cover for `49`.

So: **one genuinely new computation (item 2, the seven 2-divisibility covers)
would upgrade Theorem C+ to "every nondegenerate member"**, because the eight
residual classes would then be handled by Theorem A + item 1 + item 2 + item 3
without any reduction argument at all.  That is the recommended next move for
this lane and it is a self-contained, well-posed job.

Two remarks that put the size of the gap in perspective.

* The elementary-abelian 2-part test of §5 is a *reduction-based* substitute for
  item 2, and it works for every class that any prime can see — it certified 58
  classes at `g = 112` in run 2.  Item 2 is needed only for the eight classes
  reduction cannot see.
* Those eight classes are `k ≡ 0, ±1, ±2 (mod 3960)`, and their smallest
  nondegenerate members are `|k| ≈ 3958`.  Coefficient heights on this family
  grow like `log10 ht ≈ 1.31 k^2` (fitted on `k = 8, 14, 25, 30` → 80, 253, 815,
  1177 digits), so the smallest undecided member has roughly `2 x 10^7` decimal
  digits.  The gap is real as a *theorem* but contains no curve anyone will ever
  write down.

**(b) The surface.**  Whether `X` has a *second* infinite component is open.
§6 rules out any other fibre of the `u`-pencil; a second component would be a
multisection.  Two concrete attacks, both untouched: (i) search for rational
curves on `Xbar` of low bidegree other than `(0,*,*)`-type fibres, e.g. by
looking for `S3`-invariant curves; (ii) attack the covering condition on the
rational surface `X_e` directly (§7).

---

## 9.  Negatives, with their bounds

| statement | bound / method | result |
|---|---|---|
| order 112 (`[2,2,2,14]`) on this family | **Chabauty on the genus-2 curve `H` (rank 1)** — complete, no height bound | **impossible for every member** (only degenerate `m in {0,±1,±2}`) |
| order 112 along `E(Q)` | direct, `\|k\| ≤ 400` | 0 hits (consistent) |
| order 112 over family members | `\|n\| ≤ 500`, 1994 members | 0 hits (consistent) |
| SPLITALL over the whole surface (not just the family) | Lane 2's C sieve, `1.895e13` pairs, `H ≤ 3200` | 0 candidates |
| extra torsion on a member | `gcd_p #J(F_p)` over 24–38 primes, all 57 members `k ≤ 30` | always exactly 56 |
| extra torsion on a residue class | 42 primes `p ≤ 5000`, all 7920 classes mod 3960, gcd + elementary-abelian 2-part test + strict pair | **7898 CERTIFIED**; 22 open, of which 8 are permanently invisible to reduction (§5) |
| RM member of the family | real-subfield-disc screen, 57 members, 24–38 primes each | none (all scatter) |
| second *fibre* with infinitely many points | discriminants of the `(p,q)`-cubic pencil (deg 28) and of the branch sextic (deg 26); exhaustive genus scan `\|n\|,d ≤ 40` (1957 fibres) | none besides `c = -1/2` |
| second component of `X` in general | — | **not attempted** (multisection question open) |

---

## 10.  What I did NOT check (honest list)

* **No exact `TorsionSubgroup(J)` beyond `k ≤ 5`** (7 members).  Exactness for
  everything else rests on `containment (Theorem A) + gcd_p #J(F_p) = 56`, which
  *is* a complete proof, but is not the same computation.  It was validated
  against exact `TorsionSubgroup` on those 7 members.
* **No `Sage`/Lombardo endomorphism-algebra test** anywhere in this lane.
  Genericity rests entirely on the strict two-prime Frobenius certificate.
* **No claim that every member is `End = Z` or exactly `[2,2,14]`** — 22 of the
  7920 classes mod 3960 are still open (8 of them permanently, as far as
  reduction goes), and the route to closing those eight is §8(a).  The other 14
  should fall to a larger prime pool; that was not run.
* **No claim that the family exhausts the `[2,2,14]` curves.**  Six of the
  eleven known points of `R` are off it, and a second (multisection) component
  is not excluded.
* **The 14 singular points of `Xbar` were not classified**, so "general type" is
  asserted from ampleness of `K_Xbar` plus 0-dimensionality of `Sing` and is
  conditional on those singularities being canonical.
* **`X_e` was not explicitly parametrised**; "singular degree-2 del Pezzo, hence
  rational over `Qbar`" is read off the branch-curve degree and singularity, not
  computed.
* **Chabauty on `H`** was taken from Magma's `Chabauty` with the Mordell–Weil
  group supplied; `MW(Jac H) = Z/2 + Z/6 + Z` was not independently **saturated**
  (the rank-1 bound `RankBounds = [1,1]` and the isogeny split were both
  verified, so only index-saturation of the free part is unaudited).
* **No search for rational curves on `Xbar` other than fibres of the `u`-pencil.**
* The RM screen is a *screen*: it does not prove no member is RM.

---

## 11.  Exact resume commands

All Magma **must** go through the global slot limiter.

```bash
cd /home/claude/torsion_jac

# 1. exact torsion + certificates for members k = 1..KMAX  (forks 10 ways)
code/claude_magma_slot.sh -b KMAX:=30 PMAX:=200 NCH:=10 MemGB:=6 \
    code/claude_ov_lane1_exact.m > results/claude_ov_lane1_exact.log 2>&1 &
#    to go further:  KMAX:=60  (heights ~4700 digits; still minutes)

# 2. the congruence/infinite-family theorem (base-member + prime-set search)
code/claude_magma_slot.sh -b KSCAN:=12 PMAX:=500 MemGB:=8 \
    code/claude_ov_lane1_theorem.m > results/claude_ov_lane1_theorem.log 2>&1 &
#    KSCAN:=0 short-circuits the scan and installs the known optimum.

# 3. THE COVERING (Theorem C+).  THIS IS THE RUN THAT PRODUCED THE 7898/7920.
#    12 min wall (605 s tabulate + 104 s sweep) on one core.
code/claude_magma_slot.sh -b PMAX:=5000 LTARGET:=3960 NMAX:=200 MemGB:=8 \
    code/claude_ov_lane1_cover.m > results/claude_ov_lane1_cover.log 2>&1 &
#    To close the 12 prime-pool-limited classes, raise PMAX (the pool of primes
#    with n_p | 3960 continues to 21961, 23399, 25343, 25633, 28753, 32401,
#    34319, 35069, 37441 with n_p = 9,30,99,180,120,90,55,55,36): try
#      PMAX := 40000, NMAX := 200        (~8100 Euler factors, ~1 h)
#    To split the two k = 1980 classes, use an L divisible by 16, e.g.
#      LTARGET := 7920, PMAX := 40000    (15840 classes; the sweep is a lookup)
#    Cost = 2*sum_p n_p Euler factors (tabulate) + 2*L*#S lookups (sweep).
#    NOTHING will ever close the eight degenerate classes -- see section 8(a).

# 4. the pencil: which fibres can have infinitely many points  (forks)
code/claude_magma_slot.sh -b H:=40 NCH:=12 MemGB:=6 \
    code/claude_ov_lane1_pencil.m > results/claude_ov_lane1_pencil.log 2>&1 &

# 5. surface geometry (the RationalPointsGeneric crash is fixed via SegreEmbedding)
code/claude_magma_slot.sh -b MemGB:=24 \
    code/claude_ov_lane1_geometry.m > results/claude_ov_lane1_geometry.log 2>&1 &

# from the overnight session, unchanged and still valid:
code/claude_magma_slot.sh -b KLIST:="2,3,4,5,6,7" code/claude_ov_lane1_certify.m   # exact TorsionSubgroup
code/claude_magma_slot.sh -b code/claude_ov_lane1_order112.m                        # the Chabauty kill
code/claude_magma_slot.sh -b code/claude_ov_lane1_splitall.m                        # Delta = 4 g(-m)/(...)
```

**Not yet written** (the §8(a) job): a script that, for each of the 7 nonzero
rational 2-torsion classes of the family, writes down the 2-divisibility
square-class conditions over `Q(m)[w]/(w^2-g(m))`, forms the resulting cover of
`E`, and determines its rational points.  That is the single computation that
would turn Theorem C+ into "every nondegenerate member".

---

## 12.  Files

Code (this session): `code/claude_ov_lane1_cover.m`.
Code (resumed session): `code/claude_ov_lane1_exact.m`,
`code/claude_ov_lane1_theorem.m`, `code/claude_ov_lane1_pencil.m`,
and the `SegreEmbedding` repair of `code/claude_ov_lane1_geometry.m`.
Code (overnight): `code/claude_ov_lane1_family.m`,
`code/claude_ov_lane1_certify.m`, `code/claude_ov_lane1_fibration.m`,
`code/claude_ov_lane1_order112.m`, `code/claude_ov_lane1_splitall.m`,
`code/claude_ov_lane1_genusscan.m`, and the `claude_ov_lane1x_*` variants
(PARI/GP prototypes and their Magma verifications).

Logs: `results/claude_ov_lane1_cover.log`,
`results/claude_ov_lane1_exact.log`,
`results/claude_ov_lane1_theorem.log` + `..._theorem_fixed.log`,
`results/claude_ov_lane1_pencil.log`, `results/claude_ov_lane1_geometry.log`,
`results/claude_ov_lane1_order112.log`, `results/claude_ov_lane1_splitall.log`,
`results/claude_ov_lane1_certify.log`, `results/claude_ov_lane1_family.log`,
plus the overnight `results/claude_ov_lane1x_*.log`.

Commits: `981cdc7` (57 members), `d36bb18` (Theorem C), `e824c15` (no second
vertical component), and this session's covering + note commit.
