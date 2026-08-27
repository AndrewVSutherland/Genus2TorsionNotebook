# Geometrically simple genus-2 torsion frontier after the order-96 record

Date: 2026-07-18

## Executive summary

The collaborator's curve is a genuine new record:

```text
J(Q)_tors = [2,2,2,12],    #J(Q)_tors = 96.
```

It is geometrically simple.  Exact torsion is independently bounded by

```text
#J(F_31) = 864,   #J(F_37) = 1248,   gcd(864,1248) = 96,
```

and geometric simplicity is certified at both `p=37` and `p=73` by an
irreducible Frobenius quartic whose root powers retain degree four for every
`n=2,...,12`.

This improves Elkies' presently known simple order-80 record by 20 percent.
The larger HLP torsion groups of orders 96 and 128 are geometrically split,
so they do not compete in the simple category.

The repository audit also produced three immediate new realizations.  Exact
geometrically simple examples now exist for

```text
[2,4,4],       order 32,
[2,2,2,8],     order 64,
[2,2,4,4],     order 64.
```

These had already been found as exact torsion groups in the project, but the
necessary geometric-simplicity calculation had never been done.  All three
pass the strict root-power test at two independent good primes.

## Scope of the audit

The workspace contains 3,289 non-git files, including 534 code files, 157
research notes, 604 data files, and 99 result logs in the principal
`torsion_jac` tree.  The audit began with all nine PDFs (505 pages) and
their matching TeX sources, then followed the target-specific references into
the remaining TeX, Markdown, code, data, and result files.

The important PDF/TeX conclusions are:

- `notes.pdf` and `notes.tex` explicitly construct the
  `A(2,2,2,12)` order-96 cover but say that no example had been found.  The
  collaborator's curve fills exactly this documented gap.
- `rational_parametrizations_catalog_2026_07_12.pdf/.tex` is the best
  repository-wide status catalogue.  Its distinction between a visible
  subgroup and generic exact torsion must be retained.
- Elkies, *Families of genus-2 curves with 5-torsion* (ConM 796, 2024),
  gives the previous geometrically
  simple `[2,2,2,10]` record of order 80 and an absolute-simplicity
  certificate.
- `hlp_2_2_24_examples.pdf/.tex` gives ten exact `[2,2,24]` examples, all on
  a split/glued locus.
- `a12_parameterization.pdf/.tex` and `a2_12_resolvent.pdf/.tex` give the
  rational `A(12)` chart and the `A(2,12)` resolvent, but not `[2,24]`.
- the Stoll and Coleman PDFs are computational background rather than new
  frontier realizations.
- the remaining chapters of the 386-page LUCANT volume do not add further
  genus-2 rational-torsion realizations beyond the included Elkies article.

The primary external checks were Elkies' paper,
[Howe--Leprévost--Poonen](https://math.mit.edu/~poonen/papers/large.pdf),
[Howe's split-family table](https://ewhowe.com/genus2.html), and the current
bounded LMFDB collection.  The novelty language in this note means
"presently known after this audit," not a classification theorem.  In
particular, absence from LMFDB is supporting evidence, not a proof of global
novelty.

## Certificate standard

Several older notes call an irreducible Frobenius quartic a geometric-
simplicity certificate.  That is not enough: a power of Frobenius can acquire
a quadratic minimal polynomial, corresponding to geometric splitting after a
finite extension.

Every new claim in this audit uses both:

1. exact `TorsionSubgroup(Jacobian(C))`; and
2. a good prime for which, for a Frobenius root `pi`,

   ```text
   Degree(MinimalPolynomial(pi^n)) = 4,  n=2,...,12.
   ```

Whenever convenient, two independent good primes are recorded.  The
Frobenius polynomial reconstructed from point counts is also checked against
Magma's built-in `LPolynomial`.

## The new order-96 record

The generalized input model is

```text
y^2 + (x^2+1)y =
  756900*x^6 + 737595570*x^5 + 150572203590*x^4
  - 15854483576121*x^3 - 530648977741620*x^2
  + 32014154874551031*x + 830742747091037849.
```

Completing the square gives exactly

```text
y^2 = 3027600*x^6 + 2950382280*x^5 + 602288814361*x^4
      - 63417934304484*x^3 - 2122595910966478*x^2
      + 128056619498204124*x + 3322970988364151397.
```

The sextic factors as

```text
1740^2
*(x-123109/1740)*(x-59)*(x+299/12)
*(x+39)*(x+519)*(x+75593/145),
```

as expected from full rational 2-torsion.

The two root-power witnesses are

```text
p=37: T^4 - 4*T^3 + 30*T^2 - 148*T + 1369,
p=73: T^4 + 8*T^3 + 30*T^2 + 584*T + 5329.
```

The reusable verifier is `code/verify_record_22212_order96.m`.

The previous record is Elkies' `[2,2,2,10]` example of order 80; see
[Elkies, Families of genus-2 curves with 5-torsion](https://doi.org/10.1090/conm/796/16001).
HLP's order-96 `[2,2,24]` and order-128 `[2,2,4,8]` constructions are
`(2,2)`-isogenous to products and hence geometrically nonsimple.

## Three additional realizations found by the audit

### Exact `[2,2,4,4]`, order 64

A compact model is

```text
y^2 = x*(x+1296)*(x+3249)*(x+4096)*(x+17424)
    = x^5 + 26065*x^4 + 173387808*x^3
      + 414985109760*x^2 + 300512487407616*x.
```

Magma gives exact torsion `[2,2,4,4]`.  The witnesses are

```text
p=37: T^4 - 4*T^3 - 2*T^2 - 148*T + 1369,
p=47: T^4 + 8*T^3 + 30*T^2 + 376*T + 2209.
```

This is tuple `(36,57,64,132)`, row 26629 of
`paper/scripts_and_data/tor2244.txt`.

### Exact `[2,2,2,8]`, order 64

The smallest primitive-height geometrically simple member of the current
bank is

```text
y^2 = x*(x+1)*(x+3025)*(x+9801)*(x+15625)
    = x^5 + 28452*x^4 + 230082726*x^3
      + 463480444900*x^2 + 463250390625*x.
```

Magma gives exact torsion `[2,2,2,8]`.  The witnesses are

```text
p=41: T^4 + 4*T^3 + 6*T^2 + 164*T + 1681,
p=47: T^4 + 4*T^3 + 30*T^2 + 188*T + 2209.
```

The four smaller primitive tuples in the bank have automorphism group of
order four and two rational degree-2 elliptic subcovers; they are bielliptic.

### Exact `[2,4,4]`, order 32

The first regenerated `M(2,4,4)` pilot curve already works:

```text
y^2 = x*(x+16)*(x+(644/799)^2)*(x^2+x+16).
```

An integral square-twist model is

```text
y^2 = 407555836801*x^5 + 7193217102753*x^4
      + 17542840688944*x^3 + 112806866289408*x^2
      + 67780576546816*x.
```

Magma gives exact torsion `[2,4,4]`.  The witnesses are

```text
p=97:  T^4 + 8*T^3 + 46*T^2 + 776*T + 9409,
p=103: T^4 - 4*T^3 + 46*T^2 - 412*T + 10609.
```

Its automorphism group has order two and it has no rational degree-2
subcovers.  Thus the elliptic-fiber-product parameterization is not itself a
geometrically split construction.

Full certificates and commands are in
`notes/certify_frontier_existing_banks_20260718.md`; the reusable scanner is
`code/certify_frontier_existing_banks.m`.

## The pre-test top ten and what happened

| initial rank | group | executed test | outcome |
|---:|:---|:---|:---|
| 1 | `[2,4,4]` | regenerated the exact bank and applied two-prime root-power tests | **realized; remove from open list** |
| 2 | `[2,2,4,4]` | scanned `tor2244.txt` with exact torsion first | **realized; remove from open list** |
| 3 | `[2,2,2,8]` | scanned `tor2228.txt`, including bielliptic diagnostics | **realized; remove from open list** |
| 4 | `[2,24]` | 54,834 `A(8)` triples | 101 gained the extra 2-rank; none passed the rational-3 prefilter |
| 5 | `[2,6,6]` | 129,030 contact-6 parameters plus exact split controls | exact split targets reproduced; no simple target; good-open locus empty mod 5 and 7 |
| 6 | `[3,12]` | 7,569 `M(2,12)` parameters | exact split target reproduced; twelve simple-screen survivors all had `[12]`; good-open locus empty mod 5 |
| 7 | `[6,12]` | reran the Prym verifier | Prym rank is exactly one with trivial torsion; route remains in Prym-Chabauty range |
| 8 | `[2,2,2,24]` | audited source halves and closed the record's rational Richelot graph | order-12 generator is not divisible; no target in the 18-vertex graph |
| 9 | `[2,2,4,12]` | same halving/Richelot audit | exactly one source 2-class already has a half; no additional target |
| 10 | `[2,2,4,8]` | HLP normalization plus three bounded nonsplit searches | split positive controls pass; nonsplit searches remain cold |

There was one useful correction in the `[6,12]` rerun: the curve described in
`contact6_m612_nearmiss_recertify.m` as a `[2,12]` near miss actually has exact
torsion `[2,2,6]`.  Its `p=7` root-power simplicity certificate is valid, but
the torsion label in the old prose is not.

## The record's complete rational Richelot component

All 15 rational Richelot kernels of the record source were enumerated.  Every
manual pairing has nonzero determinant, gives a smooth codomain, and matches
exactly one Magma codomain over `Q`, including its twist.

The immediate layer consists of

```text
12 curves with exact [2,12],
 3 curves with exact [2,2,6].
```

Depth two adds two new `[2,2,6]` curves; depth three adds nothing.  The
component is therefore closed with 18 vertices:

```text
1 x [2,2,2,12],   12 x [2,12],   5 x [2,2,6].
```

The source, the two depth-two vertices, and the three immediate `[2,2,6]`
vertices form `K_(3,3)`; the twelve `[2,12]` vertices are leaves attached to
the source.  Every vertex is geometrically simple because simplicity is
isogeny invariant.

No vertex has `[4,12]`, `[2,4,12]`, `[2,24]`, `[2,2,24]`,
`[2,2,2,24]`, or `[2,2,4,12]`.  The fixed record point is therefore exhausted;
future Richelot work must vary the source in the `A(2,2,2,12)` family.

See `notes/record_22212_richelot_bfs_2026_07_18.md` and
`code/record_22212_richelot_bfs.m`.

## Revised ranked top ten that remain open

The ranking balances arithmetic evidence, maturity of the available
construction, and probability of leaving a split/Humbert locus.  The last
two entries are deliberately ambitious because they would more than double
the pre-existing record.

### 1. `[2,24]`, order 48

**Why it is plausible.**  Exact split examples exist, the `A(8)` chart makes
the extra 2-part abundant, and the finite-field work finds no global local
obstruction.

**Strategy.**  On the `W`-split `A(8)` sublocus, eliminate the cubic-contact
3-torsion equations after saturating singular and split components.  Compute
the residual components and their genera rather than extending the height
box.  In parallel, apply the same 2-primary test to Richelot codomains while
varying the new `A(2,2,2,12)` source.

**Test run.**  Among 54,834 triples, 101 acquired the required extra 2-rank;
none survived the ten-prime rational-3 prefilter.  The fixed record's complete
Richelot component also contains no `[2,24]` vertex.

### 2. `[2,6,6]`, order 72

**Why it is plausible.**  Exact split examples are easy, and a geometrically
simple exact `[6,6]` core is already known.

**Strategy.**  Adjoin the discriminant square root of one contact-6 quadratic
to the cubic-contact `[6,6]` core, eliminate contact variables, and saturate
the bielliptic factors.  Projectivize and resolve the simultaneous 5- and
7-adic boundary components; determine whether every surviving component is
forced back onto the split locus.

**Test run.**  Exact `[2,6,6]` was reproduced, but every known target is split.
A simple-screen run checked 129,030 parameters and exact-tested 20 survivors
without a hit.  The good-open target has no residues at either 5 or 7.

### 3. `[3,12]`, order 36

**Why it is plausible.**  Exact split examples exist, and `M(12)` is a
rational surface with an explicit marked order-12 class.

**Strategy.**  Construct the exact independent 3-contact cover over
`M(2,12)`, saturate the discriminant and Humbert factors, decompose the
residual cover, and carry the mod-5 boundary analysis one blowup deeper.  Use
Mordell--Weil or Chabauty on any surviving genus-one or genus-two component.

**Test run.**  Of 7,569 pairs, twelve simple-screen survivors had exact
torsion `[12]`; the exact `[3,12]` control is split.  There are no good-open
target residues modulo 5.

### 4. `[2,2,24]`, order 96

**Why it is plausible.**  HLP supplies ten exact rational examples, so the
group is arithmetically feasible and already matches the new record order.
Only geometric simplicity is missing.

**Strategy.**  Use an HLP point as a marked boundary point of the full
`A(12)` halving/extra-2 incidence variety.  Compute tangent directions
transverse to all visible Humbert components, then choose a two-parameter
plane on which one halving block has a rational section.  Search the remaining
low-degree cover and certify exact torsion plus root-power simplicity.

**Test run.**  `hlp_2_2_24_points.m` reverified all ten exact split
`[2,2,24]` examples.  The fixed record's complete Richelot component has no
`[2,2,24]` vertex, so a deformation rather than another graph step is needed.

### 5. `[6,12]`, order 72

**Why it is plausible.**  Exact split examples exist and this is the most
developed global nonsplit route in the repository.

**Strategy.**  Certify the characteristic-zero map from the degree-24 signed
cover to `E8`, certify the bigonal correspondence, then run a precision-safe
Prym-Chabauty computation together with an Abel--Prym Mordell--Weil sieve.

**Test run.**  The proposed Prym

```text
D: y^2 = -3*x^6 + 24*x^3 - 75
```

has exact rank `1`, trivial rational torsion, and an explicit generator
`(x^2+2*x+4, 5*x+5)`.  This confirms the Chabauty inequality.  The old Coleman
run must not be reused without fixing its absolute-precision and infinity-disk
handling.

### 6. `[8,8]`, order 64

**Why it is plausible.**  Split examples exist, and the reduced second-
halving cover has genuine good-open local points already at `p=7`.

**Strategy.**  Eliminate the auxiliary variables from the reduced square-
subcover or compute the genus-one pair fibers explicitly.  Determine their
ranks and rational points, then exact-test any second-halving lift.

**Test run.**  The reduced equations exactly matched the intrinsic
`T_x in 4J(F_p)` test: 8 solutions at `p=7` and 24 at `p=11`, with no
mismatches.  A fresh height-8 run checked 7,569 pairs, found 24 first-cover
solutions, and no second-cover solution.  Earlier target-specific work is
cold through height 80, so elimination is the appropriate next step.

### 7. `[4,12]`, order 48

**Why it is plausible.**  The full `M(2,12)` surface has explicit algebraic
equations for halving an independent rational 2-class, and the new record
offers a second route through Richelot redistribution of the 2-primary part.

**Strategy.**  On the full `M(2,12)` surface, analyze the degree-eight reduced
`s=m^2` equation on the mod-7 boundary rather than continuing the affine
search.  In parallel, derive the `[4,12]` condition on each of the 15
Richelot codomains over the full `A(2,2,2,12)` parameter space.

**Test run.**  The fixed record's 15 immediate codomains and its complete
depth-three component have no `[4,12]`.  Existing full-surface searches are
negative through height 50 and force any solution to a mod-7 boundary.  The
current LMFDB API also has no `[4,12]` row; this is only a bounded-database
check.

### 8. `[2,2,4,8]`, order 128

**Why it is plausible.**  HLP proves exact arithmetic feasibility on a split
locus, while the project already has a geometrically simple exact `[2,4,8]`
seed.

**Strategy.**  Work on the three genus-one simultaneous pair fibers through a
normalized HLP point, or impose the additional rational 2-direction on the
Richelot graph of a simple `[2,4,8]` seed.  This is preferable to another
direct K3 height search.

**Test run.**  Both HLP examples are recognized by the normalized full cover,
with 32 full witnesses each.  Fresh bounded runs found no nonsplit point:
3,908 direct cover triples, 7,224 full-split tangent pairs, and 84 exact
`[2,2,2,8]` base fibers all gave zero target lifts.

### 9. `[2,2,2,24]`, order 192

**Why it is plausible.**  It is the most direct next-record group: halve the
order-12 generator on the now-proven `A(2,2,2,12)` locus.  The new source
proves that the base locus itself has rational points off the split locus.

**Strategy.**  Construct the explicit 2-division/Kummer cover of the marked
order-12 class over the `A(2,2,2,12)` parameter space.  Use the new record as
a rational base point for local expansions, but search transverse directions
rather than the fiber over that point.  Start with finite-prime images and CRT
lifting before attempting a high-genus global model.

**Test run.**  On the record fiber, the exact order-12 generator is not
divisible by two and the complete rational Richelot component has no
`[2,2,2,24]` vertex.  This closes the fixed fiber but leaves the global finite
cover open.

### 10. `[2,2,4,12]`, order 192

**Why it is plausible.**  This is the other natural order-192 neighbor: keep
the order-12 direction and halve an additional independent rational 2-class.

**Strategy.**  Write the 15 rational `J[2]` classes on the
`A(2,2,2,12)` chart and impose Stoll's squareclass halving criterion for a
second class.  Quotient by the level-2 symmetry before searching.  The
Richelot formulation may lower the cover degree by choosing the desired
kernel first.

**Test run.**  On the record source exactly one of the 15 nonzero `J[2]`
classes is divisible by two,

```text
u = (x+519)*(x+39) = x^2 + 558*x + 20241,
```

and that half is already responsible for the source's existing order-4
direction.  No second half and no `[2,2,4,12]` Richelot vertex occur.

## Near misses just outside the top ten

The next five are `[2,2,16]`, `[4,16]`, `[48]`, `[60]`, and `[35]`.

- `[2,2,16]`: an explicit halving cover exists, but open searches are forced
  into simultaneous 7/11-adic boundary charts.
- `[4,16]`: the known example is bielliptic; the non-bielliptic component is
  locally and computationally colder than `[8,8]`.
- `[48]`: the `A(16)+3` searches produce exact `[16]` survivors but no
  rational 3-class; derive the symbolic 3-contact cover.
- `[60]`: both `M(12)+5` and order-20-plus-3 routes exist, but the easiest
  order-30 halving lane is locally obstructed.
- `[35]`: the simultaneous contact-5/contact-7 equations are structurally
  clear but currently collapse to difficult boundary components.

## Recommended next work

The immediate publication task is to package four theorems/examples together:

1. the order-96 `[2,2,2,12]` record;
2. the exact simple `[2,2,4,4]` example;
3. the exact simple `[2,2,2,8]` example; and
4. the exact simple `[2,4,4]` example.

For discovery work, the best balanced next attack is `[2,24]`: its 2-primary
conditions are abundant and it has no detected global local obstruction.
The highest-leverage deep attack remains the corrected `[6,12]`
Prym-Chabauty/Mordell--Weil sieve.  The highest-upside speculative attack is
the global `[2,2,2,24]` halving cover over `A(2,2,2,12)`.

## Reproducible artifacts created in this audit

```text
code/verify_record_22212_order96.m
code/certify_frontier_existing_banks.m
code/record_22212_richelot_bfs.m

notes/certify_frontier_existing_banks_20260718.md
notes/geometrically_simple_torsion_frontier_2026_07_18.md
notes/record_22212_richelot_bfs_2026_07_18.md
notes/top10_middle_targets_smoke_2026_07_18.md
notes/top10_root_smoke_2026_07_18.md

results/certify_frontier_existing_banks_20260718.log
results/record_22212_richelot_bfs.log
results/record_22212_richelot_bfs_depth3.log
```
