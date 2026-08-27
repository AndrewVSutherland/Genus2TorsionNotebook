# Census audit: infinite families, endomorphisms, and 2P−K classes

Date: 2026-08-01

## Objective

Add to the paper source the realization census of
`paper/torsion_realizations.tex` (all 72 invariant-factor types with a
documented geometrically simple exact realization over `Q`), extended by three
per-group columns:

- **F** — is an infinite family of examples known?
- **E** — do *all* known geometrically simple examples have nontrivial
  geometric endomorphisms?
- **Θ** — do the known examples carry a torsion divisor of the form
  `2P − K` (`K` canonical, `P ∈ C(Q)`, doubled points at infinity allowed)?

The table now sits in the paper draft as Section "The realization census and
three structural columns" (`tab:census`), with the audit conventions stated
there.  This note records data provenance, the new computations, and flags.

## Column F (infinite families)

Derived from `notes/infinite_families_inventory.md` (main + filter +
additional-rows tables), `paper/torsion_sources.tex`, the paper's families table
(`tab:families`), and the **new fiber certificates** below.  Grading:

- `∞` — recorded positive-dimensional family, infinitely many rational
  fibers, generic exact torsion equal to the group (marked subgroup + one
  exact fiber), generic geometric simplicity (one strict fiber):
  `[5] [6] [7] [9] [11] [12] [14] [15] [18] [20] [21] [22] [23] [30] [32]
  [2,6] [2,12] [2,20] [3,3] [4,8] [2,2,2,10]`.
- `(∞)` — positive-dimensional construction recorded, but generic
  exactness, generic *geometric* simplicity, or infinitude of rational
  fibers not established:
  `[2] [3] [4] [8] [10] [13] [16] [17] [19] [26] [28] [2,2] [2,4] [2,8]
  [2,10] [3,9] [4,4] [6,6] [2,2,2] [2,2,4] [2,2,8] [2,2,10] [2,4,4] [2,4,8]
  [2,2,2,2] [2,2,2,4] [2,2,2,6] [2,2,2,8] [2,2,2,12] [2,2,4,4]`.
- `—` — no positive-dimensional construction recorded:
  `[24] [25] [27] [29] [31] [33] [34] [36] [39] [40] [2,14] [2,16] [2,18]
  [2,22] [2,26] [2,28] [3,6] [2,2,6] [2,2,12] [2,2,14] [2,2,20]`.

### New fiber certificates (this audit)

`code/claude_census_family_fibers.m` → `data/claude_census_family_fibers.txt`.
Strict = root-power test `[Q(π^n):Q]=4` for `n ≤ 12` at the printed prime.
Control: the `[6,6]` curve reproduces its recorded strict prime `37`.

| family | sample fibers | exact torsion | strict primes |
|---|---|---|---|
| contact-7 chart `(a,b)` | (1,1),(2,1),(1,−1),(3,2) | `[7]` (all) | 5,5,7,17 |
| plain `M(12)` chart `(a,r)` | (1,3),(2,3),(1,−2),(3,2) | `[6]` (all) | 13,11,13,11 |
| `M(12)` z-subchart `a=(1−z²)/(4(r+1))` | (r,z)=(3,2),(2,3),(5,2),(3,4) | `[12]` (all) | 13,11,17,11 |
| `[2,12]` line `a=(1−r)/4` | r=3,5,−3 | `[2,12]` (all) | 17,11,7 |
| contact-6 chart `(a,b)` | (1,1),(2,−1),(3,1) | `[2,6]` (all) | 0,29,29 |

Consequences: `[6]`, `[7]`, `[12]`, `[2,6]`, `[2,12]` upgraded to `∞`.  The
`[2,12]` family previously had only `Q`-simple (irreducible-Frobenius)
certificates (`code/m12_z12x2_family.m`); the strict primes above certify
geometric simplicity of fibers, hence of the generic fiber.

Two clarifications the inventory should absorb:

- **Contact-6 is a `[2,6]` chart, not a `[6]` chart**: `h6²−(x−1)⁶` always
  factors (difference of squares) into `quadratic · x · quadratic`, so a
  rational 2-class is built in; all sampled fibers have exact `[2,6]`.
  Exact `[6]` fibers come instead from the plain `M(12)` chart.
- **Plain `M(12)` carries order 6, not 12**: the marked structure is the
  doubled class `2P−K` of order 6; rational order-12 classes appear on the
  rational-Weierstrass subchart `1−4a(r+1)=z²` (consistent with
  `notes/m12_simple_route.md`, step 2).  The inventory line "the marked
  point gives a Jacobian class of order 12" holds only on that subchart.

## Column E (endomorphisms)

From the endomorphism refinement of `paper/torsion_realizations.tex`:
"yes" exactly for `[2,22]` (extended-database geometrically simple witnesses
all RM) and `[31]` (single known witness, RM by `Z[√2]`).  Every other group
has a documented witness with `End(J_Qbar) = Z` (for `[2,2,14]` the ten
generic witnesses are new to this project).

**Flag**: `torsion_realizations.tex` cites
`data/claude_endz_certificates.txt` for the non-database witnesses, but that
file does not exist in the checkout.  The E column does not depend on it
(it uses the audit statements themselves), but the reference should be fixed
or the file restored before the paper is finalized.

## Column Θ (doubled-point / 2P−K classes)

`code/claude_census_2pk_audit.m` → `data/claude_census_2pk_audit.txt`.

Criterion (isomorphism-invariant; unique effective `D ∈ |K+g|` for `g ≠ 0`):
in the reduced Mumford representation `⟨a,b,d⟩` of a torsion class `g`,
`g = [2P−K]` iff `d = 2` and either `deg a = 2` with `disc a = 0` (finite
doubled point) or `deg a ≤ 0` (doubled point at infinity).  Validated on
Flynn's order-11 family at `t=1`, where `∞₊−∞₋` is represented as
`(1, −x³−x², 2)` and is caught by the `deg a = 0` branch.

All 72 displayed witnesses plus the other two known `[2,2,2,12]` fibers were
audited.  **Cross-check: the recomputed exact torsion matches the census
group in every row.**  Θ orders per row are in the data file; summary:

- `•` (doubled class of maximal invariant-factor order):
  `[7] [9] [11] [12] [13] [15] [17] [19] [21] [23] [24] [25] [27] [29]
  [33] [39] [3,3]`.
- `◦` (doubled classes of strictly smaller order only):
  `[8] [10] [14] [16] [18] [20] [22] [26] [28] [30] [32] [34] [36] [40]
  [2,6] [2,8] [2,10] [2,12] [2,14] [2,16] [2,18] [2,20] [2,22] [2,26]
  [2,28] [4,8] [6,6] [2,2,6] [2,2,8] [2,2,10] [2,2,12] [2,2,14] [2,2,20]
  [2,4,8] [2,2,2,6] [2,2,2,10] [2,2,2,12]`.
- `—` (no doubled class at all):
  `[2] [3] [4] [5] [6] [31] [2,2] [2,4] [3,6] [3,9] [4,4] [2,2,2] [2,2,4]
  [2,4,4] [2,2,2,2] [2,2,2,4] [2,2,2,8] [2,2,4,4]`.

Structural notes:

- An order-2 class is never `2P−K` (it would force `2P ∼ W₁+W₂` with
  distinct Weierstrass points), so exponent-2 groups are `—` necessarily.
- All three `[2,2,2,12]` fibers carry exactly the order-3 doubled class
  (`D₃ = 2∞₊ − K = ∞₊−∞₋` of the paper's symmetric construction) and nothing
  else of doubled type — no order-4/6/12 doubled class.
- The `[5]` witness (277.a.2) has **no** doubled class, even though
  `[2P−K] = 2[P−W]` would inherit odd order from any point-supported
  order-5 class; so its 5-torsion is not point-supported at all.
- **Flag**: the displayed `[2,2,2,8]` witness
  `y² = x(x+1)(x+55²)(x+99²)(x+125²)` (bank `tor2228.txt` row 3, see
  `notes/certify_frontier_existing_banks_20260718.md`) has **no doubled
  class**, hence does not lie on the strict `M(2,2,2,8)` K3 for any
  normalization — consistent with `(±1,±55,±99,±125)` failing
  `s₂² = 4abcd` under every sign choice.  Its order-8 point lives on the
  branch-square cover `A(2,2,2,8)`.  The `tab:families` row
  "`[2,2,2,8]` K3 … exact strict fiber" should be re-examined: if the only
  exact strict fiber meant is this curve, the K3 attribution is wrong.

Scope caveat (stated in the table caption): a `—`/`◦` entry *refutes*
"every known example has a maximal-order `2P−K` class" (the witness is a
counterexample); a `•` entry verifies the witness only — for groups with
many database examples the universal statement remains unaudited.

## Reproduction

```sh
cd ~/torsion_jac
magma -b code/claude_census_2pk_audit.m      > data/claude_census_2pk_audit.txt
magma -b code/claude_census_family_fibers.m  > data/claude_census_family_fibers.txt
```

Both runs complete in well under a minute (Magma V2.29-4, memory cap set in
the scripts).

## Files touched

- the paper source — new census section +
  `longtable`, intro count 70→72, reproducibility bullets, 16 new bibitems
  (mirrored to the paper-repo checkout).
- `code/claude_census_2pk_audit.m`, `data/claude_census_2pk_audit.txt`
- `code/claude_census_family_fibers.m`, `data/claude_census_family_fibers.txt`
- this note.
