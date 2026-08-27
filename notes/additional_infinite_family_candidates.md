# Additional Infinite Torsion Families Not In The Main Table

This note records a follow-up audit for positive-dimensional genus-2 torsion
families not explicitly listed in the main inventory table.

The convention is the same as in `notes/infinite_families_inventory.md`:
`[n1,n2,...]` means a subgroup isomorphic to
`Z/n1Z x Z/n2Z x ...`, and "torsion" means "torsion contained generically"
unless exactness is stated.

## Locally usable inherited rows

These are not new moduli constructions.  They are immediate subgroups of
families already in the repository, but they are legitimate infinite simple
families under the inventory convention.

| Torsion | Source family | Reason it is valid |
|---|---|---|
| `[3]` | `[30]`, `[12]`, or `[6]` | the cyclic order-30 family is one-parameter, generically simple, and contains a unique order-3 subgroup |
| `[10]` | `[30]`, `[20]`, or `[2,20]` | the cyclic order-30 and order-20 families are generically simple and contain order-10 subgroups |
| `[15]` | `[30]` | the cyclic order-30 family contains an order-15 subgroup |
| `[2,6]` | `[2,12]` | `Z/2 x Z/12` contains `Z/2 x Z/6`; the `[2,12]` family has simple specializations |
| `[2,10]` | `[2,20]` or `[2,2,2,10]` | both source families contain `Z/2 x Z/10`; the `[2,20]` loci and the Elkies `[2,2,2,10]` source have simple evidence |
| `[2,2,10]` | `[2,2,2,10]` | `Z/2 x Z/2 x Z/2 x Z/10` contains `Z/2 x Z/2 x Z/10`; simple evidence is inherited from the source family |
| `[4,4]` | `[4,8]` | `Z/4 x Z/8` contains `Z/4 x Z/4`; the tangent-cover `[4,8]` family has many simple specializations |
| `[16]` | reconstructed Elkies `[32]` | the marked class `D` has order `32` generically, so `2D` gives a cyclic order-`16` subgroup on the same one-parameter simple component |

The main local files are:

```text
notes/contact5_contact6_order30_family.md
notes/contact5_order40_family.md
notes/m12_simple_route.md
notes/m18_m14_halving.md
notes/elkies22210_richelot.md
notes/elkies32_reconstruction.md
```

These rows are useful for completeness, but they do not create new high-torsion
search routes.  For example, the `[15]` row is just the odd part of the existing
`[30]` family, not a separate construction that might combine independently
with more 2-power torsion.

## Reconstructed local row: `[32]`

Elkies' page

```text
https://people.math.harvard.edu/~elkies/g2_tors.html
```

announces a one-parameter family of absolutely simple genus-2 Jacobians with a
rational point of order `32`.  The page explicitly does not print the
parameterized formula because the coefficients are complicated rational
functions of the parameter.

The page gives one member:

```text
y^2 = (15*x - 1)*(1056*x^4 + 156183*x^3 + 26297*x^2 + 649*x - 121).
```

The repository originally only recorded a halving analysis of the visible
order-32 points on this printed member:

```text
notes/elkies32_halving.md
code/elkies32_halving_conditions.m
data/elkies32_halving_conditions.txt
```

That gap is now closed for local use.  The reconstructed contact model is

```text
y^2 = (a*x^3 + b*x^2 + c*x + 1)^2 - a^2*x^5*(x+1),
```

with `(z,r)` on an explicit genus-0 plane curve and
`p=4*Pnum(z,r)/Pden(z,r)`.  The class `D=(0,1)-infinity` satisfies
`16D=(r,0)-infinity`, so generically `D` has order `32`.

Main local files:

```text
notes/elkies32_reconstruction.md
code/elkies32_reconstruct_conditions.py
data/elkies32_reconstruct_conditions.txt
code/elkies32_component_genus.m
data/elkies32_component_genus.txt
code/elkies32_simple_certificate.m
data/elkies32_simple_certificate.txt
```

The printed specialization has a local `Q`-simplicity certificate at `p=7`:

```text
L_7(T) = 49*T^4 + 7*T^3 + 6*T^2 + T + 1
```

which is irreducible.  A local Lombardo/geometric-simplicity audit is still a
useful cleanup step, but the family is no longer literature-only.

## Subgroup audit after `[32]`

I normalized subgroup types from the established simple-backed source groups

```text
[2,12], [2,20], [2,2,2,10], [4,8], [2,4,8],
[30], [32], [22], [23], [20], [18], [14].
```

After converting Sage's elementary-factor output to invariant-factor notation,
all inherited subgroup families were already in the inventory except `[16]`.
For example, apparent rows such as `[2,3]`, `[4,5]`, or `[2,2,5]` normalize to
already-listed `[6]`, `[20]`, and `[2,10]`.

## Large but not simple

Several large infinite torsion families exist in the split/gluing world, but
they are not useful for the current simple-Jacobian target.

- Howe's order-48 family is parametrized by a rank-2 elliptic curve, but the
  construction is by gluing elliptic curves along torsion, so the resulting
  abelian surface is isogenous over `Qbar` to a product of elliptic curves.
- Howe-Leprevost-Poonen construct large split families, including genus-2
  Jacobians with many rational torsion points, but the Jacobians split
  completely.

These are worth remembering as non-simple comparison examples, not as
candidates for the simple-family inventory.

## Isolated large simple examples

Elkies' page gives isolated simple curves with rational torsion orders `34`,
`39`, and `40`, and Howe gives isolated large examples including order `70`.
These are valuable benchmarks, but they do not add infinite families unless a
positive-dimensional parametrization is available.

## Best follow-up

The subgroup audit adds only `[16]`; it does not create a new route to larger
torsion.  The next constructive cleanup for the `[32]` branch is to find a
smaller rational parametrization of the reconstructed genus-0 curve.  The next
large-torsion move should be a new native construction, not another inherited
subgroup row.
