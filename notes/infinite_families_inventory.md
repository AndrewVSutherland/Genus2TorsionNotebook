# Infinite torsion families with local equations/code

This is a workspace inventory, not a literature classification.  It lists the
positive-dimensional torsion families for which this repository currently has
usable equations or code.

For the low-level 2-primary spaces there are several natural components
(rational Weierstrass support, Richelot/quadratic-pair support, and marked
point covers).  The entries below record the components whose normal forms are
explicit in the local paper/notes; they are not meant to be a complete
component classification.

Notation: `[n1,n2,...]` means
`Z/n1Z x Z/n2Z x ...`.  The column "torsion" means "torsion contained
generically" unless the row explicitly says "generic exact".  Some rows are
fully rationally parametrized; others are algebraic covers with local equations
and verified rational points but no compact parametrization recorded in the
notes.

## Summary table

| Torsion | Family/source | Status | Main local files |
|---|---|---|---|
| `[2]` | rational-Weierstrass chart | parametrized baseline | `paper/main.tex`, this note |
| `[2,2]` | Richelot/quadratic-pair chart | parametrized baseline | `paper/main.tex`, this note |
| `[2,2,2]` | four-rational-Weierstrass chart | parametrized baseline | `paper/main.tex` |
| `[2,2,2,2]` | full-split genus-2 chart | parametrized baseline | `notes/seven_torsion_hybrids.md`, `code/m2222_plus7_search.m` |
| `[4]` | `A_1(4)` normal form | weighted-projective parametrization | `paper/NotesAndTodo.tex` |
| `[2,4]` | `A_1(4,2)` normal form | three-parameter chart | `paper/NotesAndTodo.tex` |
| `[8]` | strict `M_1(8)` and rational-Weierstrass cover | rational surface; two-parameter rational-Weierstrass cover | `paper/NotesAndTodo.tex` |
| `[2,8]` | `M_1(8,2)` quadratic-factor chart | weighted-projective chart | `paper/NotesAndTodo.tex` |
| `[2,2,4]` | rational-Weierstrass component of `A(2,2,4)` | four-parameter chart | `paper/main.tex` |
| `[2,2,2,4]` | square branch-point chart | parametrized baseline | `notes/m2224_plus3.md`, `code/m2224_plus3_search.m` |
| `[2,2,2,6]` | split `M(2,2,2,6)` family | parametrized | `notes/m2226_summary.txt`, `notes/m2226_order6_doubling.md` |
| `[2,2,2,8]` | `M(2,2,2,8)` K3 surface | algebraic K3 family | `notes/m2228_three_torsion.md`, `code/m2228_three_torsion_surface_sieve.m` |
| `[2,2,8]` | odd `M_1(8,2,2)` family | parametrized | `notes/m3222_plus3.md`, `notes/m3222_halving_boundary.md` |
| `[2,2,4,4]` | `A(2,2,4,4)` squareclass threefold | algebraic threefold | `notes/a2244_local_obstructions.md`, `code/m2248_sieve.m` |
| `[2,4,4]` | `M(2,4,4)` elliptic fiber-product family | parametrized by fiber-product points | `notes/m244_to_248_route.md`, `code/m244_to_248_sample_search.m` |
| `[4,8]` | `M_1(8,4)` tangent-cover family | algebraic cover, many simple examples | `notes/m18_m14_halving.md`, `code/m18_m14_tangent_search.m` |
| `[2,4,8]` | one-split subfamily inside `[4,8]` | algebraic subcover with verified simple examples | `notes/m18_m14_halving.md`, `code/m18_m14_one_split_search.m` |
| `[5]` | quintic-contact family | two-parameter | `notes/m10_quintic_contact5.md`, `code/m10_quintic_contact5_search.m` |
| `[6]` | `M(6)` model and contact-6 family | rational/infinite; contact chart also available | `paper/main.tex`, `notes/contact6_m36.md`, `code/contact6_m36_search.m` |
| `[6,6]` | independent cubic-contact class on the contact-6 `[1,2,2]` core | positive-dimensional algebraic cover; exact geometrically simple example | `notes/contact6_m36_66_example.md`, `code/contact6_m36_66_package.sage` |
| `[7]` | contact-7 family | two-parameter | `notes/contact7_family.md`, `code/contact7_family_search.m` |
| `[9]` | contact-9 family | one-parameter | `notes/contact9_family.md`, `code/contact9_family_search.m` |
| `[11]` | Flynn and Daowsud-Schmidt divisor-at-infinity families | one-parameter, exact and geometrically simple samples | `notes/internet_infinite_families_sweep.md`, `code/order11_family_check.m`, `code/order11_family_geom_check.sage` |
| `[12]` | `M(12)` and split `M(2,12)` charts | two-parameter charts | `paper/main.tex`, `notes/m12_simple_route.md`, `notes/m212_three_torsion.md` |
| `[2,12]` | one-parameter extra-root line in `M(12)` | parametrized, simple examples | `notes/m12_simple_route.md`, `code/m12_z12x2_family.m` |
| `[14]` | rational-root subfamily of contact-7 | parametrized subfamily | `notes/contact7_family.md`, `code/contact7_root_even_search.m` |
| `[18]` | rational-root subfamily of contact-9 | parametrized subfamily | `notes/contact9_family.md`, `code/contact9_family_search.m` |
| `[20]` | contact-5 plus explicit 4-torsion | one-parameter | `notes/contact5_order40_family.md`, `code/contact5_order40_family_search.m` |
| `[2,20]` | extra-2 loci in the contact-5/order-20 family | two one-parameter loci | `notes/contact5_order40_family.md`, `code/contact5_extra2_param_large_search.m` |
| `[21]` | Leprévost divisor-contact family, normalized onto contact-7 plus cubic contact | one-parameter, generic exact; geometrically simple exact sample | `notes/z21_leprevost_and_hlp_2026_07_18.md`, `code/z21_leprevost_family_verify.m`, `code/contact7_plus3_leprevost_bridge.m`, `code/leprevost_z21_geom_check.sage` |
| `[22]` | rational-branch subfamilies of order-11 infinity-torsion families | explicit one-parameter families, exact and geometrically simple samples | `notes/order22_from_order11.md`, `code/order22_from_order11_check.m` |
| `[23]` | Kuru-Sadek quadratic-order family, genus-2 specialization | one-parameter, exact and geometrically simple sample | `notes/internet_infinite_families_sweep.md`, `code/order23_kuru_sadek_check.m`, `code/order23_kuru_sadek_geom_check.sage` |
| `[30]` | simultaneous contact-5/contact-6 family | one-parameter, generic simple | `notes/contact5_contact6_order30_family.md`, `code/contact5_contact6_order30_family.m` |
| `[32]` | reconstructed Elkies contact family | implicit genus-0 one-parameter family; simple printed specialization at `p=7` | `notes/elkies32_reconstruction.md`, `code/elkies32_reconstruct_conditions.py`, `code/elkies32_simple_certificate.m` |
| `[2,2,2,10]` | Elkies Clebsch-Klein full-level-2 plus 5 family | algebraic family, generic exact in searches | `notes/elkies22210_richelot.md`, `code/elkies22210_richelot_sweep.m` |

## Additional Rows From Follow-Up Audit

The main table intentionally emphasizes the largest/native constructions rather
than every subgroup inherited from them.  A later audit found the following
additional positive-dimensional rows that were not explicitly listed above.

| Torsion | Status | Source |
|---|---|---|
| `[3]` | inherited simple family | subgroup of `[30]`, `[12]`, or `[6]` |
| `[10]` | inherited simple family | subgroup of `[30]`, `[20]`, or `[2,20]` |
| `[15]` | inherited simple family | subgroup of the cyclic `[30]` family |
| `[2,6]` | inherited simple family | subgroup of `[2,12]` |
| `[2,10]` | inherited simple family | subgroup of `[2,20]` and `[2,2,2,10]` |
| `[2,2,10]` | inherited simple family | subgroup of `[2,2,2,10]` |
| `[4,4]` | inherited simple family | subgroup of `[4,8]` |
| `[16]` | inherited simple family | subgroup of the reconstructed Elkies `[32]` family; if `D` has order `32`, then `2D` has order `16` |

Details and caveats are recorded in
`notes/additional_infinite_family_candidates.md`.

## Strict `M`-Space Status

The table above is organized by visible torsion group.  The local paper also
uses stricter `M` notation for loci where the marked torsion class comes from a
point on the curve.  The following `M`-spaces are explicitly infinite in the
local sources.

| Space | Infinite status | Normal form / condition |
|---|---|---|
| `M_1(2)` | all of `M_2`, dimension 3 | no extra condition |
| `M_1(4)` | empty | no `P-W` of order 4 |
| `M_1(8)` | rational surface | `q=a*x^2+b*x+c`, `C: y^2=q*(x^4+q)` |
| `tilde M_1(8)` | rational surface with rational Weierstrass support | `q=B*(x-B)*(x-C)`, `C: y^2=q*(x^4+q)` |
| `M_1(8,2)` | open subset of `P(1,2,2)` | product of three quadratics below |
| `M_1(8,2^w)` | open subset of `P(1,1,2)` | specialization `v=(u^2-s^2)/4` in `M_1(8,2)` |
| `M_1(8,2,2)` | open subset of `P^2` | formula below |
| `M_1(8,2,2,2)` / `M(2,2,2,8)` | K3 surface, `Q`-points Zariski dense | `s2^2=4*s4` in square branch coordinates |
| `M_1(8,4)` | infinite tangent-cover family | see the `[4,8]` tangent-cover section |
| `M(6)` | rational, birational to `P(1,2,3)` | `y^2=x*(x^2+h1*x+h2)*(x^3+h1*x^2+h2*x+2*h3)` |
| `M(12)` | rational, birational to `A^2_{a,r}` | `y^2+(x-r)*(T+1)*y=a*x^2*T*(T+1)`, `T=a*x^2-x+r` |
| `M(2,12)` | two-parameter split chart | `1-4*a*(r+1)=z^2`, so `a=(1-z^2)/(4*(r+1))` |
| `M(2,2,2,4)` | rational `P^3` chart | full-split square branch model `x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2)` |
| `M(2,2,2,6)` | rational surface, birational to `P^2` | split `M(6)` model in `s,m,n` |
| `M(2,2,2,8)` | K3 surface, `Q`-points Zariski dense | same as `M_1(8,2,2,2)` |

The nearby notation `M(2,2,4)` is a trap in this convention: the strict `M`
condition would try to realize a nonzero order-2 class as `2P-K`, and the same
argument that kills `M_1(4)` rules this out.  The useful nonempty moduli space
for torsion `[2,2,4]` is `A(2,2,4)`, whose rational-Weierstrass component is
listed in the equations below.

Key formulas not repeated elsewhere:

```text
M_1(8):
  q = a*x^2 + b*x + c
  C: y^2 = q*(x^4 + q)
  marked point: (0,c)

M_1(8,2):
  C: y^2 =
    (x^2 - u*x + v)
    *(x^2 + u*x + t)
    *((v + t - u^2)*x^2 + u*(v-t)*x + v*t)
  marked point: (0, v*t)

M_1(8,2,2):
  C: y^2 = (x-1)*(x-u)*(x-v)*(x+u+v+1)*q(x)
  q = -(u^2+u*v+v^2+u+v+1)*x^2
      +(v+1)*(u+1)*(u+v)*x
      -u*v*(u+v+1)
  marked point: (0, u*v*(u+v+1))
```

## Simple-Jacobian Filter

For the current search, the useful notion is not merely "positive-dimensional
with prescribed torsion", but "positive-dimensional and not trapped in the
product locus".  The repository uses two levels of local evidence:

- an irreducible degree-4 Frobenius/L-polynomial at one good prime, which is
  the standard fast certificate used by the scripts for a simple Jacobian over
  `Q`;
- the stronger Lombardo/12th-power transform check, when recorded, which rules
  out geometric splitting over `Qbar`.

One certified simple specialization on an irreducible positive-dimensional
component shows that the generic member is outside the corresponding product
locus; the nonsimple members then lie in a proper closed exceptional subset.
Rows below marked only by irreducible Frobenius should still receive the
12th-power audit before being advertised as geometrically simple.

The following rows are the infinite families in this inventory that are
supported by local simple-specialization evidence, either directly or by
containing a certified simple higher-torsion subfamily.

| Torsion contained | Simple-family evidence | Strength / caveat |
|---|---|---|
| `[5]` | the quintic-contact surface contains the simple order-20 subfamily | simple generic on the contact-5 component; exact `[5]` occurs off torsion-jump divisors |
| `[6]` | the contact-6 chart contains the simple `[30]` family, and the core-cover search has a simple `[6,6]` point | simple generic on the contact-6 component |
| `[6,6]` | the core point `(a,b)=(133/39,-7/13)` has exact torsion `[6,6]`, irreducible 12th-power Frobenius transform at `p=37`, and Sage geometric endomorphism ring `ZZ` | geometrically simple exact example on a positive-dimensional algebraic cover; no rational `Q(t)` subfamily is known |
| `[7]` | contact-7 contains the known simple `[28]` point and the rational-root `[14]` curve | simple generic on the contact-7 component |
| `[9]` | contact-9 samples have simple certificates; the rational-root subfamily has a simple `[18]` point | simple generic for `[9]` and for the `[18]` subfamily |
| `[11]` | Flynn and Daowsud-Schmidt one-parameter families have exact `[11]` specializations with Sage geometric endomorphism ring `ZZ` | genuine one-parameter simple odd-prime families |
| `[12]` | `code/m12_simple_search.m` finds irreducible Frobenius certificates immediately | simple generic on `M(12)` |
| `[2,12]` | the one-parameter extra-root line has a verified simple sample, with certificate at `p=11` in `code/m12_z12x2_family.m` | genuine one-parameter simple family |
| `[14]` | the contact-7 rational-root subfamily contains a simple `[28]` specialization | generic member of the `[14]` subfamily is simple |
| `[18]` | `notes/contact9_family.md` records a simple exact `[18]` hit and many exact `[18]` root-subfamily hits | generic member of the `[18]` subfamily is simple |
| `[20]` | the order-20 contact-5 subfamily contains the simple order-40 specialization `t=-1/3` | generic order-20 family is simple |
| `[2,20]` | parametrized extra-2 loci have many exact `[2,20]` hits, and the locus contains a Lombardo-certified geometrically simple `[2,2,20]` specialization | genuine one-parameter simple families on the `1+3` and `2+2` loci |
| `[21]` | Leprévost's `t=1` fiber has exact torsion `[21]`; Sage gives geometric endomorphism ring `ZZ` | genuine one-parameter, generic-exact odd-torsion family; the simplicity certificate is local, not claimed in the 1991 source |
| `[22]` | forcing a rational branch point in the Flynn/Daowsud-Schmidt order-11 families gives exact `[22]` samples with Sage geometric endomorphism ring `ZZ` | explicit one-parameter simple families; likely related parameterizations |
| `[23]` | Kuru-Sadek genus-2 specialization at `t=2` has exact torsion `[23]`, irreducible Frobenius at `p=13`, and Sage geometric endomorphism ring `ZZ` | best new internet find; one-parameter simple odd-prime family |
| `[30]` | simultaneous contact-5/contact-6 family has many exact `[30]` samples with simple certificates | one of the best large odd/even simple families here |
| `[32]` | reconstructed Elkies component contains the printed member; `code/elkies32_simple_certificate.m` gives irreducible `L_p` at `p=7` | local `Q`-simple certificate for the component; Elkies' page states absolute simplicity, but a local Lombardo audit is still useful |
| `[4,8]` | tangent-cover family has many exact `[4,8]` examples with simple certificates | best 2-primary simple family currently recorded |
| `[2,4,8]` | one-split subcover has a verified exact `[2,4,8]` hit with irreducible 12th-power transform at `p=47` | positive equation-defined family; compact rational parametrization not recorded |
| `[2,2,2,10]` | Elkies Clebsch-Klein source samples have exact `[2,2,2,10]` and simple certificates throughout `data/elkies22210_source_h100.txt` | strong full-level-2 plus odd torsion family |

Several lower-torsion ambient spaces are also simple in this sense by
inheritance.  For example, the simple `[4,8]`, `[2,4,8]`, and
`[2,2,2,10]` families give simple specializations on the ambient spaces with
contained torsion `[2]`, `[2,2]`, `[4]`, `[2,4]`, `[8]`, `[2,8]`,
`[2,2,4]`, and full rational `2`-torsion.  These are not the bottleneck for
large torsion, so the table above emphasizes the largest native constructions
rather than every inherited subgroup.

I would not yet count the following as simple infinite families, despite useful
positive examples or equations:

- `[2,22]`: the natural extra-2 test on the order-22 family is recorded in
  `notes/order222_from_order11.md`; the quadratic component has only boundary
  rational points and the remaining cubic-root condition has genus 3.

- `[28]`, `[40]`, `[2,2,20]`, and `[6,6]`: simple examples exist, but the
  repository does not currently prove a positive-dimensional family with that
  exact torsion.
- `[2,2,2,6]`, `[2,2,2,8]`, `[2,2,8]`, `[2,2,4,4]`, and `[2,4,4]`: there are
  useful parametrizations or algebraic spaces, but this note should not mark
  them as simple-only search bases until a clean local simple certificate is
  attached to the relevant component.
- `[8,8]`, `[4,16]`, `[35]`, `[60]`, `[72]`, `[80]`, and `[120]`: current
  material is searches, local obstructions, or isolated/boundary evidence, not
  infinite simple families.

## Equations and construction details

### Rational 2-torsion: `[2]`

The simplest rational-Weierstrass chart is

```text
C: y^2 = x*(x^4 + a*x^3 + b*x^2 + c*x + d),
```

with `d != 0` and squarefree right-hand side.  The class
`(0,0)-infinity` has order `2`.

This is only the rational-Weierstrass component of rational `2`-torsion; a
general rational 2-torsion class may be represented by a rational quadratic
factor rather than two rational branch points.

### Richelot/quadratic-pair chart: `[2,2]`

A useful computable model for a rational `(2,2)` subgroup is a sextic split
into three rational quadratic factors:

```text
C: y^2 = q1(x)*q2(x)*q3(x),
```

with the `q_i` pairwise coprime.  The divisor classes cut out by the three
quadratic factors are the three nonzero points of a rational `(2,2)` subgroup.

A rational-pair affine subchart is

```text
C: y^2 = x*(x^2 + a*x + b)*(x^2 + c*x + d),
```

where the pair `{0,infinity}` is one of the three quadratic-pair blocks.

### Four rational Weierstrass points: `[2,2,2]`

The paper records:

```text
(2,2,2) subset J(Q)  iff  C has 4 rational Weierstrass points.
```

A normal form is

```text
C: y^2 = x*(x+u)*(x+v)*(x^2 + a*x + b).
```

Here the rational branch points are `infinity`, `0`, `-u`, and `-v`; they give
torsion containing `[2,2,2]`.

Main source: `paper/main.tex`, section `(2,2,2) and beyond`.

### One point of order 4: `[4]`

The `A_1(4)` normal form in `paper/NotesAndTodo.tex` is

```text
Q = x^2 + d,
q = a*x^2 + b*x + c,
C: y^2 = q*(Q^2 + q).
```

The point of order `4` is represented by the Mumford divisor

```text
D = (Q, 0),
```

with `2D` represented by `(q,0)`.  The parameter space is an open subset of
the weighted projective space `P(2,3,4,2)` in `(a,b,c,d)`.

### One point of order 4 plus independent 2-torsion: `[2,4]`

The `A_1(4,2)` chart is birational to an open subset of `A^3_{t,u,v}`.  Put

```text
beta = t*(v-u),
Q = x^2 + beta,
q = (u+v-2*beta-t^2)*x^2 + beta*x + (u*v-beta^2).
```

Then

```text
C: y^2 = q*(Q^2 + q)
```

has the marked order-`4` class `(Q,0)` and an independent rational `2`-torsion
class, giving torsion containing `[2,4]`.

Main source: `paper/NotesAndTodo.tex`, theorem `(4,2) parameterization`.

### Rational-Weierstrass cover of `M_1(8)`: `[8]`

The rational-Weierstrass double cover of `M_1(8)` is recorded as

```text
q = B*(x-B)*(x-C),
C: y^2 = q*(x^4 + q).
```

The rational point

```text
P = (0, B^2*C)
```

has order `8` relative to a rational Weierstrass point in the support of `q`.
Thus this gives a two-parameter source of rational `8`-torsion.

Main source: `paper/NotesAndTodo.tex`, corollary on `tilde M_1(8)`.

### `M_1(8,2)`: `[2,8]`

The local model is

```text
C: y^2 =
  (x^2 - u*x + v)
  *(x^2 + u*x + t)
  *((v + t - u^2)*x^2 + u*(v-t)*x + v*t).
```

The point of order `8` is recorded as `(0, v*t)`, and the extra quadratic
factorization supplies an independent rational `2`-torsion class.  Hence the
torsion contains `[2,8]` on the smooth open chart.

Main source: `paper/NotesAndTodo.tex`, theorem `M_1(8,2)`.

### Rational-Weierstrass component of `A(2,2,4)`: `[2,2,4]`

The normal form from `paper/main.tex` is

```text
C: y^2 = x*(x+u^2)*(x+v^2)*(x^2 + (t^2+2*s)*x + s^2).
```

It comes from halving the rational `2`-torsion class `(0,0)-infinity` in the
four-rational-Weierstrass model.  Thus the smooth open chart has torsion
containing `[2,2,4]`.

This is one component: the paper notes that `A(2,2,4)` has another component
where the order-`4` point doubles to the rational `2`-torsion divisor coming
from the non-rational Weierstrass pair.

Strict-notation warning: in the paper's `M` convention, an `M(2,2,4)`
condition would ask an order-`2` class to be represented as `2P-K`.  The same
argument used there to show `M_1(4)` is empty rules this out for a nonzero
2-torsion class: it would force `2P` to be linearly equivalent to a sum of two
distinct Weierstrass points.  Thus the useful nonempty moduli space for torsion
`[2,2,4]` is `A(2,2,4)`, not strict `M(2,2,4)`.

### Full rational 2-torsion: `[2,2,2,2]`

Use the normalized full-split model

```text
C: y^2 = x*(x-1)*(x-a)*(x-b)*(x-c),
```

with distinct branch points away from `0,1,infinity`.  This gives full
rational `J[2]`, i.e. torsion containing `[2,2,2,2]`.

Local code mainly uses this as a base for adding `7`-torsion:
`code/m2222_plus7_search.m` and `code/m2222_plus7_boundary_search.m`.

### Square branch-point family: `[2,2,2,4]`

The chart used in the `M(2,2,2,4)+3` computations is

```text
C: y^2 = x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2).
```

On the nonsingular open chart, full rational `2`-torsion is visible and the
square branch-point structure gives a half of one rational `2`-torsion class.
Thus the visible torsion contains `[2,2,2,4]`.

Main files: `notes/m2224_plus3.md`, `code/m2224_plus3_search.m`.

### Split `M(2,2,2,6)`: `[2,2,2,6]`

The even-degree model recorded in `notes/m2226_order6_doubling.md` is

```text
y^2 = x*(x + 2*s^2 - s*n)
        *(x + 2*s^2 + s*m - 2*s*n - m*n)
        *(x + 2*s^2 + s*m - s*n - m*n)
        *(2*x - m*n)
        *(2*x + 4*s^2 - 4*s*n - m*n).
```

After moving `x=0` to infinity, the point `P=(0,2)` gives
`g=P-infinity` of order `6`, while the split branch points give full rational
`2`-torsion.  The generic visible torsion is `[2,2,2,6]`.

Main files: `notes/m2226_summary.txt`,
`code/m2226_order6_doubling.m`, `code/m2226_order6_halving_direct_search.m`.

### `M(2,2,2,8)` K3 surface: `[2,2,2,8]`

The model used throughout the `M(2,2,2,8)+3` work is

```text
C: y^2 = x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2),
```

with the K3 condition

```text
s2(a,b,c,d)^2 = 4*a*b*c*d.
```

This is the full-rational-2 chart with a further halving, giving torsion
containing `[2,2,2,8]` on the nonsingular open surface.

Main files: `notes/m2228_three_torsion.md`,
`notes/m2228_three_torsion_summary.txt`,
`code/m2228_three_torsion_surface_sieve.m`.

### Odd `M_1(8,2,2)`: `[2,2,8]`

The odd two-parameter family is

```text
C_{u,v}: y^2 = f_{u,v}(X),
Q = (-1, u*v*(u+v+1)),
```

where

```text
f_{u,v}(X) =
((1-u)*X+1)*((1-v)*X+1)*((u+v+2)*X+1)
*(-X^2 + (u^2*v-u^2+u*v^2-u-v^2-v-2)*X
        -(u^2+u*v+v^2+u+v+1)).
```

The visible rational `2`-torsion together with the distinguished order-`8`
point gives torsion containing `[2,2,8]`.

Main files: `notes/m3222_plus3.md`, `notes/m3222_halving_boundary.md`,
`code/m3222_plus3_search.m`.

### `A(2,2,4,4)` squareclass threefold: `[2,2,4,4]`

Start from the square branch-point full-2 model

```text
C: y^2 = x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2).
```

The additional halving conditions for a second independent order-`4` class can
be written as

```text
y1^2 = (a^2 - c^2)*(a^2 - d^2),
y2^2 = (b^2 - c^2)*(b^2 - d^2),
w^2  = (a^2 - c^2)*(b^2 - c^2).
```

This algebraic threefold is the `A(2,2,4,4)` model used in the local
obstruction and boundary searches.

Main files: `notes/a2244_local_obstructions.md`,
`paper/main.tex`, `code/m2248_sieve.m`.

### `M(2,4,4)` fiber product: `[2,4,4]`

The model is

```text
C: y^2 = x*(x+u^2)*(x+v^2)*(x^2 + (t^2+2*s)*x + s^2).
```

The parameters come from the elliptic fiber product

```text
E2: Y^2 = X*(X+t^2)*(X+4*s+t^2),
E3: Y^2 = X*(X^2 - 4*(t^2+2*s)*X + 16*s^2).
```

For `R in E3` and `P=(x1,y1) in E2`, the note uses the explicit dual
2-isogeny formula `phi(R)` and sets

```text
Q = P + phi(R) = (x2,y2),
u = y1/(2*x1),
v = y2/(2*x2).
```

This gives curves with generic torsion `[2,4,4]`.

Main files: `notes/m244_to_248_route.md`, `code/m244_isogeny_probe.m`,
`code/m244_to_248_sample_search.m`.

### `M_1(8,4)` tangent cover: `[4,8]`

Start from

```text
C: y^2 = x*A(x)*B(x),
```

where

```text
A = n^4*x^2
  + (m^3*n + 4*m^2*t + m*n^3 - 8*m*n*t + 4*n^2*t)*x
  + m^4

B = (m*n + 2*n^2 + 4*t)*x^2
  + (m^2 + 4*m*n + n^2 + 8*t)*x
  + (2*m^2 + m*n + 4*t).
```

The class `[x,0]` is halved by solving

```text
h(x) - x*(M*x+N)^2 = c4*(x^2 + U*x + V)^2,
```

where `h=A*B` and `c4=lc(h)`.  On the chart `n=1`, the first square condition
is parametrized by

```text
t = (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1)).
```

The remaining tangent equations define the `[4,8]` cover.  The notes record
many exact `[4,8]` examples and simple Frobenius certificates.

Main files: `notes/m18_m14_halving.md`, `code/m18_m14_tangent_search.m`,
`code/m18_m14_torsion_structure_search.m`.

### One-split inside `[4,8]`: `[2,4,8]`

Inside the `[4,8]` tangent-cover family above, require exactly one of the
quadratics `A` or `B` to split over `Q`.  This adds one independent rational
`2`-torsion class and gives exact `[2,4,8]` in the verified examples.

This is recorded as a positive equation-defined subfamily, but the notes do
not give a compact rational parametrization.  The exact height-20 run found
eight verified hits, all with torsion `[2,4,8]`, including a geometrically
simple example.

Main files: `notes/m18_m14_halving.md`,
`code/m18_m14_one_split_search.m`,
`code/m18_m14_one_split_verify_248.m`.

### Quintic contact: `[5]`

The elementary contact family is

```text
h(x) = 1 + a*x + b*x^2,
f(x) = h(x)^2 - (1+a+b)^2*x^5.
```

On `C: y^2=f(x)`, the point `P=(0,1)` satisfies

```text
div(y-h(x)) = 5*P - 5*infinity,
```

so `[P-infinity]` has order dividing, and generically exact, `5`.

Main files: `notes/m10_quintic_contact5.md`,
`code/m10_quintic_contact5_search.m`.

### Contact-6: `[6]`

Use

```text
h6 = 1 + a*x + b*x^2 + x^3,
f  = h6^2 - (x - 1)^6.
```

Then `h6-y` has divisor

```text
6*P - 6*infinity,    P = (1, a+b+2),
```

on the smooth nonboundary open chart.  This gives a marked class of generically
exact order `6`.

Main files: `notes/contact6_m36.md`, `code/contact6_m36_search.m`.

### Contact-7: `[7]`

Take

```text
h = 1 - (7/2)*x + a*x^2 + b*x^3,
f = (h^2 + (x - 1)^7)/x^2.
```

Then the function `x*y-h(x)` has divisor `7*P-7*infinity`, where
`P=(1,h(1))`, on the smooth nonboundary open chart.

Main files: `notes/contact7_family.md`, `code/contact7_family_search.m`.

### Contact-9: `[9]`

Take

```text
h = 1 - (9/2)*x + (63/8)*x^2 - (105/16)*x^3 + a*x^4,
f = (h^2 + (x - 1)^9)/x^4.
```

Then `x^2*y-h(x)` gives a rational divisor class of generically exact order
`9`.

Main files: `notes/contact9_family.md`, `code/contact9_family_search.m`.

### Flynn and Daowsud-Schmidt order-11 families: `[11]`

Two published one-parameter genus-2 families give rational `11`-torsion.
Flynn's family is

```text
C_t: y^2 = x^6 + 2*x^5 + (2*t+3)*x^4 + 2*x^3
             + (t^2+1)*x^2 + 2*t*(1-t)*x + t^2.
```

Daowsud-Schmidt's continued-fraction family is

```text
C_u: y^2 = x^6 - 4*x^5 + 8*(1+u)*x^4 - (10+32*u)*x^3
             + 8*(1+6*u+2*u^2)*x^2
             - 4*(1+6*u+16*u^2)*x + 64*u^2 + 1.
```

Local checks after the internet sweep:

```text
Flynn t=1: TorsionSubgroup(J)(Q) = [11]
Daowsud-Schmidt u=1: TorsionSubgroup(J)(Q) = [11]
```

Both have an irreducible Frobenius polynomial at `p=3`, and Sage/Lombardo
checks give geometric endomorphism ring `ZZ` for several small parameters.

Main files: `notes/internet_infinite_families_sweep.md`,
`code/order11_family_check.m`, `code/order11_family_geom_check.sage`.

### Rational-branch subfamilies of order-11 families: `[22]`

For an even sextic with two rational points at infinity, if
`D_inf = infinity_+ - infinity_-` has order `11` and `W=(r,0)` is a rational
finite Weierstrass point, then

```text
2*(W - infinity_+) = infinity_- - infinity_+ = -D_inf,
```

so `W - infinity_+` has order `22`.

In Flynn's order-11 family, the rational-root condition has discriminant
`16*r^5`.  Setting `r=s^2` gives

```text
t_eps(s) = (-s^2*(s^2+1)*(s^4-s^2+1) + 2*eps*s^5)/(s^2-1)^2,
eps = +/-1.
```

In the Daowsud-Schmidt order-11 family, the rational-root condition has
discriminant `256*(r-1)^5`.  Setting `r=1+s^2` gives

```text
u_eps(s) = (-s^2*(s^2+1)*(s^4-s^2+1) + 2*eps*s^5)/(4*(s^2-1)^2),
eps = +/-1.
```

At `s=2`, all four branches have exact torsion `[22]` and Sage/Lombardo
geometric endomorphism ring `ZZ`.

Main files: `notes/order22_from_order11.md`,
`code/order22_from_order11_check.m`,
`code/order22_from_order11_geom_check.sage`.

### Kuru-Sadek genus-2 quadratic-order family: `[23]`

Kuru-Sadek's quadratic-order construction gives, after the rationalizing
specialization in their corollary, a one-parameter genus-2 family with a
rational point of order `23` on the Jacobian.  For `t != 0, +/-1`, set

```text
beta  = (t^2 + 1)^2/(4*t^2)
sbeta = (t^2 + 1)/(2*t)
s      = (t^2 - 1)/(2*t)
alpha = beta - s^5/(beta*sbeta)
lambda = (alpha - 1)^4/((alpha - beta)^2*alpha)
```

and

```text
A(x) = (x^3*(x-alpha)^2
        - (x-1)*((x-1)^4 - lambda*(x-beta)^2*x))
       /(2*(x-alpha)*(x-beta)),
C_t: y^2 = A(x)^2 - lambda*x^4*(x-1).
```

The TeX source line for `alpha` appears to omit the exponent on `beta`; the
formula above is the interpretation that reproduces the published genus-2
example at `t=2`.

At `t=2`, the primitive integral model is

```text
y^2 = -299054816676000*x^5
      + 937313042871529*x^4
      - 1165161421194050*x^3
      + 677279473485625*x^2
      - 132825168000000*x
      + 8294400000000.
```

Magma gives `TorsionSubgroup(J)(Q) = [23]`, and an irreducible Frobenius
polynomial at `p=13`.  Sage/Lombardo gives geometric endomorphism ring `ZZ`.

Main files: `notes/internet_infinite_families_sweep.md`,
`code/order23_kuru_sadek_check.m`,
`code/order23_kuru_sadek_geom_check.sage`.

### `M(12)`: `[12]`

The base `M(12)` model is

```text
y^2 + (x-r)*(T+1)*y = a*x^2*T*(T+1),
T = a*x^2 - x + r,
P = (0,0).
```

The marked point gives a Jacobian class of order `12`.  Completing the square
gives

```text
Y^2 = ((x-r)*(T+1))^2 + 4*a*x^2*T*(T+1).
```

Main files: `notes/m12_simple_route.md`, `code/m12_simple_search.m`,
`code/m212_construct_3torsion.m`.

### One-parameter `Z/12 x Z/2`: `[2,12]`

Impose a split `T+1` and an extra rational root of the residual quartic.  The
simple component recorded in the notes is

```text
a = (1-r)/4,
u = -2*r/(r-2).
```

Then

```text
T+1 = (x-2)*(x + (2*r+2)/(r-1)),
```

and the completed-square polynomial factors as

```text
W = (x-2)*(x + 2*r/(r-2))*(x + (2*r+2)/(r-1))
    *(x^3 + 4/(r-1)*x^2 + (-r^2-5*r-2)/(r-1)*x
       + (2*r^2+2*r)/(r-1)).
```

This gives torsion containing `Z/12Z x Z/2Z`, and the notes include simple
examples.

Main files: `notes/m12_simple_route.md`, `code/m12_z12x2_family.m`.

### Contact-7 rational-root subfamily: `[14]`

Start from the contact-7 family and force a rational Weierstrass point by

```text
r = 1 - s^2,
h(r) = eps*s^7,    eps = +/-1.
```

For `r != 0`, this determines

```text
a = (eps*s^7 - 1 + (7/2)*r - b*r^3)/r^2.
```

The marked `7`-torsion class and the rational `2`-torsion class give torsion
containing a point of order `14`.

Main files: `notes/contact7_family.md`, `code/contact7_root_even_search.m`.

### Contact-9 rational-root subfamily: `[18]`

Start from the contact-9 family and force a rational Weierstrass point by

```text
r = 1 - s^2,
h(r) = eps*s^9,    eps = +/-1.
```

This determines

```text
a = (eps*s^9 - (1 - (9/2)*r + (63/8)*r^2 - (105/16)*r^3))/r^4.
```

The marked `9`-torsion class and the rational `2`-torsion class give torsion
containing a point of order `18`.

Main files: `notes/contact9_family.md`, `code/contact9_family_search.m`.

### Contact-5 plus 4-torsion: `[20]`

Inside the quintic-contact family, the halving component is

```text
b = (a^2 - 1)/2.
```

Writing `t=a`,

```text
h = 1 + t*x + ((t^2 - 1)/2)*x^2,
f = h^2 - ((t + 1)^4/4)*x^5.
```

For every nonsingular specialization,

```text
H = [x^2 + 2*x/(t+1), (t+2)*x + 1]
```

satisfies `2H=[x-1,0]`, so the contact `5`-torsion and this order-`4` class
combine to give a point of order `20`.

Main files: `notes/contact5_order40_family.md`,
`code/contact5_order40_family_search.m`.

### Extra-2 contact-5/order-20 loci: `[2,20]`

The independent extra `2`-torsion condition is that the residual quartic
`f/(x-1)` be reducible.  With

```text
u = t + 1,
y = u*x,
```

the residual quartic is equivalent to

```text
u*(y^4 + 4*y^3 + 8*y^2 + 8*y + 4) - 4*y*(y+1)^2.
```

The `1+3` locus is parametrized by

```text
t = -(z^4 + 4*z + 4)/(z^4 + 4*z^3 + 8*z^2 + 8*z + 4).
```

The `2+2` locus is parametrized by

```text
t = -(r^6 - 2*r^5 + 2*r^4 - 4*r^3 + 4*r^2 - 8*r + 8)
     /((r^2 - 2)^2*(r^2 - 2*r + 2)).
```

On smooth nonboundary specializations, these give torsion containing
`Z/2Z x Z/20Z`.  The exact searches found many `[2,20]` curves and a special
simple `[2,2,20]` point, but `[2,2,20]` itself is not currently recorded as an
infinite family.

Main files: `notes/contact5_order40_family.md`,
`code/contact5_order20_extra2_search.m`,
`code/contact5_extra2_param_large_search.m`.

### Simultaneous contact-5/contact-6: `[30]`

Seek

```text
f = h6^2 - (x-1)^6 = h5^2 - K*x^5,
```

where

```text
h6 = x^3 + A*x^2 + B*x + C,
h5 = e*x^2 + d*x + c.
```

The coefficient equations reduce to a genus-zero branch.  The implemented
parametrization uses

```text
t = (5*R^2 - 20*R + 19)/(R^2 - 5),
Y = -2*(5*R^2 - 22*R + 25)/(R^2 - 5),
u = t^3,
s = t^5 + t^4 + (5/2)*t^3 + (1/2)*t
    +/- t*(t-1/2)*(t+1)*Y,
```

then recovers `q,A,B,C,e,d,c` as in
`notes/contact5_contact6_order30_family.md`.  The curve is

```text
C_R: y^2 = (x^3 + A*x^2 + B*x + C)^2 - (x-1)^6.
```

For tested smooth nonboundary fibers,

```text
Order(D5) = 5,
Order(D6) = 6,
Order(D5+D6) = 30,
TorsionSubgroup(J)(Q) = [30].
```

The family has simple specializations, so it is not contained in the
decomposable locus; the generic member is geometrically simple outside a
proper exceptional set.

Main files: `notes/contact5_contact6_order30_family.md`,
`code/contact5_contact6_order30_family.m`.

### Elkies Clebsch-Klein family: `[2,2,2,10]`

Elkies' full-level-2 atypical `5`-torsion family is

```text
C: y^2 = x * prod_i (x - r_i^2),
```

where

```text
sum_i r_i = 0,
sum_i r_i^3 = 0.
```

For primitive integral representatives with nonzero, pairwise distinct
squares, the curve has full rational `2`-torsion and a rational `5`-torsion
point.  The generic exact rational torsion in the source family is therefore
`[2,2,2,10]`.

Main files: `notes/elkies22210_richelot.md`,
`code/elkies22210_richelot_sweep.m`.

## Positive covers/examples not counted as established infinite families

These are useful, but I would not list them in the main "known infinite
families" table without further parametrization or a rational-point argument.

- `[28]`: the contact-7 rational-root first-halving surface has explicit
  equations and recovers the known simple `[28]` example
  `(a,b)=(11/2,-7/2)`.  See `notes/contact7_family.md` and
  `code/contact7_halving_surface_search.m`.  The notes do not record a
  parametrized infinite `[28]` family.

- `[3,12]`: extra independent `3`-torsion on `M(2,12)` has split/nonsimple
  examples, including `z=-5/3, r=-3/5`, but no simple infinite family is
  recorded.  See `notes/m212_three_torsion.md`.

- `[2,2,20]`: there is a geometrically simple exact example at
  `z=-1/7` in the contact-5 extra-2 locus, but the double-linear analysis
  found only this nondegenerate point in the searched range.  It is not
  currently an infinite family.  See `notes/contact5_order40_family.md` and
  `notes/how_we_found_2220_examples.md`.

- `[40]`: the contact-5/order-20 family has the isolated order-40
  specialization `t=-1/3`.  The order-40 cover appears high genus and no
  infinite order-40 family is recorded.  See `notes/contact5_order40_family.md`.

- `[6,6]`, `[8,8]`, `[4,16]`, `[35]`, `[60]`, `[72]`, `[80]`, `[120]`: the
  repository has searches, local obstructions, or algebraic cover equations
  for several of these targets, but no confirmed infinite family is currently
  recorded.
