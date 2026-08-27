# Bounded smoke tests for `[2,24]`, `[2,6,6]`, and `[3,12]`

Date: 2026-07-18

This note records fresh bounded reruns of the strongest existing repository
routes for three middle-ranked targets.  None of the reruns found a new
geometrically simple target.  They do, however, cleanly identify the current
wall in each construction and the next calculation that has more value than a
larger blind height box.

## Certificate convention

Several older search scripts call an irreducible quartic local `L`-polynomial
a `SimpleCertificate`.  That is a useful screen, but by itself it is not the
final geometric-simplicity certificate.  For any future target hit, after
computing the exact rational torsion invariants, choose a good prime and check

```text
[ Degree(MinimalPolynomial(pi^n)) : n in [2..12] ]
```

is `[4,4,4,4,4,4,4,4,4,4,4]`, or run an independent geometric endomorphism
test.  No fresh target below reached that final verification stage.

## 1. Target `[2,24]`

### Strongest existing route

The current direct construction is the coprime-composition route on `A(8)`:
the chart supplies a marked order-8 class, one asks for an additional rational
2-class and an independent rational 3-class.  The main files are

```text
code/agent_a2_24_composite8x3.m
code/agent_a2_24_wsplit_3tors.m
notes/agent_a2_24_composite.md
```

The historical dense `W`-split run in the note produced about one million
2-rank-two `A(8)` curves and 1,708 exact `[2,8]` outcomes after the
3-divisibility prefilter, but no rational 3-torsion and no `[2,24]`.

### Fresh smoke test

Command:

```text
magma -b H:=5 NParts:=1 Part:=0 progress:=1000000000 \
  code/agent_a2_24_composite8x3.m
```

Output summary:

```text
SEARCH_DONE tested=54834 smooth=49823 rank2=101 pre3=0
            tors=0 hits24=0 hits224=0
```

Thus the small rerun already reproduces the established separation: the extra
rational 2-class occurs, but none of the 101 rank-two curves survives the ten
good-prime necessary filter for rational 3-torsion.

### Current obstruction

The order-8 and extra-2 conditions are abundant.  The wall is the intersection
with the rational 3-torsion locus.  Existing finite-field work finds no global
local obstruction, so the negative height data should be read as a sparse
codimension-two intersection, not as nonexistence.

### Next high-leverage step

The fastest new experiment is the rational Richelot graph of the new exact
`[2,2,2,12]` record.  Enumerate its 15 rational maximal isotropic kernels and
then breadth-first-search smooth codomains to depth two.  A Richelot isogeny is
prime to 3, so it preserves the rational 3-primary Galois module; geometric
simplicity is also isogeny invariant.  Therefore a neighbor whose 2-primary
torsion is `[2,8]` gives a geometrically simple exact `[2,24]` after only an
exact torsion computation.

If the Richelot graph is cold, eliminate the cubic-contact 3-torsion equations
on the `W`-split `A(8)` sublocus, saturate by all discriminant and split loci,
and compute the irreducible components and their genera.  This is more
informative than extending the existing million-sample search.

## 2. Target `[2,6,6]`

### Strongest existing route

The contact-6 polynomial factors generically as

```text
f = x*((b+3)*x^2 + (a-3)*x + 2)
      *(2*x^2 + (b-3)*x + (a+3)).
```

The `[6,6]` core route adds an independent cubic-contact 3-class.  To obtain
`[2,6,6]`, one additionally splits one quadratic, giving factor type
`[1,1,1,2]`.  The relevant files are

```text
code/contact6_m36_extra_root_curve_search.m
code/contact6_m36_verify_hits.m
notes/contact6_m36.md
```

### Fresh exact hit, but on the nonsimple locus

Command:

```text
magma -b height:=14 prime_bound:=31 max_exact:=100 max_hits:=1 \
  progress_interval:=1000000000 simple_only:=false \
  code/contact6_m36_extra_root_curve_search.m
```

Output:

```text
HIT66 eps 1 r 4/3 a -1/42 b -13/7
factor_degrees [ 1, 1, 1, 2 ]
invs [ 2, 6, 6 ] simple false pcert 0
f = 7112448*x^5 - 36091440*x^4 + 68732496*x^3
    - 58231404*x^2 + 18522000*x

checked 13627  smooth 13535  survivors 11  exact 11  hits 1
```

The exact torsion target is therefore real on this chart, but this known point
is not geometrically simple.  The verifier

```text
magma -b code/contact6_m36_verify_hits.m
```

reconfirmed exact `[2,6,6]` for both known curves.  Every displayed good-prime
`L`-polynomial through `p=67` factors into quadratics; for the curve above,

```text
p=11: L_p = (11*x^2 + 1)^2
p=13: L_p = (13*x^2 - 2*x + 1)^2.
```

This uniform behavior agrees with the previously identified split locus.  The
finite list of factorizations is a strong diagnostic; an explicit elliptic
quotient or extra-involution certificate should be retained as the formal
proof of nonsimplicity.

### Fresh simple-screen run

Command:

```text
magma -b height:=14 prime_bound:=31 max_exact:=20 max_hits:=1 \
  progress_interval:=1000000000 simple_only:=true \
  code/contact6_m36_extra_root_curve_search.m
```

Output summary:

```text
checked 129030  smooth 128035  survivors 92  exact 20  hits 0
```

The `max_exact=20` cap was reached, with no simple-screened exact target.

### Current obstruction

The finite tables printed by both reruns have

```text
p=5: good 12, allowed66 0
p=7: good 36, allowed66 0.
```

Consequently every rational target in this extra-root chart must lie
simultaneously in bad/boundary residue disks at 5 and 7.  The known exact
points land on a split component.

### Next high-leverage step

Work on the cubic-contact `[6,6]` core before specializing blindly:

1. adjoin a square root of the discriminant of exactly one of the two quadratic
   factors;
2. eliminate the cubic-contact variables and saturate by singular, repeated
   root, and bielliptic factors;
3. projectivize the remaining cover and perform simultaneous 5- and 7-adic
   boundary blowups;
4. compute the genera/Jacobian ranks of the surviving components.

The decisive question is whether every live 5/7-boundary component is forced
onto the bielliptic locus.  Answering that is much more valuable than another
height-30 direct run.

## 3. Target `[3,12]`

### Strongest existing route

The marked order-6 contact class is halved on the standard `M(2,12)` chart

```text
m = (1-z^2)/(4*(r+1)),
T = m*x^2 - x + r,
W = (x-r)^2*(T+1)^2 + 4*m*x^2*T*(T+1).
```

The marked class has order 12; the remaining condition is an independent
rational 3-class.  The relevant files are

```text
code/contact6_m36_halveD_m312_search.m
code/contact6_m36_halveD_p5_boundary_analysis.m
notes/contact6_m36.md
```

### Fresh smoke test

Command:

```text
magma -b height:=8 prime_bound:=19 max_exact:=20 max_hits:=5 \
  max_print:=12 progress_interval:=1000000000 \
  code/contact6_m36_halveD_m312_search.m
```

Key output:

```text
prime 5 allowed_312 0 bad 19 rank_counts [ <1, 6> ]
...
checked 7569 smooth 7140 unique_ar 3570
residue_survivors 13 simple_survivors 12
exact_tests 12 hits 0
```

Every one of the twelve simple-screened exact survivors had torsion `[12]`.
The remaining printed survivor was the known point

```text
z=-5/3, r=-3/5, m=-10/9,
```

which has exact `[3,12]` but lies on the nonsimple split locus.

The existing height-40 component analysis is consistent with this smoke test:
45 simple-screened boundary survivors all had exact `[12]`, while the three
exact `[3,12]` points all lay on the split `Rinf+Z0` component.

### Current obstruction

There are no compatible good-open classes modulo 5.  Every rational
`[3,12]` point in this chart must be 5-adically boundary, and all exact points
found so far on that boundary are bielliptic.  Larger affine height bounds do
not address this structural concentration.

### Next high-leverage step

Construct the exact cubic-contact 3-torsion cover over `M(2,12)`, saturate by
the discriminant and the explicit bielliptic/Humbert factor, and decompose the
residual cover.  On each live residual component, carry the existing mod-5
boundary classification one blowup deeper and compute a low-genus model.  If a
genus-1 or genus-2 component remains, determine its Mordell--Weil group and use
a Mordell--Weil sieve or Chabauty.  If the saturated residual cover is empty,
this proves that this particular halved-contact chart produces `[3,12]` only
on the split locus and the project should pivot to a different order-12
family.

## Bottom line

- `[2,24]`: the extra 2-part is easy; rational 3-torsion is the wall.  Test the
  new record's Richelot graph immediately.
- `[2,6,6]`: exact examples are easy, but current ones lie on a split component;
  the target is forced to simultaneous 5/7-boundary disks.
- `[3,12]`: exact examples exist only on the split mod-5 boundary so far;
  decompose the saturated 3-contact cover rather than extending height.
