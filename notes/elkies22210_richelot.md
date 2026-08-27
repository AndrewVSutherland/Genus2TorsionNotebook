# Elkies [2,2,2,10] family and Richelot search

This records the test of Elkies' full-level-2 atypical `5`-torsion family as
a possible source of torsion larger than the known `[2,2,20]` examples.

Elkies' Proposition 2.11 gives the split model

```
C: y^2 = x * prod_i (x - r_i^2)
```

where `(r1:...:r5)` lies on the Clebsch-Klein cubic

```
sum_i r_i = 0,     sum_i r_i^3 = 0.
```

For primitive integral representatives with nonzero, pairwise distinct
squares, this gives full rational `2`-torsion and a rational `5`-torsion
point.  The generic exact rational torsion is therefore `[2,2,2,10]`.

The script is

```
code/elkies22210_richelot_sweep.m
```

It enumerates primitive bounded Clebsch-Klein points, computes exact source
torsion, and optionally checks rational Richelot codomains and
`TwoPowerIsogenies` codomains.

## Faster enumeration

The original loop fixed four coordinates and checked the cubic condition.  The
current version fixes `(r1,r2,r3)` and solves for `r4`.  If

```
s = r1+r2+r3,     A = r1^3+r2^3+r3^3,
```

then, with `r5 = -(r1+r2+r3+r4)`, the cubic condition is

```
3*s*r4^2 + 3*s^2*r4 + (s^3-A) = 0.
```

Thus each triple gives at most two integral candidates for `r4`.  This reduces
the source enumeration from an `H^4` loop to an `H^3` loop.

## Checks run

Smoke test:

```
magma -b height:=9 max_sources:=2 max_twopower_sources:=2 \
    code/elkies22210_richelot_sweep.m \
    > data/elkies22210_richelot_h9_m2.txt
```

This rediscovered Elkies' example

```
rs = [1,-8,-7,5,9]
C: y^2 = x(x-1)(x-25)(x-49)(x-64)(x-81)
torsion [2,2,2,10]
```

All 15 Richelot codomains and all 15 `TwoPowerIsogenies` Jacobian codomains
had torsion `[2,10]`.

Height `30`, source-only:

```
magma -b height:=30 max_sources:=100000 run_isogenies:=false \
    code/elkies22210_richelot_sweep.m \
    > data/elkies22210_source_h30.txt
```

This found 11 distinct source square-classes, all with exact torsion
`[2,2,2,10]`.

Height `30`, Richelot and two-power neighborhoods:

```
magma -b height:=30 max_sources:=11 max_twopower_sources:=11 \
    code/elkies22210_richelot_sweep.m \
    > data/elkies22210_richelot_h30_m11.txt
```

All 11 sources had exact torsion `[2,2,2,10]`.  All 165 Richelot codomains and
all 165 `TwoPowerIsogenies` Jacobian codomains had torsion `[2,10]`.

Height `100`, source-only with the quadratic enumerator:

```
magma -b height:=100 max_sources:=100000 run_isogenies:=false \
    code/elkies22210_richelot_sweep.m \
    > data/elkies22210_source_h100.txt
```

Summary:

```
triple_checked 8120601
tuple_checked 514680
clebsch_points 11280
unique_sources 94
best_order 80
best_exponent 10
interesting_count 0
```

All 94 sources again had exact torsion `[2,2,2,10]`.

Height `100`, Richelot-only:

```
magma -b height:=100 max_sources:=100000 max_twopower_sources:=0 \
    code/elkies22210_richelot_sweep.m \
    > data/elkies22210_richelot_h100_no2power.txt
```

This checked all 94 sources and all 1410 rational Richelot codomains.  The
best torsion remained the source torsion `[2,2,2,10]`; every Richelot codomain
listed had torsion `[2,10]`, and `interesting_count` was `0`.

## Exact source-halving audit

The source model has the six finite rational Weierstrass roots

```text
0, r1^2, r2^2, r3^2, r4^2, r5^2.
```

Its `15` nonzero rational `J[2]` classes are therefore the unordered branch
pairs.  For roots `beta_i,beta_j`, the audit constructs the exact Jacobian
class

```text
T_ij = [(x-beta_i)*(x-beta_j), 0]
```

and calls Magma's exact

```text
IsDivisibleBy(T_ij, 2).
```

The implementation is the opt-in `run_source_halving` mode of

```text
code/elkies22210_richelot_sweep.m.
```

It asserts that all `15` constructed classes are distinct, nonzero, and killed
by `2`.  If a half is found, it verifies `2*H=T_ij`, verifies that `H` has
order `4`, and prints the already exact source torsion invariants.  Since the
source also has full rational `J[2]` and a rational `5`-part, any such half
would force a subgroup

```text
[2,2,2,20]
```

of order `160`.

The known Elkies example was audited with

```text
magma -b height:=9 max_sources:=2 run_isogenies:=false \
  run_source_halving:=true code/elkies22210_richelot_sweep.m
```

It found the source tuple

```text
[1,-8,-7,9,5]
```

with exact torsion `[2,2,2,10]`; all `15` source classes were nondivisible.

The complete height-`30` source run was

```text
magma -b height:=30 max_sources:=100000 run_isogenies:=false \
  run_source_halving:=true code/elkies22210_richelot_sweep.m
```

and gave

```text
unique_sources 11
SOURCE_HALVING_SUMMARY sources 11 classes 165 divisible 0 records []
```

The complete height-`100` source run used concise per-class output,

```text
magma -b height:=100 max_sources:=100000 run_isogenies:=false \
  run_source_halving:=true source_halving_verbose:=false \
  code/elkies22210_richelot_sweep.m
```

and gave

```text
triple_checked 8120601
tuple_checked 514680
clebsch_points 11280
unique_sources 94
SOURCE_HALVING_SUMMARY sources 94 classes 1410 divisible 0 records []
```

All `94` sources again had exact torsion `[2,2,2,10]`.  Thus no source in the
height-`100` Clebsch--Klein enumeration has a rationally divisible nonzero
`2`-class, and this source-halving route produces no `[2,2,2,20]` candidate in
the stated box.  This is a bounded exact computation, not a global
nonexistence result on the family.

A compact run record is in

```text
data/elkies22210_source_halving_audit_summary.txt.
```

## Conclusion

This actual Elkies `[2,2,2,10]` family is easy to enumerate and reliably gives
the expected source torsion, but the bounded evidence is negative for finding
larger rational torsion by either accidental specialization or rational
Richelot neighbors.  Through height `100`, no source exceeded order `80`, and
no rational Richelot codomain exceeded order `20`.

The practical conclusion is that this route is less promising than searches
that start from the known geometrically simple `[2,2,20]` example or impose an
independent odd-prime condition on a `[2,20]` family.
