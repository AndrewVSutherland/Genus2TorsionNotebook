# Target `[2,2,4,8]`: family-level Richelot sieve, second stage

Date: 2026-07-18

## Outcome

This stage did **not** find a geometrically simple genus-2 Jacobian with exact
torsion `[2,2,4,8]`.  It did complete the intended family-level control and
two bounded reverse-source experiments, and it found two new exact
`[2,4,8]` one-split curves.  Both new curves are unusable for the simple
target: each has a rational Richelot quotient which Magma exhibits explicitly
as a Cartesian product of two elliptic curves.  Hence their Jacobians, their
row-51 base, and every Jacobian in that isogeny class are not geometrically
simple.

The exact bounded results are:

| lane | bases/rows | Richelot edges | one-split | reduction gate | exact torsion | usable simple `[2,4,8]` sources |
|---|---:|---:|---:|---:|---:|---:|
| recorded one-split controls | 4 polynomials, 2 Q-isomorphism classes | 3 forward | 2 forward non-full | 0/1 full | 1 full | 1 |
| primitive square branches | 100 bases | 1,500 | 293 | 77 | 77, all `[2,2,2]` | 0 |
| `tor2244`, rows 1--100 | 100 rows | 1,500 | 600 | 558 | 556 `[2,2,4]`; 2 `[2,4,8]` | 0 |

In addition, the direct full-cover check of all 26,653 `tor2244` rows and
all 799,590 ordered Weierstrass charts found no full-cover witness.

The two `tor2244` hits occur at the same base, row 51, and are distinct over
`Q`.  Sending all of their rational Richelot kernels forward gives, for each
source, one elliptic-product quotient, one one-split Jacobian with exact
torsion `[2,4,8]`, and one full-Weierstrass Jacobian with exact torsion
`[2,2,4,4]`.  Thus the hoped-for full-2 branch again loses the 8-chain.

## 1. Family-level sieve

The main implementation is
[the family sieve](../code/target_2248_family_richelot_sieve.m).  It has two
stages.

1. It constructs fully split curves

   ```text
   y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2)
   ```

   from primitive integer tuples, enumerates all 15 rational Richelot
   splittings, and retains one-split quotients.  A two-good-prime reduction
   gate requires `v_2(gcd(#J(F_p))) >= 6` before exact torsion is computed.
   Exact `[2,4,8]` candidates must then pass the strict Frobenius root-power
   simplicity test.
2. It sends every accepted source through all rational Richelot kernels.
   Full rational `J[2]` is read from the factorization of the completed
   hyperelliptic polynomial.  The target gate is
   `v_2(gcd(#J(F_p))) >= 7`; exact torsion is computed after the gate.  The
   `exact_all_full` switch computes exact torsion even on rejected full-2
   neighbors and was enabled in the control run.

The script normalizes curves, deduplicates over `Q` using genus-2 invariants
and `IsIsomorphic`, records exact torsion, and prints a strict simplicity
certificate consisting of a good prime, irreducible Frobenius polynomial,
and degree 4 for the minimal polynomial of `pi^n` for every `n=2,...,12`.

## 2. Recorded one-split controls

The eight sign/reciprocal hits represented in the earlier height-20 search
collapse first to four displayed polynomials and then to two Q-isomorphism
classes.  Both classes have exact torsion `[2,4,8]`.  Exactly one has the
required strict simplicity certificate, namely

```text
y^2 =
  7061463847622250*x^5
  +104632219276049025*x^4
  +135735215960638800*x^3
  +188573481843278400*x^2
  +51200550567936000*x.
```

At `p=47` its Frobenius polynomial is

```text
T^4 - 4*T^3 + 30*T^2 - 188*T + 2209,
```

and all root-power degrees for `n=2,...,12` are 4.  Its three rational
Richelot neighbors have factor-degree patterns

```text
[1,1,1,1,1,1], [1,1,1,1,2], [1,1,1,1,2].
```

The unique full-Weierstrass neighbor fails the target reduction gate already
at `p=19`, where `#J(F_19)=448` has 2-adic valuation 6.  Because
`exact_all_full=true`, it was nevertheless computed exactly and has torsion
`[2,2,2,4]` of order 32.  No target occurs.  The complete output is in
[the control log](../results/target_2248_family_richelot_sieve_controls.log).

The reproducing command, run from `torsion_jac`, was

```text
magma -b tuple_bound:=4 max_bases:=1 max_sources:=10 \
  include_controls:=true exact_all_full:=true write_log:=true \
  log_file:=results/target_2248_family_richelot_sieve_controls.log \
  code/target_2248_family_richelot_sieve.m
```

## 3. Clean 100-base square-branch experiment

The clean bounded run examined exactly 100 Q-isomorphism-distinct primitive
square-branch bases.  It enumerated exactly 1,500 Richelot edges; 1,493 had
Jacobian codomain, 293 were one-split, and 77 passed the two-prime
`v_2 >= 6` gate.  Exact `TorsionSubgroup` was run on all 77 survivors.  Every
one was `[2,2,2]` of order 8, so there was no `[2,4,8]` source.

There was an off-by-one reporting bug in the earlier file
`target_2248_family_richelot_sieve_pilot50.log`: it reported
`bases_used=51` but only enumerated 750 edges, i.e. exactly 50 bases.  The old
code incremented the next-base counter before applying `max_bases`.  The
limit is now checked before append/increment.  The clean run consequently
reports the mutually consistent values

```text
bases_used = 100, edges = 1500.
```

The clean output is in
[the 100-base log](../results/target_2248_family_richelot_sieve_square_rows100.log).
Its command was

```text
magma -b tuple_bound:=14 max_bases:=100 max_sources:=20 \
  include_controls:=false exact_all_full:=false write_log:=true \
  log_file:=results/target_2248_family_richelot_sieve_square_rows100.log \
  code/target_2248_family_richelot_sieve.m
```

This lane is a poor source generator near the origin: its surviving
one-split quotients generically have only rational `[2,2,2]`, far below the
needed `[2,4,8]`.

## 4. Bounded `tor2244` reverse-source sweep

The companion program
[the `tor2244` source probe](../code/target_2248_tor2244_source_probe.m)
reads the repository bank of full-split `[2,2,4,4]` tuples.  For each base it
enumerates all 15 Richelot neighbors, retains the one-split factor pattern,
applies the same two-good-prime `v_2 >= 6` gate, and runs exact
`TorsionSubgroup` on every survivor.

The completed chunk was precisely rows 1 through 100:

```text
rows                 100
Richelot edges      1500
one-split edges      600
reduction passes     558
exact torsion tests  558
exact [2,2,4]        556
exact [2,4,8]          2
```

There was no candidate deduplication before these exact calls: these are
edge counts, not moduli counts.  The full verbose certificate is
[the rows 1--100 log](../results/target_2248_family_richelot_sieve_tor2244_rows1_100.log).
The command was

```text
magma -b bank_file:=paper/scripts_and_data/tor2244.txt \
  start_row:=1 max_rows:=100 write_log:=true \
  log_file:=results/target_2248_family_richelot_sieve_tor2244_rows1_100.log \
  verbose_candidates:=false code/target_2248_tor2244_source_probe.m
```

An extended rows 101--500 run is separate work in progress; it is not part of
the completed counts above.

## 5. Direct full-cover check in every Weierstrass chart

The Richelot search is only one route to the target.  As an independent
check, [the direct full-cover driver](../code/target_2248_tor2244_all_charts.m)
applies the existing exact equations for

```text
M(2,2,4,8) -> A(2,2,4,4)
```

to every row of `tor2244.txt`.  A raw tuple is not necessarily in the correct
normalization: one must choose an ordered pair among all six Weierstrass
points and move it to `0,infinity`.  The driver therefore tests all 30 such
ordered choices, retains every resulting square-normalized four-tuple, and
then tests all 24 finite-root permutations and all cover sign choices.

The exhaustive bank result is

```text
bank rows                         26,653
ordered (zero,infinity) charts   799,590
square-normalized charts         160,018
full-cover source rows                 0
full-cover charts                      0
full-cover witnesses                   0
```

The candidate file is empty.  This is stronger than the standard-chart
check, which also found zero witnesses but left a genuine coordinate gap.
The complete-chart program closes that gap for this finite bank.  It is not
an exhaustion of all rational points of the threefold.

The known Howe--Poonen--Leprevost split `[2,2,4,8]` examples are positive
controls for exactly this normalization.  The existing
`m2248_hpl_normalization_check.m` finds, for each example, six square
normalizations, four full-cover charts, and 32 full witnesses.  Thus the
zero bank result is not caused by omitting nonstandard Weierstrass charts or
by a vacuous full-cover test.

Artifacts:

* [complete-chart driver](../code/target_2248_tor2244_all_charts.m)
* [complete-chart log](../results/target_2248_tor2244_all_charts.log)
* [empty candidate file](../results/target_2248_tor2244_all_charts_candidates.txt)
* [HPL positive-control summary](../data/m2248_hpl_normalization_summary.txt)

## 6. The two row-51 `[2,4,8]` sources

Both hits come from row 51,

```text
[a,b,c,d] = [55216,56550,62234,64090],
```

whose normalized base is

```text
y^2 = x^5 + 14227308012*x^4
    + 75511013334935609136*x^3
    + 177186449725214688497936641600*x^2
    + 155106885219680453949408579750144000000*x.
```

The base has exact torsion `[2,2,4,4]`.  Edges 13 and 14 have exact torsion
`[2,4,8]` and factor-degree pattern `[1,1,1,1,2]`.  They are not isomorphic
over `Q`.  Their reverse kernels are identified explicitly as

```text
edge 13: L1*L2 | L3*L4 | L5*infinity
edge 14: L1*L3 | L2*L4 | L5*infinity.
```

The strict root-power search returns no certificate for the monic base or
either source.  That failure alone would not prove nonsimplicity.  Here there
is a conclusive independent certificate: on each source, Magma returns one
rational Richelot codomain of type `SetCart` and prints it as

```text
Cartesian Product<Elliptic Curve, Elliptic Curve>.
```

Therefore each source Jacobian admits a rational `(2,2)`-isogeny to an
elliptic product.  It is not geometrically simple; by isogeny invariance, the
row-51 base and all of these neighbors are also nonsimple.

An independent cross-check in
[the row-51 audit note](target_2248_source51_audit.md) reaches the same
conclusion from explicit degree-2 elliptic subcovers of the base and also
records compact factored equations for both large source sextics.

The focused verifier
[the row-51 forward script](../code/target_2248_source51_forward.m)
reconstructs the two curves rather than trusting copied large coefficients,
performs Q-isomorphism deduplication, repeats exact torsion and strict
root-power tests, identifies the reverse kernels, prints the two elliptic
products, and identifies every nondegenerate forward kernel.  It runs exact
torsion on **every** rational Jacobian neighbor, a stronger control than
`exact_all_full`.

For each of the two sources the forward result is:

| rational kernel | codomain | full rational `J[2]` | `v_2 >= 7` gate | exact torsion |
|---|---|---:|---:|---|
| edge 1 | elliptic product (`SetCart`) | -- | -- | -- |
| edge 2, `Q \| L1*L3 \| L2*L4` | genus-2 Jacobian | no | pass | `[2,4,8]` |
| edge 3, `Q \| L1*L4 \| L2*L3` | genus-2 Jacobian | yes | pass | `[2,2,4,4]` |

Thus across both sources there are six rational Richelot edges, four
Jacobian codomains, two full-Weierstrass codomains, four reduction-gate
passes, four exact computations, and zero targets.  Full equations, quadratic
factors, determinants, formula scales, reduction orders, exact torsion, and
both elliptic-product equations are in
[the focused row-51 log](../results/target_2248_source51_forward.log).

The reproducing command was

```text
magma -b write_log:=true \
  log_file:=results/target_2248_source51_forward.log \
  code/target_2248_source51_forward.m
```

## Interpretation and next experiment

The family sieve is behaving coherently but exposes a sharp obstruction:

* the certified-simple `[2,4,8]` control reaches full rational `J[2]` only
  after its 8-chain falls to a 4-chain (`[2,2,2,4]`);
* the generic small square-branch reverse lane falls much farther, to
  `[2,2,2]`;
* the `tor2244` lane can preserve/gain enough 2-power torsion to make exact
  `[2,4,8]`, but the first two hits lie in a decomposable isogeny class, and
  their full-2 neighbors are only `[2,2,4,4]`.

The best continuation is therefore the already-separated extension of the
`tor2244` sweep, but with a mandatory simplicity audit on every exact
`[2,4,8]` hit (preferably first on the monic base, since simplicity transfers
through the Richelot isogeny), followed immediately by the focused forward
test used here.  The exact next chunk command is

```text
magma -b bank_file:=paper/scripts_and_data/tor2244.txt \
  start_row:=101 max_rows:=400 write_log:=true \
  log_file:=results/target_2248_family_richelot_sieve_tor2244_rows101_500.log \
  verbose_candidates:=false code/target_2248_tor2244_source_probe.m
```

Any new hit should be rejected at once if a `SetCart` quotient appears.  A
strict-simple hit should instead have all rational kernels pushed forward,
with exact torsion forced on every full-Weierstrass neighbor even when the
2-adic reduction gate rejects it.
