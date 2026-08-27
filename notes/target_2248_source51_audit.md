# Independent audit of the row-51 `[2,4,8]` Richelot sources

Date: 2026-07-18

## Conclusion

The two hits at row 51, Richelot edges 13 and 14, are genuine and distinct
curves with exact rational torsion `[2,4,8]`.  They are not `Q`-isomorphic
and are not quadratic twists.  They are, however, **geometrically
nonsimple**: both are `(2,2)`-isogenous to the same split row-51 base.

Enumerating all rational Richelot kernels on both sources gives six edges:
two product codomains and four genus-2 Jacobian codomains.  Exactly two of
the latter have full rational Weierstrass splitting.  Exact Magma torsion
computations give `[2,2,4,4]` on both full-splitting neighbors, so there is
no `[2,2,4,8]` target in this component.

The independent verifier and complete transcript are:

- `code/target_2248_source51_audit.m`
- `results/target_2248_source51_audit.log`

The verifier deliberately calls `TorsionSubgroup` on every Jacobian
codomain.  No `v_2(#J(F_p))` gate controls an exact test.

## Common row-51 base

The row in `paper/scripts_and_data/tor2244.txt` is

```text
(a,b,c,d) = (55216,56550,62234,64090).
```

Thus

```text
C0: y^2 = x(x+3048806656)(x+3197902500)
             (x+3873070756)(x+4107528100)
```

or, expanded,

```text
y^2 = x^5 + 14227308012*x^4
      + 75511013334935609136*x^3
      + 177186449725214688497936641600*x^2
      + 155106885219680453949408579750144000000*x.
```

Magma computes

```text
J0(Q)_tors = [2,2,4,4].
```

The strict Frobenius/root-power scan tested 34 good primes through 199 and
found no witness.  That absence by itself would not prove nonsimplicity.
Here Magma also returns two rational degree-2 elliptic subcovers.  Their
minimal models are

```text
E1: y^2 + x*y = x^3 - 10397089570*x - 408045332529100,
E2: y^2 + x*y = x^3 - 55998808590*x + 2103010558586100.
```

The full degree-2 maps are printed in the log.  Consequently

```text
J0  ~_Q  E1 x E2,
```

so `J0` is geometrically nonsimple.  Since geometric simplicity is
isogeny-invariant, every Jacobian in the Richelot component audited here is
also nonsimple.

## The two exact `[2,4,8]` sources

Writing each sextic in factored form keeps the equations readable and is
exact.  Edge 13 gives

```text
C13: y^2 = A13(x)
             (x+2022167680)(x+3645718180)
             (x+4569338020)(x+6192888520),
```

where

```text
A13(x) = 16633507796891244495522000000*x^2
       + 102441147839805531462183167718600000000*x
       + 156161626691139563033760727015816007500800000000.
```

Edge 14 gives

```text
C14: y^2 = A14(x)
             (x+2174573700)(x+3609306244)
             (x+4605749956)(x+6040482500),
```

where

```text
A14(x) = 12989551065856432869126176768*x^2
       + 82376380487178725718573052003465822208*x
       + 131715816550450032461119529717256680751759360000.
```

For both curves Magma returns

```text
J(C)(Q)_tors = [2,4,8],  order 64.
```

Their absolute genus-2 invariants differ, so they are not even the same
geometric curve; in particular they are neither `Q`-isomorphic nor
quadratic twists.  They are nevertheless `Q`-isogenous: each is linked to
`J0` by a degree-4 Richelot isogeny, so composing through `J0` gives a
degree-16 `Q`-isogeny between them.

## Exhaustive rational forward Richelot sweep

Each source has factorization type `[1,1,1,1,2]`.  Let `Q` denote its
quadratic factor and `L1,...,L4` its four linear factors, ordered as in the
factorizations above.  There are exactly three rational maximal isotropic
kernels.  The outcome is the same for both sources:

| kernel | codomain | exact torsion when a Jacobian | full Weierstrass? |
|---|---|---:|---:|
| `Q | L1*L2 | L3*L4` | product of elliptic curves (`delta=0`) | — | — |
| `Q | L1*L3 | L2*L4` | genus-2 Jacobian | `[2,4,8]` | no |
| `Q | L1*L4 | L2*L3` | genus-2 Jacobian | `[2,2,4,4]` | yes |

The verifier absorbs the sextic leading coefficient into `Q` before using
the determinant/bracket Richelot formula.  This matters: applying the
formula to monic factors without that coefficient can silently produce a
quadratic twist instead of the actual codomain.  Each nondegenerate manual
codomain is checked for `Q`-isomorphism with exactly one Magma built-in
Richelot output.

Across the two sources the four Jacobian outputs collapse to only two
`Q`-isomorphism classes:

1. the two non-full `[2,4,8]` outputs are `Q`-isomorphic to one another;
2. both full `[2,2,4,4]` outputs are `Q`-isomorphic to the original row-51
   base `C0`.

Hence the exact audit counts are

```text
raw [2,4,8] sources                  2
unique source Q-isomorphism classes  2
geometrically simple sources         0
rational forward kernels             6
product codomains                     2
Jacobian codomains                    4
full-Weierstrass Jacobians            2
exact torsion tests                   4
unique forward Q-isomorphism classes  2
[2,2,4,8] target hits                 0
```

There is therefore no target polynomial to report.  More importantly, this
row-51 component should not be used as a geometrically simple `[2,4,8]`
source bank in the continuing search.
