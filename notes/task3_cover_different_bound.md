# Different bound for the corrected degree-16 `[64]` cover

The driver is `code/task3_cover_different_bound.sage`; its compact output is
in `data/task3_cover_different_bound.txt`.  This pass does not rerun either the
180-second maximal-order/genus call or the 90-second global discriminant call.

## Boundary divisor

Write the genus-2 equation as `y^2=x q(x)`, with coefficients `f1,...,f5`
from the normalized `t`-family.  The finite bad-fiber support factors as

```text
A = num(f1),                         degree 8, irreducible and squarefree,
B = num(disc(q)),                    degree 104, irreducible and squarefree,
C = num(f5)
  = const*t^3*(t-1)*(t+1)*(t^2+1)*(t^2+t+1),
Phi = t^4+t^3+t^2+t+1.
```

The factors are pairwise coprime.  Including infinity gives 124 geometric
boundary points, but most of them do not ramify in this particular degree-16
field.

## The two large factors are unramified

The monic degree-16 polynomial `h(t,M)` has good integral coefficients at
`A` and `B`.  Two clean residue certificates are enough to prove separability
over the corresponding characteristic-zero residue fields:

* `A mod 29` has the simple root `t=2`; `h(2,M)` is a product of four distinct
  quartics.
* `B mod 41` has the simple root `t=19`; `h(19,M)` is a product of distinct
  factors of degrees `2,2,6,6`.

Thus all `8+104=112` geometric points over `A B=0` are unramified.  This is the
main improvement over treating every zero of the genus-2 discriminant as a
branch point of the halving field.

## Exploratory local Puiseux information

Exact Singular Puiseux computations over `QQ` give eight branches of index 2
at infinity and at `t=0`; the exact `t=1` computation represents the same 16
Puiseux conjugates with denominator 2.  Tame computations give the same local
type at a root of `Phi` modulo 41 and on the remaining components of `C`
modulo 61:

```text
Phi:             p=41, t=10,
t=-1:            p=61, t=60,
t^2+1:           p=61, t=11,
t^2+t+1:         p=61, t=13.
```

In every case all reported Puiseux denominators are 2.  These exploratory
computations predict
different contribution 8 at each of the seven `C` points, each of the four
`Phi` points, and infinity, hence

```text
deg Different = 12*8 = 96,
2g-2 = -32+96,
g = 33.
```

The compact driver does not reproduce these Puiseux trees.  Moreover, the
number-field Puiseux routine over `QQ(zeta_5)` hit a Singular nested-field bug.
Consequently neither `g=33` nor the intermediate `g<=57` obtained by assigning
contribution 8 at `Phi` and infinity is promoted to an unconditional
characteristic-zero theorem here.  A direct equisingularity proof for the
tame certificates, or local integral bases at `C` and `Phi`, would close that
gap.

## Rigorous bound

The compact certificate proves that all 112 points over `A B=0` are
unramified.  The only remaining possible branch points are the seven distinct
roots of `C`, the four roots of `Phi`, and infinity: 12 geometric points.  In
characteristic zero, a degree-16 cover has local different contribution
`sum(e_P-1) <= 15` at each point.  Therefore

```text
deg Different <= 12*15 = 180.
```

Riemann-Hurwitz gives

```text
g <= (180-30)/2 = 75.
```

This replaces the earlier Newton-polygon bound `g <= 1785` with the fully
reproduced unconditional bound `g <= 75`.  The values `g<=57` and `g=33` are
recorded only as sharply supported targets for the next local certification.

## Reproduction

```text
timeout 30s sage code/task3_cover_different_bound.sage \
  > data/task3_cover_different_bound.txt
```

The exploratory global generic discriminant and modular discriminant calls
were each killed at their 120-second caps and are not part of reproduction.
