# Exact transverse tangent at the cyclic `[60]` HLP seed

The calculation in `code/hlp_z60_transverse_tangent.py` is based at the
actual cyclic-`[60]` curve

```text
F = -46250000*x^6 + 1761500625*x^4
    -22332312000*x^2 + 94277468160,
```

not at the older HLP `[6,12]` seed.  Magma extracts a generator `T` and

```text
12*T = (x^2-4608/395, -(164280/79)*x, 2),  order 5,
20*T = (x^2-316/25,  -5476,             2),  order 3,
15*T = (x^2-1506/125,-24642/5,          2),  order 4,
30*T = (x^2-1728/125,0,                 2),  order 2.
```

This extraction is reproduced by `code/hlp_z60_marked_extract.m`.

## Exact marked identities

The order-5 class does **not** have an identity of the restricted form
`H^2-F=k*q^5`.  Because the curve has two non-rational points at infinity,
the relevant principal function is `A+B*y`.  Put

```text
q5 = x^2-4608/395,
B5 = 1-(125/1728)*x^2,
A5 = (48752125/110592)*x^5
     -(1558625/144)*x^3 + 66600*x,
k5 = 1778923230671875/4076863488.
```

Then exactly

```text
A5^2-F*B5^2 = k5*q5^5,
A5+B5*(-164280/79*x) = 0 mod q5.
```

For the order-3 class,

```text
A3 = 2775*x^2-29600,
A3^2-F = 46250000*(x^2-316/25)^3.
```

For the order-4 class and its double,

```text
q2  = x^2-1728/125,
u4  = x^2-1506/125,
ell = 2775*q2,
ell^2-F = 46250000*u4^2*q2,
ell = -24642/5 mod u4.
```

## Tangent ranks

Use seven free sextic coefficients.  After fixing the scale of the order-5
principal function by `B5(0)=1`, its norm identity has eleven auxiliary
variables and eleven coefficient equations.  The order-3 and order-4
identities each have seven auxiliary variables and seven equations.  Exact
rational differentiation at the seed gives

```text
auxiliary ranks       = 11, 7, 7,
combined rank         = 25 in 32 variables,
combined tangent dim  = 7.
```

All three auxiliary blocks are invertible.  Consequently projection of the
marked incidence to the seven sextic coefficients is etale at this point:
**every rational first-order sextic deformation has a unique rational lift
of the marked 5-, 3-, and 4-layers.**

The infinitesimal `PGL(2)` action together with equation/`y` scaling has
rank four.  Thus the marked moduli tangent has dimension `7-4=3`, as it
should for genus 2.

## Leaving the split locus

Magma finds that the geometric reduced automorphism group of this seed has
order two; its only nonidentity element is `x -> -x`.  Therefore only one
Humbert-4 branch passes through the point.  In affine sextic coefficient
space its tangent is the even-sextic subspace plus the infinitesimal
`PGL(2)` orbit.  It has rank six (dimension two after the gauge quotient).

An exact primitive normal, on `(df0,...,df6)`, is

```text
N_H4 = (0,81125,0,904800,0,9916416,0).
```

In particular

```text
delta F = x,       N_H4(delta F)=81125 != 0.
```

This is a genuine rational non-split tangent direction, not a coordinate
change.  Its sparse rational lifts are

```text
order 5:
  dq51 = 2209/172868558400,
  dA50 = 20011/2431344,
  dA52 = -11231495/8402724864,
  dA54 = 354441725/6453292695552,
  dB51 = 43/14004541440;

order 3:
  dq31 = -3/1620896000,
  dA31 = 137/202612,
  dA33 = -75/1620896;

order 4:
  dq21 = -1/145880640,
  du41 = 29/4668180480,
  dL41 = 8125/175056768.
```

All omitted first derivatives are zero.

## A small exact auxiliary slice, and its limitation

There is a compact exact transverse base slice obtained by integrating the
three nonzero order-3 tangent coordinates linearly.  Let

```text
q3(t) = x^2-316/25 -(3/1620896000)*t*x,
A3(t) = 2775*x^2-29600
        +(137/202612)*t*x -(75/1620896)*t*x^3,
F(t)  = A3(t)^2-46250000*q3(t)^3.
```

Then `F(0)=F`, `dF/dt at 0 = x`, and the order-3 identity holds exactly for
every `t`.  Pulling the order-5 and order-4 norm equations back to this
explicit `t`-line gives 18 equations in their 18 auxiliary variables.  Their
Jacobian at `t=0` is invertible, so the full marked pullback is an exact
algebraic curve over `Q`, smooth at the HLP point, and has a unique rational
formal branch over `Q[[t]]`.

This does **not** parametrize the full marked curve by rational functions.
Only the curve and order-3 layer above are explicitly rational in `t`; the
order-5 and order-4 layers form a finite algebraic pullback.  A rational
number `t != 0` yields a cyclic-`60` specialization only when both finite
blocks acquire rational points simultaneously.  The tangent calculation
proves that no geometric/split obstruction prevents this, but it does not
produce a second rational point.

## Reproduction

```text
python3 code/hlp_z60_transverse_tangent.py
magma -b code/hlp_z60_marked_extract.m
```

Both jobs use negligible memory; the Magma extractor is capped at 2 GB.
