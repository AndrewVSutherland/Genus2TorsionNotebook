# Contact-6 `[6,12]`: exact geometry of the weighted endpoint `R3` cover

Date: 2026-07-10.

## Conclusion

After saturation by the square-quartic open condition `s4 != 0`, the exact
endpoint `R3` halving curve has exactly two prime one-dimensional components
over `Q`.

1. A rational component projects two-to-one to an irreducible degree-`8`
   genus-zero plane curve `P8(e,mu)=0`.  Its full normalization is the conic

   ```text
   6*t^2 + y^2 = 100,
   ```

   so this component has an explicit rational parametrization and genuine
   rational open `R3` halves.

2. The component carrying the four weighted `E9` core residue disks is
   birational to an irreducible degree-`16` plane curve `P16(e,mu)=0` of
   normalized genus one.  Its normalization has the minimized binary-quartic
   model

   ```text
   Y^2 = 2*Z^4 - 6*X^4.
   ```

   This curve has no `Q_2` point and no `Q_3` point, although it has the
   previously observed `Q_5` points.  Therefore none of the smooth `E9`
   core-plus-`R3` branches can globalize to a rational `[6,12]` curve.

The rational `P8` component is real geometry, not discarded boundary.  At an
explicit point its source and dual torsion are respectively `[2,2,6]` and
`[2,12]`, and the `R3` class is exactly divisible by two.  It misses `[6,12]`
because both sides have only one rational `3`-direction.  Thus `P8` remains a
possible family to intersect with other contact-core strata, but it does not
rescue the weighted `E9` lane.


## Exact saturated decomposition

Use the notation of `contact6_m612_weighted_E9_R3.m`:

```text
S = R3*(mu*x+nu)^2 - D*A1*A2,
H1 = 8*s4^2*s1 - s3*(4*s4*s2-s3^2),
H0 = 64*s4^3*s0 - (4*s4*s2-s3^2)^2.
```

The raw equations have gcd one.  Exact saturation

```text
I_sat = <H1,H0> : <s4>^infinity
```

takes about 24 seconds and gives a one-dimensional ideal with a 32-element
Groebner basis.  A primary decomposition, externally capped at 300 seconds,
completed in 191 seconds and returned exactly two components.  Both are prime
and one-dimensional.

Eliminating `nu` is much smaller.  Up to a nonzero rational constant,

```text
Res_nu(H1,H0)
 = e^4*(e+2/3)^2*(e+3/5)^2*s4(e,mu)^8
   *P8(e,mu)^2*P16(e,mu).
```

The factors in `e` and `s4` are exceptional projection or chart-boundary
contributions, not components of `I_sat`.  The two remaining factors are
irreducible over `Q`:

```text
P8:  total degree 8,  degree_e 8,  degree_mu 4,  10 terms;
P16: total degree 16, degree_e 16, degree_mu 8, 35 terms.
```

A convenient integral multiple of `P8` is

```text
1250*e^8 + 3000*e^7 + 2700*e^6 + 1080*e^5 + 162*e^4
- 75*e^4*mu^2 - 65*e^3*mu^2 + 3*e^2*mu^2
+ 9*e*mu^2 + 3*mu^4.
```

Over the function field `Q(P8)`, the gcd of `H1,H0` as polynomials in `nu`
has degree two.  Over `Q(P16)`, it is linear.  Thus the first space component
is a quadratic cover of its plane projection, while the second is birational
to `P16`.  The projective embedding degrees of the original space components
were not separately computed; the certified projection degrees and normalized
genera are the useful invariants here.


## The rational `P8` component

The normalization of the `(e,mu)` projection is rational.  With parameter
`t`, one model is

```text
e(t) = -(25/3)*t^2/(t^4-25*t^2+1250/3),

mu(t) = (-5*t^7+(1750/9)*t^5-(6250/3)*t^3)/D(t),

D(t) = t^8-50*t^6+(4375/3)*t^4-(62500/3)*t^2+1562500/9.
```

The monic quadratic recovery polynomial for `nu` has linear coefficient

```text
B(t) = (-(500/3)*t^5+(175000/27)*t^3-(625000/9)*t)/D(t).
```

Its discriminant is

```text
-6*(t^2-50/3)*h(t)^2,

h(t) = (t^7-(500/9)*t^5+(28750/27)*t^3-(62500/9)*t)/D(t).
```

Consequently the full space component normalizes to

```text
y^2 = 100-6*t^2,
nu = (-B(t)+h(t)*y)/2.
```

Parametrize this conic through `(t,y)=(4,-2)` by

```text
t(u) = 4*(u^2+u-6)/(u^2+6),
y(u) = 2*(u^2-24*u-6)/(u^2+6).
```

Substitution gives exact identities `H1=H0=0` in `Q(u)` and a nonzero
rational function `s4(u)`.  At `u=12` this gives the compact exact point

```text
(e,mu,nu) = (-200/409, -36320/167281, 38136/167281).
```

Direct substitution verifies `H1=H0=0` and `s4 != 0`.

For the corresponding endpoint source

```text
(a,b) = (1/e,0) = (-409/200,0),
```

exact Jacobian arithmetic gives

```text
source torsion: [2,2,6],
distinguished dual torsion: [2,12],
halves of (R1,R2,R3): [false,false,true].
```

Both curves are squarefree and geometrically simple: at `p=23` their common
irreducible Frobenius polynomial is

```text
529*T^4 + 22*T^2 + 1.
```

Only one invariant factor on either Jacobian is divisible by `3`.  Hence this
is a genuine simple rational `R3`-halving specialization, but not a contact
`[3,3]` core and not a `[6,12]` hit.


## The `P16` component and the global obstruction

The plane curve `P16` is irreducible and has normalized genus one.  Its four
rational plane singularities resolve as follows:

```text
plane point       degrees of places above it
(-3/5:0:1)        [8]
(-6/17:0:1)       [2]
(0:0:1)           [2,2]
(0:1:0)           [8]
```

Thus none of the visible rational singular points hides a rational point on
the normalization.  A plane-point search through bound 500 finds only these
four points.

Let `D2` be the degree-two place over `(-6/17:0:1)`.  Exact Riemann--Roch
dimensions are

```text
deg(D2)=2,  l(D2)=2,  l(2D2)=4,  l(4D2)=8.
```

The returned basis of `L(D2)` is explicitly checked to be
`[nonconstant,1]`.  Choosing `x` from that basis and a complement `y` in
`L(2D2)`, the nine functions

```text
1,x,x^2,x^3,x^4,y,x*y,x^2*y,y^2
```

have a unique relation in `L(4D2)`.  Completing the square and applying
Magma's exact minimization gives

```text
Y^2 = 2*Z^4 - 6*X^4.
```

Its Jacobian is

```text
E: y^2 = x^3 + 3*x,
rank(E)=1,  E(Q)_tors = Z/2Z.
```

The rational quotient by `mu -> -mu` is also genus one.  Using
`z=mu^2`, it has the elliptic model

```text
E_quot: y^2 = x^3 - 12*x,
rank(E_quot)=1,  E_quot(Q)_tors = Z/2Z.
```

This agrees with the expected rational two-isogeny between the two
Jacobians.


## Elementary local proof

On the minimized quartic, scale a local point so that `X,Z` are integral and
`min(v_p(X),v_p(Z))=0`.

At `p=2`, a fourth power is `0` or `1 mod 16`.  According as only `X`, only
`Z`, or both are odd, the right side

```text
2*Z^4 - 6*X^4
```

is `10`, `2`, or `12 mod 16`.  None is a square modulo 16.  Hence there is no
`Q_2` point.

At `p=3`, if `Z` is a unit, reduction gives `Y^2=2`, impossible modulo 3.
If `3|Z`, then `X` is a unit and the right side has valuation exactly one,
again impossible for a square.  Hence there is no `Q_3` point.

At `p=5`, `(X,Z,Y)=(1,0,2) mod 5` is a smooth point and Hensel lifts.  Magma
independently reports

```text
Q_2: false,  Q_3: false,  Q_5: true.
```

This is precisely the phenomenon seen in the weighted calculation: the
endpoint core and `R3` equations are simultaneously smooth over `Q_5`, but
their component is globally impossible.


## Why this rules out the `E9` disks

At `e=0`, the two projection factors specialize as

```text
P8(0,mu)  = constant*mu^4,
P16(0,mu) = constant*mu^6*(mu^2+9).
```

The open weighted `R3` points have `mu != 0` and, modulo `5`, satisfy
`mu^2=-9`, giving `mu=1,4`.  Therefore every one of the previously certified
`E9` core-plus-`R3` residue disks lies on `P16`, not on the rational `P8`
component.  The `Q_2/Q_3` obstruction on `P16` eliminates all of them.


## Reproduction

From `torsion_jac`, run:

```text
magma -b code/contact6_m612_weighted_R3_geometry_probe.m
magma -b code/contact6_m612_weighted_R3_geometry_components.m
magma -b code/contact6_m612_weighted_R3_geometry_small_component.m
magma -b code/contact6_m612_weighted_R3_geometry_p8_param.m
magma -b code/contact6_m612_weighted_R3_geometry_p8_audit.m
magma -b code/contact6_m612_weighted_R3_geometry_p16_places.m
magma -b code/contact6_m612_weighted_R3_geometry_p16_quotient.m
magma -b code/contact6_m612_weighted_R3_geometry_p16_degree2_model.m
```

The exact saturation summary is

```text
magma -b code/contact6_m612_weighted_R3_geometry_saturation.m
```

and the heavier primary decomposition should remain externally bounded:

```text
timeout 300s magma -b mode:=decompose \
  code/contact6_m612_weighted_R3_geometry_saturation.m
```


## Next step

The focused weighted `E9` endpoint lane is closed globally.  The remaining
use of this calculation is the rational `P8` family: intersect its explicit
`(e(u),mu(u),nu(u))` parametrization with other exact contact-core equations.
The first audited point shows that `R3` halving alone naturally gives simple
`[2,12]` duals; an intersection must additionally impose the missing second
rational `3`-direction.
