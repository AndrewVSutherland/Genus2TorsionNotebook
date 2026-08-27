# Contact-7 two-root surface: the visible `[2,14]` base

## Outcome

The contact-7 family contains a rational two-parameter surface on which the
generic rational torsion subgroup is exactly

```text
Z/2 x Z/14.
```

The surface is also generically geometrically simple.  Thus it is a clean
base for imposing one additional `3`-class: the desired `[2,42]` would have
order `84`.

The formulas and all exact checks are implemented in

```text
code/contact7_two_root_surface.m
code/contact7_two_root_boundary.m
```

This note only constructs and certifies the `[2,14]` surface.  It does not
perform the subsequent `+3` search.

## Derivation

Start with

```text
h(x) = 1 - (7/2)x + a*x^2 + b*x^3,
f(x) = (h(x)^2 + (x-1)^7)/x^2.
```

The numerator defining `f` is divisible by `x^2`, and

```text
h(x)^2 - x^2*f(x) = -(x-1)^7.                 (1)
```

To force a rational Weierstrass root, put

```text
r = 1-s^2,       h(r) = eps*s^7,
```

where `eps` is `+1` or `-1`.  Then `f(r)=0`.  Force a second root by

```text
w = 1-t^2,       h(w) = delta*t^7.
```

Writing

```text
Us = eps*s^7 - 1 + (7/2)r,
Ut = delta*t^7 - 1 + (7/2)w,
```

the two linear equations in `a,b` give, away from the raw determinant
boundary,

```text
a = (Us*w^3 - Ut*r^3)/(r^2*w^2*(w-r)),
b = (r^2*Ut - w^2*Us)/(r^2*w^2*(w-r)).        (2)
```

There is a much better cancelled form.  The signs are redundant: replace
`s` by `eps*s` and `t` by `delta*t`, so from now on take both signs to be
`+1` and allow signed parameters.  Define

```text
V(u) = 2u^5 + 4u^4 + 6u^3 + 8u^2 + 10u + 5,
A(u) = V(u)/(2(u+1)^2).
```

Indeed,

```text
u^7 - 1 + (7/2)(1-u^2)
    = (1/2)(u-1)^2*V(u),
(1-u^2)^2 = (u-1)^2*(u+1)^2.
```

Dividing each root equation by `r^2` or `w^2` therefore leaves

```text
a + b*r = A(s),
a + b*w = A(t).
```

Consequently the rational surface has the compact parametrization

```text
b = (A(t)-A(s))/(s^2-t^2),
a = A(s)-b*(1-s^2).                            (3)
```

After cancellation, the common parameter denominator of `a,b` is

```text
(s+1)^2*(t+1)^2*(s+t).
```

In particular, the apparent raw boundaries `s=1`, `t=1`, and `s=t` in
(2) are removable for the rational map (although `s=t` is still a singular
curve fiber because the two roots coincide).

## Generic torsion

Let `C: y^2=f(x)` and let infinity be the unique point at infinity.  From
(1),

```text
(x*y-h(x))*(x*y+h(x)) = (x-1)^7.
```

When the fiber is smooth and `h(1) != 0`, this gives

```text
D7 = [(1,h(1))-infinity],       order 7.
```

The two forced roots give independent classes

```text
Dr = [(r,0)-infinity],
Dw = [(w,0)-infinity],
```

with

```text
ord(Dr)=ord(Dw)=ord(Dr+Dw)=2.
```

Thus the visible subgroup is

```text
(Z/2)^2 x Z/7 = Z/2 x Z/14.
```

The sample below proves that the displayed generic torsion is the *full*
generic rational torsion.  Rational torsion sections specialize injectively
on the smooth characteristic-zero family, while its specialized Jacobian has
exact torsion `[2,14]`.  Hence no additional generic torsion section exists.
In particular, the residual cubic `f/((x-r)(x-w))` is generically
irreducible; otherwise there would be another generic rational `2`-class.

## Exact simple sample

Take

```text
s=2, t=3, eps=delta=1.
```

Then

```text
r=-3, w=-8,
a=691/1440,
b=-5983/1440,
h(1)=-247/40,
```

and

```text
f = x^5 + (21281089/2073600)x^4
        + (17638547/1036800)x^3
        - (11789879/2073600)x^2
        + (3733/160)x - 5609/720.
```

The exact factorization is

```text
f = (x+3)(x+8)
    *(x^3 - (1528511/2073600)x^2
           + (464863/414720)x - 5609/17280),
```

with the cubic irreducible over `Q`.  Scaling `y` by `1440` gives the
integral odd-degree model

```text
Y^2 = 2073600*x^5 + 21281089*x^4 + 35277094*x^3
      -11789879*x^2 + 48379680*x - 16153920.
```

Magma's exact `TorsionSubgroup` computation returns

```text
J(Q)_tors = [2,14],       #J(Q)_tors = 28.
```

At the good prime `p=7`, the Frobenius characteristic polynomial is

```text
Phi_7(X) = X^4 + X^3 - 2X^2 + 7X + 49.
```

It is irreducible and ordinary, has Galois group `D(4)`, and if `pi` is a
root then

```text
[Q(pi^n):Q] = 4 for every 2 <= n <= 12.
```

This is the standard abelian-surface absolute-simplicity certificate.  Thus
the reduction, and hence the Jacobian over `Qbar`, is absolutely simple.
Specialization of geometric endomorphisms then also proves that the generic
Jacobian on the surface is geometrically simple.

For reference, Magma's reduced globally minimal generalized model is

```text
y^2 + (x^2+x)y =
    518400*x^5 + 5320272*x^4 + 8819273*x^3
    -2947470*x^2 + 12094920*x - 4038480.
```

## Exact boundary structure

The reduced parametrization (3) has only the three pole components

```text
s=-1,       t=-1,       s=-t.
```

The first two are the incompatible choices `h(0)=-1` although `h(0)=1`;
the third generically asks for opposite values of `h` at the same root.
The compatible loci `s=1` and `t=1` are removable and can give smooth
curves with a root at `x=0`.

On the parameter chart, put

```text
H(s,t) = s^2*t^2 + 2*s^2*t + s^2 + 2*s*t^2 + 2*s*t
         + (1/2)*s + t^2 + (1/2)*t,
```

and

```text
R(s,t) = s^4*t^2 + 2*s^4*t + s^4
         + 2*s^3*t^3 + 7*s^3*t^2 + 8*s^3*t + 3*s^3
         + 6*s^2*t^3 + 13*s^2*t^2 + 8*s^2*t + s^2
         + 6*s*t^3 + 7*s*t^2 + 2*s*t + 2*t^3 + t^2.
```

The symbolic calculation gives

```text
h(1) = -s^2*t^2*H(s,t)
       / ((s+1)^2*(t+1)^2*(s+t)).
```

The full discriminant is

```text
Disc_x(f) = 8 *
  s^14*t^14*(s-t)^2*H(s,t)^7*R(s,t)^2*R(t,s)^2*K(s,t)
  / ((s+1)^26*(t+1)^26*(s+t)^13),                 (4)
```

where `K` is an irreducible symmetric polynomial of bidegree `(8,8)`.
Its exact coefficients are printed by `contact7_two_root_boundary.m`; they
are left out here only because the expanded polynomial is long.

Thus the singular fibers on the finite parameter chart are exactly supported
on

```text
s=0, t=0, s=t, H=0, R(s,t)=0, R(t,s)=0, K(s,t)=0.
```

Their meanings are transparent:

- `s=0` or `t=0`: a forced root collides with the marked contact point;
- `s=t`: the two forced Weierstrass roots coincide;
- `H=0`: `h(1)=0`, so the contact point degenerates;
- `R(s,t)=0` or `R(t,s)=0`: one forced root becomes multiple;
- `K=0`: the residual cubic develops an internal repeated root.

The boundary script verifies (4) exactly, including the denominator
exponents.  Magma's raw factor list reports scalar `32` because it normalizes
the swapped root factor as `R(t,s)/2`; with `R` normalized above, the scalar
is the displayed `8`.

## Reproduction

```text
magma -b mode:=symbolic code/contact7_two_root_surface.m
magma -b mode:=sample   code/contact7_two_root_surface.m
magma -b code/contact7_two_root_boundary.m
```

The first command verifies the rational formulas, the second verifies the
contact identity, visible divisor orders, full rational torsion, and absolute
simplicity, and the third computes the complete discriminant boundary.
