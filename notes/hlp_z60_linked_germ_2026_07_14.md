# The linked `B5=q0` germ is locally split

At the exact HLP cyclic-`[60]` seed, the quadratic coefficient `B5` in the
order-5 norm identity equals the support `q0` of the rational 2-torsion class
doubled from the marked order-4 class. This note tests whether retaining the
link `B5=q0` can deform away from the Humbert-4 locus.

## Reduced linked incidence

Use the order-4 identity to eliminate the sextic:

```text
f = (q0*L)^2-k4*q0*u^2.
```

After imposing `B5=q0`, the remaining equations are

```text
A5^2-f*q0^2 = k5*q5^5,
H3^2-f       = k3*q3^3.
```

There are 23 normalized auxiliary variables and 18 coefficient equations.
Exact rational differentiation at the HLP seed gives rank 18, so the linked
germ is smooth of dimension 5. Its projection to the seven sextic
coefficients also has tangent rank 5.

## Exact even subgerm and the dimension argument

There is an exact parity-preserving subgerm obtained by taking

```text
q0,u,q5,q3,H3 even; A5 odd; L constant.
```

On this locus the odd coefficient equations vanish identically. The ten
remaining even equations in thirteen variables have Jacobian rank 10, hence
this subgerm is smooth of dimension 3. Every resulting sextic is even, so
every curve on it has the non-hyperelliptic involution `x -> -x` and lies on
Humbert 4.

The two odd infinitesimal `PGL(2)` directions, together with the three sextic
directions from the even subgerm, have rank 5. This equals the full linked
sextic tangent rank. Equivalently, the `PGL(2)`/equation-scaling saturation
of the exact even subgerm has a surjective differential onto the smooth
linked germ. The formal inverse-function theorem therefore shows that this
saturation contains a neighbourhood of the HLP point in the linked germ.

Thus the linked germ is locally contained in Humbert 4 to **all formal
orders**. It cannot leave the split locus at second order. This is a local
statement at the HLP seed; no claim is made about remote components of the
global linked incidence.

As a direct two-jet check, a genuine non-coordinate tangent has primitive
sextic derivative

```text
(11012603904,0,-1737666275,0,68458750,0,0).
```

The forced quadratic lift remains even, and the preserving involution extends
with `b1=c1=b2=c2=0` in the trace-zero matrix chart `[-1,b;c,1]`.

## Consequence for cyclic `[60]`

Keeping `B5=q0` is a dead end for obtaining a geometrically simple example
near this seed. A transverse search must allow the order-5 quadratic `B5`
and the order-2 support `q0` to separate.

## Reproduction and memory

```text
python3 code/hlp_z60_linked_second_order.py
```

The exact `Fraction` calculation uses only matrices of size at most `23 x 24`.
Measured peak resident memory was 12,288 kB.
