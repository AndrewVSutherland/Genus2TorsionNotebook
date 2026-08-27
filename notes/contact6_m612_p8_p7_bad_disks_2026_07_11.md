# The four forced `p=7` disks on the `P8` extra-3 cover (2026-07-11)

## Result

The four projective parameter residues left by the good-reduction sieve split
into two very different pairs.

```text
tau = 2,-3 mod 7     double zeros of the P8 conic coordinate t;
tau = 1,-1 mod 7     parameter-pole disks, tau^2+6 = 0 mod 7.
```

On the primitive weighted endpoint signature

```text
E3: (v(e),v(L),v(U),v(v),v(N),v(R))
    = (2,-1,-2,-2,-2,-4),
```

the zero disks have no point, while each punctured pole disk has two signed
smooth branches.  An explicit Weil-pairing calculation shows that all these
pole-disk branches are orthogonal to the marked rational 3-line.  Thus they
belong to the degree-12 orthogonal support factor and to its connected
degree-24 signed lift `L^2=M`; they are not branches from the irrelevant
degree-27 nonorthogonal factor.

Consequently `p=7` does **not** obstruct the full `P8` route: the two pole
disks contain genuine local branches of the desired relative 3-cover.

For the zero disks, the scaled E3 signature is excluded at every depth by
the same nonsquare leading unit.  The other signatures in the existing
endpoint fan through `v(e)=8` are also excluded.  Valuation signatures beyond
that bounded fan have not been exhausted, so this is not a proof that the
entire zero disks are empty.

## Exact endpoint parameter

Write

```text
A = tau^2+tau-6,
B = tau^2+6.
```

The endpoint parameter is exactly

```text
e = -400*A^2*B^2 /
    (768*A^4-1200*A^2*B^2+1250*B^4).
```

Put `tau=tau0+7*s` and divide by `7^2`.  Exact reduction gives

```text
tau0 =  2:   e/7^2 = 3*s^2                       mod 7,
tau0 = -3:   e/7^2 = 6*s^2                       mod 7,
tau0 =  1:   e/7^2 = 2*(1+2*s)^2                 mod 7,
tau0 = -1:   e/7^2 = 4*(1-2*s)^2                 mod 7.
```

The nonzero squares modulo 7 are `{1,2,4}`.  Hence the first two leading
units are always nonsquares, whereas the last two are always squares away
from their unique deeper subdisks.

The deeper pole residues are

```text
tau = 1+7*3 = 22 mod 49,
tau = -1+7*4 = 27 mod 49.
```

They are the first Hensel approximations to the two roots
`alpha^2+6=0` in `Q_7`.  Near either root,

```text
e = (tau-alpha)^2 * unit,
```

and that analytic unit has square residue.  Thus every point of a pole disk
except the boundary point `tau=alpha` reaches the same E3 initial system at
some finite scaled depth.

In contrast, `tau=2,-3` are exact rational zeros of `A`.  For every `n>=1`
and unit `s`, substitution `tau=tau0+7^n*s` gives respectively

```text
e/7^(2*n) = 3*s^2 or 6*s^2 mod 7.
```

The leading unit remains nonsquare at every depth.

## Exact E3 system

Retain the two cubic coefficients and write

```text
q = x^2+U*x+v^2,
H = x^3+N*x^2+R*x+v^3,
H^2-q^3 = L^2*f.
```

For `e=7^2*E`, make the exact substitution

```text
L=l/7, U=u/7^2, v=w/7^2, N=n/7^2, R=r/7^4.
```

The diagnostic constructs the five original coefficient equations over
`Q`, applies these substitutions and the exact equation-dependent powers of
7, clears only powers of the nonzero unit `E`, and asserts that their
reductions are the following E3 forms:

```text
2*n-3*u-6*l^2 = 0,
E*(n^2+2*r-3*u^2-3*w^2)-2*l^2 = 0,
2*w^3+2*n*r-u^3-6*u*w^2 = 0,
r^2+2*n*w^3-3*u^2*w^2-3*w^4 = 0,
2*r*w^3-3*u*w^4 = 0.
```

After imposing `E*l*w*(u^2-4*w^2) != 0`, solutions exist precisely for

```text
E = 1,2,4.
```

For each such `E` there is one support and two signs of `l`:

```text
E=1: (l,u,w,n,r)=(+/-1,3,1,4,1),
E=2: (l,u,w,n,r)=(+/-2,5,4,2,2),
E=4: (l,u,w,n,r)=(+/-3,6,2,1,4),
```

with the signs interpreted modulo 7.  At every point the `5 x 5` Jacobian
with respect to `(l,u,w,n,r)` has rank 5.  Thus, after fixing any 7-adic
parameter with the indicated leading `E`, multivariate Hensel gives a unique
lift of each sign.

Composing with the four parameter disks gives:

```text
zero disks:  six primitive s-residues, no E3 point;
pole disks:  six noncentral s-residues, two rank-5 signed points each;
```

and one deeper residue in each disk, as described above.

## Separating the orthogonal support

The degree-12 resolvent uses `M=L^2` and therefore parametrizes supports
`{Q,-Q}`.  The signed rational-3 cover is degree 24; modular lex recovery
shows that `M` is nonsquare in the degree-12 support field over `F_7(u)`.

To determine which 40-support orbit contains an E3 branch, let

```text
P = [(x-1)^2, h6 mod (x-1)^2]
```

be the marked order-3 class, and let

```text
Q = [q,(H/L) mod q]
```

be a cubic-contact class.  The explicit pairing test is

```text
e_3(P,Q) =
  Res(q,H-L*h6) /
  (L^2*Res((x-1)^2,L*h6-H)).
```

The formula was independently checked against Magma's `WeilPairing` on a
smooth `P8` fiber over `F_31`: 20 signed contact classes were tested, with
2 orthogonal and 18 nonorthogonal, and the formula classified all 20
correctly.  The observed pairing values were all three cube roots of unity.

Under the E3 weighted substitution, the pairing numerator and denominator
have the same initial form

```text
E^2*l^2*w^6.
```

It is nonzero on the open E3 points, so every lifted pole-disk branch has
pairing congruent to 1 modulo 7.  The three cube roots of unity in `Q_7`
have distinct reductions; therefore the exact pairing is 1.  This certifies
membership in the degree-12 orthogonal support factor.  Since `l` is already
a nonzero residue and both signs have full rank, the signed degree-24 lift
also has `Q_7` branches.

## What is and is not closed

Certified:

- the primitive E3 layer is empty in both zero disks;
- every scaled E3 layer in a zero disk is empty, since its `E` remains
  nonsquare;
- every punctured pole disk has two signed, orthogonal, smooth E3 branches;
- the branches lie on the target degree-12 support orbit and degree-24 sign
  lift, not merely on the ambient degree-40 quotient.

For the zero disks, the previously enumerated endpoint signatures through
`v(e)=8` add no point: the `v(e)=4` and `8` signatures are empty, while the
surviving `v(e)=6` signatures require square leading `E`.  Still unresolved
are possible new valuation signatures beyond that bounded fan.  No global
rational point is produced by the pole-disk local branches; globalizing the
degree-24 cover remains the arithmetic problem.

## Reproduction

From `torsion_jac`:

```bash
magma -b code/contact6_m612_p8_weil_pairing.m
magma -b code/contact6_m612_p8_p7_bad_disks.m
```

Both runs finish in a few seconds and contain assertions for the exact disk
expansions, derivation of the E3 initial equations, complete torus-point
enumeration, fixed-parameter Jacobian ranks, and the initial Weil pairing.
