# Internet sweep for additional infinite torsion families

Date: 2026-07-02.

Goal: find positive-dimensional genus-2 Jacobian torsion families not already
recorded locally, and keep only families that look useful for simple Jacobians.

## Actionable new families

### Order 11

Source: K. Daowsud and T. A. Schmidt, *Continued fractions for rational
torsion*, arXiv:1708.05511 / J. Number Theory 189 (2018), 115--130.
The paper also quotes Flynn's older order-11 family.

Flynn family:

```text
C_t: y^2 = x^6 + 2*x^5 + (2*t+3)*x^4 + 2*x^3
             + (t^2+1)*x^2 + 2*t*(1-t)*x + t^2.
```

Daowsud-Schmidt family:

```text
C_u: y^2 = x^6 - 4*x^5 + 8*(1+u)*x^4 - (10+32*u)*x^3
             + 8*(1+6*u+2*u^2)*x^2
             - 4*(1+6*u+16*u^2)*x + 64*u^2 + 1.
```

Local checks:

```text
Flynn t=1:
  f = x^6 + 2*x^5 + 5*x^4 + 2*x^3 + 2*x^2 + 1
  Magma TorsionSubgroup(J)(Q) = [11]
  irreducible Lp at p=3: 9*x^4 + 6*x^3 + 4*x^2 + 2*x + 1
  Sage geometric_endomorphism_algebra_is_field(B=100) = True
  Sage geometric_endomorphism_ring_is_ZZ(B=100) = True

Daowsud-Schmidt u=1:
  f = x^6 - 4*x^5 + 16*x^4 - 42*x^3 + 72*x^2 - 92*x + 65
  Magma TorsionSubgroup(J)(Q) = [11]
  irreducible Lp at p=3: 9*x^4 + 6*x^3 + 4*x^2 + 2*x + 1
  Sage geometric_endomorphism_algebra_is_field(B=100) = True
  Sage geometric_endomorphism_ring_is_ZZ(B=100) = True
```

Conclusion: `[11]` should be added to the simple-family inventory.  This is a
new odd-prime base for attempts at `[22]`, `[44]`, `[55]`, or `[66]`, but the
first realistic target is probably `[22]` by forcing a rational Weierstrass
point.

### Order 23

Source: H. Kuru and M. Sadek, *Quadratic torsion orders on Jacobian varieties*,
arXiv:2410.14455.  The abstract gives a one-parameter family with order
`2g^2+7g+1`; for `g=2` this is `23`.

The genus-2 specialization of the corollary can be written as follows.  For
`t != 0, +/-1`, set

```text
beta  = (t^2 + 1)^2/(4*t^2)
sbeta = (t^2 + 1)/(2*t)
s      = (t^2 - 1)/(2*t)
alpha = beta - s^5/(beta*sbeta)
lambda = (alpha - 1)^4/((alpha - beta)^2*alpha)
```

and

```text
A(x) = (x^3*(x-alpha)^2
        - (x-1)*((x-1)^4 - lambda*(x-beta)^2*x))
       /(2*(x-alpha)*(x-beta)),
C_t: y^2 = A(x)^2 - lambda*x^4*(x-1).
```

The TeX source line for `alpha` appears to omit the exponent on `beta`; using
`beta^((g+1)/2)` is the interpretation that makes `A(x)` polynomial and exactly
reproduces the published `t=2` example.

Published `t=2` primitive model:

```text
y^2 = -299054816676000*x^5
      + 937313042871529*x^4
      - 1165161421194050*x^3
      + 677279473485625*x^2
      - 132825168000000*x
      + 8294400000000.
```

Local checks:

```text
Magma TorsionSubgroup(J)(Q) = [23]
irreducible Lp at p=13: 169*x^4 - 52*x^3 + 24*x^2 - 4*x + 1
Sage geometric_endomorphism_algebra_is_field(B=100) = True
Sage geometric_endomorphism_ring_is_ZZ(B=100) = True
```

Conclusion: this is the best new find.  It is a genuine one-parameter simple
Jacobian family with a large prime torsion point.  The obvious next algebraic
move is to force a finite rational branch point in this family and see whether
an infinite `[46]` subfamily or at least a simple `[46]` example exists.

## Follow-up: order 21 extracted

The 1991 Leprévost scan has now been transcribed.  It gives a direct
one-parameter family with generic exact torsion `[21]`; the local `t=1` fiber
has full torsion `[21]` and geometric endomorphism ring `ZZ`.  The formulas,
bad-fiber factorization, and exact checks are in
`notes/z21_leprevost_and_hlp_2026_07_18.md` and
`code/z21_leprevost_family_verify.m`.

## References still not equation-extracted

The Kuru-Sadek introduction lists older genus-2 infinite families with rational
torsion orders `13`, `15`, `17`, `19`, `22`, `23`, `24`, `25`, `26`,
`27`, and `29`, citing Leprévost's 1991 and 1995 papers.  The 1995 paper title
explicitly includes "torsion et simplicité", so it is highly relevant for our
simple-Jacobian filter.  Their remaining machine-usable equations still need
to be extracted from the original sources.

## Excluded or lower priority

- Howe's order-48 family is constructed by gluing elliptic curves, so the
  resulting Jacobians are isogenous to products and are not appropriate for the
  simple-Jacobian-only goal.
- Howe-Leprévost-Poonen split-Jacobian constructions are similarly useful for
  large torsion examples but not for simple Jacobians.
- Elkies' `N=32` web page remains interesting, but the page omits the
  parametrizing formula; locally we only have the printed example and halving
  equations.
