# The split-complement cover of the first orbit-12 quotient

Fix

```text
s = e2(r3,r4,r5),       q = r2/r1,
```

after scaling `r1=1`.  The complementary roots are the roots of

```text
f_s,q(U) = U^3+(1+q)U^2+sU-(1+q)(q-s),
```

and the first radicand quotient is

```text
v^2 = G_s(q) = q(s-q)(2-q^2+(q+2)s).
```

The calculations below are verified by

```text
code/elkies22210_orbit12_split_cover_geometry.m.
```

## 1. Adjoining one complementary root

If `x` is one root of `f_s,q`, then

```text
q^2-(x^2+s-1)q-(x+1)(x^2+s) = 0.
```

Put

```text
d = 2q-(x^2+s-1).
```

The curve parametrizing a rational value of `q` together with one
rational complementary root is therefore

```text
R_s: d^2 = H_s(x)
             = x^4+4x^3+2(s+1)x^2+4sx+(s+1)^2,

q = (x^2+s-1+d)/2.
```

Its quartic discriminant is

```text
256(16s^4-16s^3-8s^2-24s-11),
```

so `R_s` has genus one for generic `s`.

After selecting `x`, write the remaining quadratic factor as

```text
U^2+(q+x+1)U+s+x(q+x+1).
```

Its discriminant is

```text
delta_s(q,x)
  = (q+x+1)^2-4(s+x(q+x+1))
  = q^2-2qx+2q-3x^2-2x+1-4s.
```

Thus the complete ordered splitting curve is obtained from `R_s` by

```text
w^2 = delta_s(q,x),
```

and the other two roots are

```text
(-(q+x+1)+w)/2,  (-(q+x+1)-w)/2.
```

## 2. Exact genus of each cover

Let

```text
D_s(q) = Disc_U(f_s,q)
       = (1+q)^2s^2-4s^3+4(1+q)^4(q-s)
         -27(1+q)^2(q-s)^2-18(1+q)^2s(q-s).
```

Generically `D_s` has five finite simple roots and a branch point at
infinity.  Hence `R_s -> P^1_q` is a degree-three cover with six simple
branch points.  The branch loci of `D_s` and `G_s` are generically
disjoint; up to a nonzero constant their resultant is

```text
s^3(s-1)^2(s+1)^3(s^2-s-1)(3s^2+12s+17)^2.
```

Riemann--Hurwitz then gives the following table.

| imposed conditions | equations over `R_s` | genus |
|---|---:|---:|
| one rational root | none | 1 |
| complete split cubic | `w^2=delta_s` | 4 |
| first radicand + one root | `v^2=G_s(q)` | 7 |
| first radicand + complete split | both equations | 19 |

For the last two rows, the pullback of the four branch points of
`v^2=G_s(q)` has respectively `12` and `24` points.  These give genera
`7` and `19` directly by Riemann--Hurwitz.

## 3. The known split fiber `s=59/49`

At

```text
s = 59/49,
q = 8/7,
(r3,r4,r5) = (-9/7,-5/7,-1/7),
v = 192/343,
```

the first radicand is a square and the complementary cubic splits.  The
generic branch assumptions remain valid.  Magma computes the exact
function-field genera

```text
genus(R_s)                              = 1,
genus(complete split curve)             = 4,
genus(first radicand + one root)         = 7,
genus(first radicand + complete split)   = 19.
```

There is an important arithmetic warning.  The elliptic curve `R_s` at
this fiber has

```text
rank R_s(Q) = 4,       R_s(Q)_tors = 0,
```

whereas the first-radicand elliptic quotient has

```text
rank E_s(Q) = 2,       E_s(Q)_tors = Z/4.
```

For the `S3` splitting curve of genus four, the standard group-algebra
decomposition contains two copies of the one-root elliptic quotient.
Consequently its Jacobian already has rank at least eight.  Classical
Chabauty on the direct split curve cannot work.

## 4. Recommended rigorous attack

Complete splitting forces `D_s(q)` to be a square.  It is therefore more
efficient to begin with the necessary `V4` cover

```text
A_s:  v^2 = G_s(q),    z^2 = D_s(q).
```

This curve has genus seven.  Its three quotient curves are

```text
v^2 = G_s(q)             (genus 1),
z^2 = D_s(q)             (genus 2),
y^2 = G_s(q)D_s(q)       (genus 4),
```

and

```text
Jac(A_s) ~ E_s x Jac(C_D) x Jac(C_GD)
```

over `Q`.  At `s=59/49`, `D_s(q)` is an irreducible quintic.  This gives
a concrete rigorous workflow:

1. compute 2-descent rank bounds and saturated Mordell--Weil subgroups
   for the genus-two and genus-four quotient Jacobians;
2. use the three quotient maps for a Mordell--Weil sieve on `A_s`;
3. if the total rank is below `7`, finish by classical Chabauty; if it is
   `7` or `8`, use quadratic Chabauty with the `V4` correspondences (the
   product decomposition gives Neron--Severi rank at least three);
4. factor `f_s,q` at the resulting finite list of `q` values.  A square
   discriminant can also come from an irreducible cyclic cubic, so this
   final factorization is essential.  Test the remaining three halving
   radicands only on the genuinely split values.

This discriminant-first route keeps the certification problem at genus
seven with three explicit quotient factors.  Direct work on the exact
genus-19 cover, or classical Chabauty on the rank-at-least-eight
genus-four split curve, is much less favorable.
