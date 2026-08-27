# Halving a Point on the `[2,4,8]` One-Split Family

Goal: start from the `[2,4,8]` family inside the `M_1(8,4)` tangent-cover
family and try to halve one of the additional rational `2`-torsion classes.
This would push the visible `2`-primary torsion beyond `[2,4,8]`, typically
toward `[4,4,8]` or another order-`128` group.

## Algebraic Setup

The parent family is

```text
C: y^2 = f(x) = x*A(x)*B(x),
```

with

```text
t = (2R^2 + (1-w^2)R - 2w^2)/(4(w^2-1)).
```

The known `[4,8]` tangent cover is the first halving of `T_x=[x,0]`:

```text
h(x) - x*(M*x+N)^2 = c4*(x^2+U*x+V)^2,      h=A*B.
```

The `[2,4,8]` subfamily is obtained by requiring exactly one of the two
quadratics `A,B` to split over `Q`.  A rational root `alpha` of the split
quadratic gives an extra rational Weierstrass point and hence an extra
rational `2`-torsion class

```text
T_alpha = (alpha,0) - infinity.
```

To halve this class, shift `alpha` to the origin:

```text
f(alpha+z) = z*h_alpha(z).
```

Then `T_alpha` is divisible by `2` if and only if there are rational
`U,V,M,N` such that

```text
h_alpha(z) - z*(M*z+N)^2
  = lc(h_alpha)*(z^2+U*z+V)^2.              (*)
```

This is the same tangent-square condition as for halving `[x,0]`, just applied
after translation to the chosen rational branch point.  The search script
solves this equation first and only then asks Magma to verify divisibility in
the exact Jacobian.

The split conditions themselves factor.  For example, on the reduced square
branch `w=s^2`,

```text
disc(A) ~ (R-1)^2*(s^2-R)*(s^2+R)
          *(R^2*s^4 + 2R*s^4 + R^4 - s^4 - 2R^3 - R^2)
          /((s-1)^2*(s+1)^2*(s^2+1)^2),

disc(B) ~ (R-1)^2*(R*s^2 + 3s^2 - 3R - 1)
          *(R*s^2 + 3s^2 + 3R + 1)
          /((s-1)*(s+1)*(s^2+1)),
```

up to rational square constants.

## Finite-Prime Diagnostic

I added

```text
code/m18_m14_248_halve_split_root_finite.m
```

which works over finite fields and counts good affine residues satisfying:

1. the first `[4,8]` halving,
2. exactly one of `A,B` split,
3. a split branch-point class is divisible by `2`.

The output is in

```text
data/m18_m14_248_halve_split_root_finite.txt
```

Key counts:

```text
p 7  one_split 4   root_halvable 0
p 11 one_split 4   root_halvable 0
p 13 one_split 8   root_halvable 8
p 17 one_split 24  root_halvable 24
p 19 one_split 36  root_halvable 0
p 23 one_split 84  root_halvable 0
p 29 one_split 96  root_halvable 104
p 31 one_split 152 root_halvable 0
p 37 one_split 200 root_halvable 144
p 41 one_split 280 root_halvable 224
p 43 one_split 340 root_halvable 0
```

Thus the split-root halving target has no good-open residues at `p=7` or
`p=11`, and also none at `p=19,23,31,43`.  Any rational example in this chart
must therefore reduce to bad or boundary behavior at several small primes.

## Exact Search

The exact algebraic search is

```text
code/m18_m14_248_halve_split_root_search.m
```

It enumerates rational `R,w`, applies the `[4,8]` cover filter, requires
exactly one split quadratic, solves the shifted tangent condition `(*)` for
each rational split root, and then verifies in the exact Jacobian.

The height-20 sanity run found the known one-split tangent points but no
split-root halving candidates:

```text
checked          259080
cover              2532
smooth              180
one_split            12
first_tangent         8
root_tests           16
root_candidates       0
exact_verified        0
hits                  0
```

The height-50 run is recorded in

```text
data/m18_m14_248_halve_split_root_h50.txt
```

and finished with

```text
checked          9566648
cover              15322
smooth               922
one_split             28
first_tangent          8
root_tests            16
root_candidates        0
exact_verified         0
hits                   0
```

## Conclusion

This specific way of halving a point on the `[2,4,8]` family, namely halving
the extra rational split-root `2`-torsion point, is not promising by ordinary
height search.  Algebraically it is forced onto simultaneous bad/boundary
behavior at small primes, and the exact shifted-tangent search found no
candidates through height `50`.

The other second-halving directions on the parent `[4,8]` family remain the
previously studied `[8,8]` and `[4,16]` covers.  Those are separate from
halving the new split-root class; the `[8,8]` direction remains the more
locally plausible of those two, while `[4,16]` is already forced to boundary
modulo `7`.
