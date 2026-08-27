# Starting from M(2,4,4) toward [2,4,8]

This route starts from the model in `paper/main.tex` for the component of
`A/M(2,4,4)` where the doubled `(4,4)` subgroup is a maximal isotropic in
`J[2]`.

The genus-2 curve is

```text
C: y^2 = x(x+u^2)(x+v^2)(x^2 + (t^2+2s)x + s^2).
```

The parameters come from the elliptic fiber product

```text
E2: Y^2 = X(X+t^2)(X+4s+t^2),
E3: Y^2 = X(X^2 - 4(t^2+2s)X + 16s^2).
```

For `R in E3` and `P=(x1,y1) in E2`, put

```text
phi(R) = (Y_R^2/(4X_R^2), Y_R*(16s^2-X_R^2)/(8X_R^2)) in E2,
Q = P + phi(R) = (x2,y2),
u = y1/(2*x1),
v = y2/(2*x2).
```

This gives curves with generic torsion `[2,4,4]`.  The natural target for
`[2,4,8]` is to halve the order-4 class whose double is the quadratic
2-torsion class

```text
Tq = [x^2 + (t^2+2s)x + s^2, 0].
```

This is motivated by the known split `(2,4,8)` LMFDB example mentioned in
`paper/main.tex`, where the order-8 point quadruples to the irreducible
quadratic factor.

## Scripts

```text
code/m244_isogeny_probe.m
code/m244_to_248_sample_search.m
```

The probe confirms the elliptic quotient picture and the explicit dual
2-isogeny formula above.

The sampler has two modes:

- ordinary mode computes exact `J(Q)_tors` for sampled fiber-product points;
- `targeted_only:=true` first tests whether the half of `Tq` is itself
  divisible by 2, and only computes full torsion if that target halving occurs.

## Pilot result

The small pilot

```text
magma -b height:=4 point_bound:=80 max_hits:=5 max_curves:=200 \
    code/m244_to_248_sample_search.m \
    > data/m244_to_248_sample_h4_b80.txt
```

validated the construction:

```text
tested_curves 201
hits 0
TORSION_COUNTS
  [2,2,4,4] 10
  [2,4,4]   190
```

Two broader height-6 pilots were started and interrupted after producing no
hits:

```text
data/m244_to_248_sample_h6_b120.txt
data/m244_to_248_targeted_h6_b120.txt
```

## Conclusion

Yes, this is a sensible route.  It is less restrictive than starting from
`A(2,2,4,4)`, and it has a clean geometric parameterization by an elliptic
fiber product.  The first pilots confirm that the model produces the expected
`[2,4,4]` curves, but they did not immediately find `[2,4,8]`.

The next useful step is not a broader blind sample of fiber-product points.  It
is to derive finite-prime and then symbolic conditions for the divisibility of
the quadratic-half 4-torsion class, analogous to the reduced `[8,8]` work in the
`M_1(8,4)` family.
