# M(2,2,2,6): trying to halve an order-6 class in `J`

This is the less restrictive route suggested after the theta-doubling tests: do not ask for the half to lie on the curve.  Instead, look for an order-6 class in the coset

```text
g + J[2],    g = P - infinity,
```

which lies in `2J(Q)`.  A half would have order 12.

The odd quintic model is again

```text
Y^2 = L1 L2 L3 L4 L5,    P = (0,2),
```

with the five finite rational Weierstrass points `W_i` cut out by `L_i=0`.  The tested classes are

```text
g,
g + (W_i - infinity),
g + (W_i + W_j - 2 infinity),  1 <= i < j <= 5.
```

The exceptional order-3 class is automatically discarded by checking `Order(D) eq 6`.

## Exact rational specialization search

`code/m2226_order6_halving_direct_search.m` specializes `[s:m:n]`, constructs the hyperelliptic Jacobian, forms every class `g+T`, checks it has order 6, and then calls Magma's exact

```magma
IsDivisibleBy(D, 2)
```

over `Q`.

Completed run:

```text
magma -b height:=8 max_hits:=20 progress_interval:=500 code/m2226_order6_halving_direct_search.m
```

Result:

```text
checked curves 1543
checked classes 24688
hits 0
```

A larger `height:=12` run was stopped after more than 73,000 class checks, also with no hits.

## Exact finite-field sieve

`code/m2226_order6_halving_finite_field_sieve.m` performs the same test over finite fields, using `AbelianGroup(J)` to solve divisibility by 2 in `J(F_p)`.  Magma's `IsDivisibleBy` is only available over `Q`, so over finite fields the script writes `J(F_p)` as an explicit finite abelian group and checks whether the coordinate vector lies in `2G`.

The first useful prime is `p=7`; `p=3,5` have no good points for this model.  At `p=7`:

```text
good parameter points in P^2(F_7): 10
order-6 classes tested: 150
divisible by 2: 0
```

Thus there is no example with good reduction modulo 7.  As with the 2-torsion halving route, a complete global proof would still require checking rational points that reduce to the boundary modulo 7.  But the exact rational search and the exact good-reduction finite-field test both indicate that this general order-6 halving route does not produce nondegenerate `(2,2,2,12)` examples inside the split `M(2,2,2,6)` family.
