# Agent notes: contact-5/contact-5 route for `Z/5 x Z/5`

Files:

```text
code/agent_z5x5_contact5_contact5_symbolic.m
code/agent_z5x5_contact5_contact5_search.m
```

I modeled this on the successful contact-5/contact-6 order-30 family, but the
equal-contact-order case collapses much more rigidly.

## Normalized setup

Let `C: y^2=f(x)` be an odd quintic model.  A rational 5-contact presentation
at `x=a` has the form

```text
f = h^2 - K*(x-a)^5,     deg(h) <= 2.
```

For two distinct rational contact points, an affine change of `x` puts them at
`0` and `1`.  Write

```text
h0 = e*x^2 + d*x + c,
h1 = A*x^2 + B*x + C,
f  = h0^2 - K*x^5 = h1^2 - K*(x-1)^5.
```

The same `K` is forced by the leading coefficient of `f`.  The coefficient
equations are

```text
E0 = c^2 - C^2 - K = 0,
E1 = 2*c*d - 2*B*C + 5*K = 0,
E2 = d^2 + 2*c*e - B^2 - 2*A*C - 10*K = 0,
E3 = 2*d*e - 2*A*B + 10*K = 0,
E4 = e^2 - A^2 - 5*K = 0.
```

## Quotient by the same-contact branch

Before normalizing the second point to `1`, say it is at `r`.  Then

```text
(h0-h1)*(h0+h1) = K*(x^5 - (x-r)^5).
```

The component `r=0` is exactly the same-contact branch: `h1 = h0` or
`h1 = -h0`, giving the same 5-torsion class up to sign and no independent
second class.

On the open `r != 0`, set `x=r*X`.  After dividing by `r^5`, the fixed
quartic is

```text
Phi(X) = X^5 - (X-1)^5
       = 5*X^4 - 10*X^3 + 10*X^2 - 5*X + 1.
```

Thus a nontrivial rational two-contact solution would require a factorization

```text
K*Phi = U*V,       deg(U), deg(V) <= 2,
```

where `U=h0-h1` and `V=h0+h1` have rational coefficients.

## Obstruction

Magma gives

```text
Factorization(Phi over Q) = Phi
Discriminant(Phi) = 125
IsIrreducible(Phi) = true
```

Over `Q(sqrt(5))`,

```text
Phi = (X^2 - X + (5-sqrt(5))/10)
      *(X^2 - X + (5+sqrt(5))/10).
```

Since `Q[x]` is a UFD and `Phi` is irreducible of degree `4`, the product
`U*V=K*Phi` is impossible with both `U,V` of degree at most `2`, unless
`K=0`.  The case `K=0` makes `f` a square/boundary, not a smooth genus-2
curve.

Conclusion: the literal rational contact-5/contact-5 route cannot produce
`Z/5 x Z/5`.  It is not just cold in a box; the normalized ansatz is empty off
the same-contact boundary.

## Finite-field smoke test

The obstruction is specifically rational factorization of `Phi`.  Over finite
fields where `Phi` has a degree-2 divisor, the construction reappears exactly.
For example, over `F_11`,

```text
U  = x^2 + 7*x + 1
V  = 5*x^2 + 10*x + 1
h0 = 3*x^2 + 3*x + 1
h1 = 2*x^2 + 7*x
f  = 10*x^5 + 9*x^4 + 7*x^3 + 4*x^2 + 6*x + 1
```

Magma reports

```text
#J(F_11) = 100
Order(D0) = 5, Order(D1) = 5, Order(D0+D1) = 5
no nontrivial relation a*D0+b*D1=0 for a,b in F_5
```

Similar independent pairs occur at `p=19`, `p=29`, `p=31`, and `p=41`.  This
confirms that the sum/difference algebra is the right model; the rational
failure comes from the irreducibility of the fixed cyclotomic quartic.

## Next useful branch

To keep a contact-style attack alive, the second 5-torsion class should not be
forced to be a rational point-minus-infinity class.  The next realistic variant
is: keep one rational 5-contact `f=h0^2-K*x^5`, but represent the second
5-torsion class by a degree-2 Mumford divisor/contact condition, or work with a
quadratic/conjugate 5-contact pair and descend the resulting Galois-invariant
Jacobian class.
