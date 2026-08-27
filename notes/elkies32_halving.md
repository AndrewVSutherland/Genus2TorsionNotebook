# Elkies N=32 Example: Halving Test

Elkies' page announces a one-parameter family with a rational torsion point of
order 32, but does not give the family formula.  The only usable input from the
page is the printed member

```
y^2 = (15*x - 1)*(1056*x^4 + 156183*x^3 + 26297*x^2 + 649*x - 121).
```

It has Weierstrass points at infinity and `x=1/15`, and non-Weierstrass
rational points `(0,11)` and `(-1,1440)`.  The difference between either
non-Weierstrass point and either Weierstrass point has order 32.

## Algebraic halving condition

For an odd quintic model `C: y^2 = F(X)` and a point `P=(0,y0)`, a divisor
`P - infinity` is divisible by 2 if and only if there are

```
u(X)   = X^2 + a*X + b,
ell(X) = m*X^2 + n*X + k,    k = +/- y0,
```

with

```
F(X) - ell(X)^2 = lc(F) * X * u(X)^2.
```

Writing `F(X)=F0+F1*X+...+F5*X^5`, the top two coefficients eliminate `a,b`:

```
a = (F4 - m^2)/(2*F5),
b = (F3 - 2*m*n - F5*a^2)/(2*F5).
```

The remaining conditions are the two plane equations

```
K2 = F2 - n^2 - 2*m*k - 2*F5*a*b = 0,
K1 = F1 - 2*n*k - F5*b^2 = 0.
```

For the finite Weierstrass point `w=1/15`, I used

```
X = 1/(x-w),   Y = y*X^3,
g(X) = X^6*f(w+1/X),
```

so that `w` becomes the point at infinity, and then applied the same test.

## Computation

Script:

```
code/elkies32_halving_conditions.m
```

Output:

```
data/elkies32_halving_conditions.txt
```

The script tests all four order-32 divisor classes visible on the printed
curve:

```
(0,11) - infinity
(-1,1440) - infinity
(0,11) - (1/15,0)
(-1,1440) - (1/15,0)
```

For both signs `k=+/-y0`, every resulting ideal `<K1,K2>` is zero-dimensional
and has no rational affine solution in `(m,n)`.  The resultant in `n` is an
irreducible degree-16 polynomial over `Q` in each case.

There is also an immediate local obstruction at the good prime `p=7`: for all
four divisor classes and both signs, the halving equations have zero solutions
mod 7.  This proves that none of these four marked order-32 points on Elkies'
printed curve is divisible by 2 over `Q`.

## Conclusion

The printed `N=32` member does not produce order 64 by halving any of its four
obvious order-32 divisor classes.  The obstruction is not a height-search
failure: it is already visible in the exact halving equations and modulo the
good prime `7`.

To search the actual infinite `N=32` family for a halved specialization, the
missing ingredient is Elkies' omitted parametrization.  Without that formula,
we can only test individual printed or reconstructed members, not perform a
meaningful family-level search.
