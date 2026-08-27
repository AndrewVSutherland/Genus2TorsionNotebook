# The split-quadratic `d=5/3` curve

## Outcome

No rational cyclic-order-60 example was found on this slice.  The slice is
nevertheless a genuine, everywhere-tested one-dimensional candidate: it has
smooth open points at every prime tested, and it is the only constant value of
`d` of rational height at most 100 surviving the independent hyperplane audit.

The strongest negative result is an exact 7-adic graph search.  It excludes a
cyclic-order-60 point on the unique good/open branch above `p=7` whenever

* `c` has height at most 1000 and is 7-integral;
* the other affine coordinates `b,w,s,t,a` are 7-integral and have reduced
  numerator and denominator at most 9,787,268,040.

This is not a global nonexistence proof: points of larger height, points bad or
nonintegral at 7, and other nonconstant sections of the full surface remain.

## Exact projection

Use `x=d Z-c`, with `d=5/3`, and write the two roots of the full quadratic
Padé denominator as `Z=0,1`.  The Hermite construction is

```
Q = Z^2 + (t^2-s^2-1) Z + s^2,
S = (A^2-Q^5)/(Z^2 (Z-1)^2),
```

where `A` is monic of degree five, has the prescribed values and derivatives
at 0 and 1, and has free `Z^4` coefficient `a`.  Exact matching requires

```
S(Z) = tau F_{b,w}(d Z-c),                 tau in Q^{*2}.
```

Let `rho=-b/(2b-1)` be the rational root of the linear factor of the compact
M(12) quintic.  Then `u=(rho+c)/d` is a rational root of `S`.  At that root,
`Q(u)=v^2` and `A(u)=v^5`.  The conic `v^2=Q(u)` is parametrized by the slope
`k` through `(0,s)`:

```
u = (t^2-s^2-1-2 s k)/(k^2-1),
v = s+k u,
a = (v^5-A_0(u))/(u^2 (u-1)^2).
```

Thus `(s,t,k)` are intrinsic coordinates before the remaining M(12)
conditions are imposed.

Expand the M(12) quintic at its linear-factor root.  If `f_j` is the
coefficient of `(x-rho)^j`, then

```
f_1/f_5 = rho^2 (rho+1)^2.
```

If `C_j` is the coefficient of `(Z-u)^j` in `S`, comparison after
`x=rho+d(Z-u)` gives the exact square cover

```
R = C_1/lc(S) = (rho(rho+1)/d^2)^2.
```

For each of the two signs `h=+/-sqrt(R)`, fixed `d` gives

```
m^2 = 1+4 h d^2,       rho=(-1+/-m)/2,
b = rho/(2 rho+1),     c=d u-rho.
```

The normalized second and third Taylor coefficients recover `w` linearly
(the fourth handles the exceptional denominator).  The implementation then
checks all six coefficients of `S=tau F` over `Q`, including that `tau` is a
nonzero rational square.  The recovery formulas were independently checked
against every local point through `p=31`: 35 of 35 ambient points recovered.

## Exact elimination audit

An independent exact elimination audit used five proportionality minors. In
the raw affine chart their `(total degree, term count)` pairs are `(12,158)`,
`(12,172)`, `(12,183)`, `(9,116)`, and `(8,64)`. The last minor is linear in
`c`; its leading coefficient has degree 4 with 8 terms and its constant
coefficient has degree 8 with 56 terms. Eliminating `c` modulo 7 gives
successive resultant sizes `(44,45535)`, `(36,15975)`, `(28,5120)`, and
`(20,983)`. Their only common boundary factor is `b(2b-1)lc(S)`. After
saturation, the smallest open factor is irreducible of degree 16 with 427
terms. This explains why direct characteristic-zero elimination is unattractive
and motivates the intrinsic Taylor projection and 7-adic graph search.

## Local data and CRT

Smooth/open local counts, after deduplicating intrinsic `(s,t,k)`, are

| p | ambient points | intrinsic triples |
|---:|---:|---:|
| 7 | 1 | 1 |
| 11 | 2 | 2 |
| 13 | 1 | 1 |
| 17 | 3 | 3 |
| 19 | 3 | 3 |
| 23 | 8 | 7 |
| 29 | 10 | 10 |
| 31 | 7 | 7 |

The unique `p=7` point has

```
(b,w,c,d,s,t,a,lambda)=(3,5,6,4,1,4,5,2),
F=3 Z^5+2 Z^4+3 Z^3+5 Z^2+2 Z+2,
q=Z^2+5 Z+3,
v=5 Z+5,
B=2 Z^2+2 Z+3,
A=Z^5+Z^4+Z^3+2 Z+4.
```

Direct arithmetic verifies `A^2-B^2 F=q^5`. Cantor arithmetic gives
`ord(D_5)=5`, `ord(D_12)=12`, and `ord(D_5+D_12)=60`; moreover
`#J(F_7)=60`. Thus the local point has exactly the intended cyclic-order-60
control, not merely the proportionality equations.

CRT reconstruction was performed in `(s,t,k)`, rather than in the much larger
ambient coordinates.  The eight-prime product and every seven-prime subset
were used, so one prime may be omitted.  Across 40,152 local combinations,
2,367 coordinate triples reconstructed and 2,366 were distinct.  Exact
Taylor recovery and the full polynomial identity gave zero hits.

A diagnostic rational-height sieve found no triple of simultaneous good/open
reductions with `height(s),height(t),height(k) <= 30`.  This diagnostic is not
claimed as a rigorous exclusion of points bad at one of the small primes.

## 7-adic graph search

Modulo 7 the unique point is

```
(b,w,c,d,s,t,a) = (3,5,6,4,1,4,5),   tau=4.
```

The residual Jacobian in `(b,w,s,t,a)` is

```
[2 3 2 5 6]
[2 3 3 4 5]
[6 2 2 1 4]
[5 3 2 6 1]
[1 3 3 3 1]
```

and has determinant 2 modulo 7.  Hence `c` is an etale parameter.  For every
reduced rational `c=n/e` with `|n|,e <= 1000`, `7` not dividing `e`, and
`c=6 (mod 7)`, the other five coordinates were Newton-lifted to `7^24`.

There are 152,249 such values of `c`.  At modulus

```
7^24 = 191581231380566414401,
```

symmetric rational reconstruction is unique for
`|numerator|,denominator <= floor(sqrt((7^24-1)/2)) = 9,787,268,040`.
There were 6,428 simultaneous five-coordinate reconstructions; every one
failed the complete exact identity.  Therefore there are zero exact hits in
the good/open 7-adic branch within the bounds stated above.

## Reproduction

```
python3 code/m12_general5_fullquad_d53_search.py --mode crt \
  --primes 7,11,13,17,19,23,29,31

python3 code/m12_general5_fullquad_d53_search.py --mode height --height 30 \
  --primes 7,11,13,17,19,23,29,31

python3 code/m12_general5_fullquad_d53_padic.py --height 1000 --precision 24
```

Both programs use elementary polynomial arithmetic and constant/small memory;
all long runs were made with a 1.9 GB virtual-memory cap.
