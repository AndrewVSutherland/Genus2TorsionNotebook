# Exact rational-root projection on the `c=2` slice

This was derived before the search pivoted to the globally stronger constant
section `d=5/3`.  It records the useful exact geometry and the negative CRT
test already completed on `c=2`.

Let `rho=-b/(2b-1)` be the root of the rational linear factor of the compact
M(12) quintic.  For `x=d Z-2`, the Hermite quotient `S` has rational root

```
u=(rho+2)/d.
```

At this root `Q(u)=v^2`, `A(u)=v^5`.  Using the line of slope `k` through the
known conic point `(0,s)` gives

```
u=(t^2-s^2-1-2 s k)/(k^2-1),   v=s+k u,
a=(v^5-A_0(u))/(u^2(u-1)^2).
```

If `C_j` is the `(Z-u)^j` Taylor coefficient of `S`, comparison at the M(12)
linear-factor root gives

```
R=C_1/lc(S)=(rho(rho+1)/d^2)^2.
```

Write `h=+/-sqrt(R)`.  Since `rho=d u-2` on `c=2`, one obtains the second
square cover and the complete recovery of `d`:

```
ell^2=u^2+8h,
d=(3u+/-ell)/(2(u^2-h)).
```

Then `b=rho/(2rho+1)`.  The normalized second and third Taylor coefficients
recover `w` linearly, after which the full six-coefficient identity and the
square multiplier are checked exactly.  Thus the `c=2` curve can be searched
in only the intrinsic variables `(s,t,k)` followed by two rational-square
tests; no Groebner basis is needed.

Before this projection was found, direct ambient CRT used the seven primes
`7,11,13,17,19,23,29`.  It combined 1,920 local choices, found three tuples
whose six ambient coordinates all had canonical rational reconstructions,
and found zero exact identities.  Every six-prime subset was subsequently
superseded by the `d=5/3` priority.  Accordingly this is a meaningful negative
bounded search, not a proof that the `c=2` curve has no rational point.
