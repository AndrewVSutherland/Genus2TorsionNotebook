# Nonconstant-`B` projection and normalization for cyclic order `49`

## Rational root normalization

On the chart `u*v != 0`, put

```text
w=u^2,  tau=v/u,  B=u*(1+tau*x),  q=-1/tau.
```

At the rational zero `x=q` of `B`, the norm identity gives

```text
A(q)^2=(q-1)*(q-r)^7.
```

Consequently `(q-r)/(q-1)` is a rational square.  Writing it as `t^2`
gives the exact normalization

```text
r=q-(q-1)*t^2,
A(q)=sign*(q-1)^4*t^7,
A'(q)=sign*(q-1)^3*t^5*(t^2+7)/2.
```

The other choice of sign is obtained by `t -> -t`, so it does not require
a separate projection.

In `X=x-q` coordinates, the `X^2` and `X^3` coefficients solve the next
two Taylor coefficients of the monic quartic `A`.  The only residual
equations are the coefficients `X^4,...,X^7`.  This removes the square
condition and the sign involution before elimination, without a Groebner
basis.  The derivation and exact identity check are in
`code/z49_nonconstant_root_normalization.py`.

## Exact `T=1` slice

For `T=tau=1`, hence `q=-1` and `r=2*t^2-1`, the four normalized equations
in `(a,b,w,t)` have profiles

```text
X^4: degrees (4,4,2,28), 39 terms
X^5: degrees (4,4,2,30), 62 terms
X^6: degrees (4,4,2,32), 90 terms
X^7: degrees (2,2,1,16), 19 terms.
```

The `X^7` equation is linear in `w`.  On its regular chart, solving for
`w` and stripping the open-boundary factors `t=0` and `t=1` leaves three
quartics with `t`-degrees `23,25,27`.  The two staged `a`-resultants factor
as follows, with tuples `(degree in b, degree in t, multiplicity)`:

```text
R_45: [(0,1,18),(0,1,10),(0,2,8),(2,14,1),(10,78,1)]
R_46: [(0,1,16),(0,1,10),(0,2,8),(12,102,1)].
```

Eliminating `b` between the `b`-degree `2` and `b`-degree `12` factors
gives

```text
t^40 * P_248(t),
```

where `P_248` is irreducible over `Q`.  The only rational linear root is
`t=0`, for which `r=q` and `B(r)=0`.  Thus this small regular component has
no open rational point.  The `b`-degree `10` versus `b`-degree `12`
resultant reached the imposed `300 MB` cap, so that large component and the
exceptional denominator chart remain unresolved.  This is deliberately
reported as a partial exact exclusion, not as an exclusion of the full
`T=1` slice.

## The `q=infinity` / integral `3`-adic chart

In slope coordinates the square relation is

```text
t^2=(1+tau*r)/(1+tau).
```

The branch at `q=infinity` has `t=1+tau*s`, giving the polynomial formula

```text
r=1+2*s+tau*(s^2+2*s)+tau^2*s^2.
```

The only integral open disk at `3` has

```text
a=1+9*A, b=1+3*B, w=1+9*W, tau=9*L, s=4+9*S.
```

After division by the common factor `9`, the four initial forms modulo `3`
are

```text
A+B^2+B+1,
A+B^2+L-W-1,
-A-B^2+L-S+W-1,
A+B^2-B-L-S+W-1.
```

They cut out one exceptional line:

```text
A=-(B^2+B+1), L=-1, W=-B, S=-1.
```

One more strict transform again gives four independent linear equations in
the four transverse variables.  Continuing once more yields the adapted
coordinates used by the search:

```text
a   = -8-36*B-36*B^2-27*B^3+81*A2,
b   = 1+3*B,
w   = 28+18*B-27*B^2-27*B^3+81*W2,
tau = 18+27*B+81*L2,
s   = 22-27*B-27*B^2+81*S2.
```

Their initial forms are

```text
A2-B^4-B^3-B^2+B-1,
A2-B^4-B^3+L2-W2-1,
-A2+B^4+B^3-B^2+B+L2-S2+W2,
A2-B^4-B^3-B^2+B-L2-S2+W2.
```

The Jacobian in `(A2,W2,L2,S2)` is invertible modulo `3`.  Hence the
nonconstant integral locus is a single smooth `3`-adic disk parametrized by
`B`; it is not a collection of hidden boundary components.

The coordinate-aware search enumerated all `5,927` reduced rational `B`
of height at most `80` with denominator prime to `3`.  It lifted the four
adapted coordinates to `3^17`, reconstructed with bound `5000`, and exact
tested the quotient equations and the additional requirement that `w` be
a rational square:

```text
lifts                 5,927
full reconstructions     10
exact quotient points     0
rational-square scales    0
exact open hits            0
```

This search used `84,720 KB` peak RSS and `24.1` seconds.  It is a bounded
negative result, but it is targeted at the only integral `3`-adic disk and
uses coordinates adapted to its normalization.  The complementary direct
root-parameter search in `code/z49_structural_root_parameter_hensel.py`
tested `3,095` rational root parameters of height `50`, with `10,816`
transverse lifts. Sixty cases reconstructed all five coordinates within
bound `5000`, but none passed exact substitution. The result is recorded in
`data/z49_structural_root_parameter_h50.txt`.

## Verdict

The nonconstant incidence curve has a clean rational square normalization
and a unique smooth integral `3`-adic disk, which supports the point-count
evidence for one curve-sized component.  No rational point was found.  The
new exact information is the exclusion of the small regular component of
the genuinely nonconstant `T=1` slice; the main unresolved algebraic target
is its large `(b-degree 10, b-degree 12)` projection component, preferably
handled modularly rather than by another characteristic-zero resultant.

## Artifacts

```text
code/z49_nonconstant_root_normalization.py
code/z49_nonconstant_t1_projection.py
code/z49_nonconstant_infinity3_search.py
code/z49_structural_root_parameter_hensel.py
data/z49_nonconstant_t1_projection.txt
data/z49_nonconstant_infinity3_h80.txt
data/z49_structural_root_parameter_h50.txt
```
