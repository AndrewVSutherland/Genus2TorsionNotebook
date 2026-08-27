# Split-`B` Hermite search for cyclic order `60`

## Family and exact marked `5`-class

Normalize the quadratic norm factor to

\[
B=x(x-1),\qquad q=x^2+(t^2-s^2-1)x+s^2.
\]

Let `A` be the monic quintic with free `x^4` coefficient `m` and endpoint
Hermite data

\[
 A(0)=s^5,\quad A'(0)=\frac52s^3q'(0),\qquad
 A(1)=t^5,\quad A'(1)=\frac52t^3q'(1).
\]

Then

\[
 F_{s,t,m}=\frac{A^2-q^5}{x^2(x-1)^2}
\]

is a quintic.  On

```text
s*t*Disc(q)*Disc(F)*Res(q,F) != 0
```

the Mumford class `[q,A mod q]` is a nonzero class of exact order `5`.
The formulas are implemented independently in
`code/hermite5_order60_modular_sieve.c` and
`code/m12_general5_fullquad_rootpair.py`.

The exact coprimality boundary is cheaper than a direct resultant suggests.
Since `A^2-B^2 F=q^5` and `Res(q,B)=s^2t^2`, on `s*t != 0` one has

\[
 \operatorname{Res}(q,F)=0\quad\Longleftrightarrow\quad
 \operatorname{Res}(q,A)=0.
\]

Symbolic reconstruction gives

\[
 \operatorname{Res}(q,A)=-\frac{s^4t^4}{4}R(s,t,m),
\]

where

```text
R = -4*m^2 - 2*m*s^5 + 6*m*s^4*t - 4*m*s^3*t^2 + 12*m*s^3
    - 4*m*s^2*t^3 - 12*m*s^2*t - 4*m*s^2 + 6*m*s*t^4
    - 12*m*s*t^2 + 6*m*s - 2*m*t^5 + 12*m*t^3 + 4*m*t^2
    + 6*m*t - 20*m + s^8 - 3*s^7*t - 2*s^7 - 2*s^6*t^2
    + 6*s^6*t - 2*s^6 + 19*s^5*t^3 - 2*s^5*t^2 - 3*s^5*t
    + 6*s^5 - 30*s^4*t^4 - 10*s^4*t^3 + 30*s^4*t^2 - 15*s^4
    + 19*s^3*t^5 + 10*s^3*t^4 - 50*s^3*t^3 - 20*s^3*t^2
    + 15*s^3*t + 26*s^3 - 2*s^2*t^6 + 2*s^2*t^5 + 30*s^2*t^4
    - 18*s^2*t - 12*s^2 - 3*s*t^7 - 6*s*t^6 - 3*s*t^5
    + 30*s*t^4 + 15*s*t^3 - 42*s*t^2 - 9*s*t + 18*s
    + t^8 + 2*t^7 - 2*t^6 - 16*t^5 - 15*t^4 + 34*t^3
    + 8*t^2 + 12*t - 24.
```

## Modular filter

At every good prime `p` outside `{2,3,5}`, cyclic rational order `60`
requires the reduction to contain an element of order `12`.  Thus

1. `60 | #J(F_p)`, and
2. the `2`-Sylow is not elementary abelian.

For a squarefree odd quintic, if `k` is its number of irreducible factors,
then `dim J(F_p)[2]=k-1`.  Therefore the second condition is exactly

```text
v_2(#J(F_p)) > k-1.
```

The C sieve obtains `k` without factoring: while computing `#C(F_p)` and
`#C(F_(p^2))`, it also counts the roots `r1,r2` of `F` in those two fields.
Then

```text
nquad = (r2-r1)/2
k = r1 + nquad + (5-r1-2*nquad > 0).
```

Independent Magma counts validate the refinement:

```text
p     good Hermite residues     #J divisible 60     also has order 4
7              82                       10                    6
11            560                       56                   50
13           1056                      144                  108
17           2840                      528                  400
```

All non-open and bad-reduction residues pass conservatively.

## Bounded results

Before the order-`4` refinement, the following searches were run with the
weaker necessary condition `60 | #J(F_p)`:

- Integer height `5`, primes through `31`: no survivor, agreeing with the
  independent exact Magma screen.
- Integer height `100`, primes through `31`: `25,812` conservative survivors.
  Only `39` were constrained at six or seven genuinely good primes; every
  one was killed at the first available fresh good prime among
  `37,41,43,47`.
- Adding `p=37` left `8,282` conservative height-`100` survivors.  Only six
  had at least seven good constraints; adding `p=41` killed all six.
- Through `p=41`, `2,859` conservative height-`100` survivors remain, but
  every one is non-open/bad at at least three of the ten primes.  There is
  no survivor with eight good constraints.
- Rational height `10` has `127` parameter values and `1,969,008` exact-open
  triples after the cheap discriminant test.  Through `p=41`, `842` pass
  conservatively, but none is constrained at eight good primes.

These are negative search results, not nonexistence proofs: a rational point
may have bad reduction at several selected small primes.

## Fixed root `c=2` matching slice

For an affine match to the compact `M(12)` model, put `x=dZ-c` and fix
`c=2`.  This was the strongest fixed-coordinate slice in the local census.
Its smooth open point counts are

```text
p                 7   11   13   17   19   23   29
# c=2 matches     1    2    1    4    6    5    8
```

`code/m12_general5_fullquad_c2_crt.py` combines all local points, rationally
reconstructs `(b,w,d,s,t,m)`, and accepts only an exact identity

\[
 S_{s,t,m}(Z)=\tau F_{b,w}(dZ-2),\qquad \tau\in\mathbf Q^{\times2}.
\]

Results:

```text
primes 7..23: modulus 7,436,429; 240 local combinations;
               reconstruction bound 2,000; exact hits 0.

primes 7..29: modulus 215,656,441; 1,920 local combinations;
               reconstruction bound 15,000; exact hits 0.
```

The second run reconstructed all six coordinates for `220` CRT combinations,
but none satisfied the exact polynomial identity.  Points bad at one of the
chosen primes are not covered by this CRT statement.

## Current conclusion

The full split-quadratic Hermite cover is locally alive, but neither its
low-height parameter box nor the strongest fixed-root slice has produced a
cyclic order-`60` curve.  The implementation is memory-light (well below
`100 MB`; the residue masks themselves use under `1 MB` through `p=43`).
The next useful bounded run is the refined order-`4` sieve, not a larger run
with only the divisibility condition.

## Refined production runs

The order-`4` exponent refinement was then enabled.  The allowed mask counts
for `p=7,11,13,17,19,23,29,31,37,41` were

```text
6, 50, 108, 400, 500, 1496, 2714, 3178, 5406, 8882.
```

The integer box `|s|,|t|,|m| <= 1000` examined `4,000,000` nonzero `(s,t)`
pairs and `2,157,158,316` CRT candidate lifts.  It left `2,450,907`
conservative survivors, but only `72` were constrained at eight genuinely
good primes and none at nine.  Fresh exact masks killed all `72`: `58` at
`p=43`, `12` at `p=47`, and `2` at `p=53`.  Runtime was about 86 seconds and
peak memory about 32 MB.

The rational height-`30` box contained `1111` parameter values and checked
`1,365,170,136` exact-open triples.  Nine survivors were constrained at eight
good primes and none at nine; fresh masks killed six at `p=43` and three at
`p=47`.  Thus this refined production run also produced no cyclic-`60` hit.
