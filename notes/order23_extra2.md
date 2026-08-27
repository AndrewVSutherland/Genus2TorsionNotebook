# Kuru-Sadek [23] family: rational 2-torsion condition

This note records the first attempt to add an independent rational `2`-torsion
point to the Kuru-Sadek genus-2 family with rational `23`-torsion.  A success
would give a rational point of order `46` on a simple Jacobian.

## Family

For `t != 0, +/-1`, set

```text
beta  = (t^2 + 1)^2/(4*t^2)
sbeta = (t^2 + 1)/(2*t)
s      = (t^2 - 1)/(2*t)
alpha = beta - s^5/(beta*sbeta)
lambda = (alpha - 1)^4/((alpha - beta)^2*alpha)
```

and

```text
A(x) = (x^3*(x-alpha)^2
        - (x-1)*((x-1)^4 - lambda*(x-beta)^2*x))
       /(2*(x-alpha)*(x-beta)),
C_t: y^2 = f_t(x) = A(x)^2 - lambda*x^4*(x-1).
```

The odd quintic model already has the rational Weierstrass point at infinity.
An independent rational `2`-torsion point is therefore equivalent to a finite
rational branch point, i.e. a rational root `r` of `f_t(x)`.

## Algebraic condition

After clearing denominators in `f_t(r)` and using the symmetry `t -> -t`, the
condition descends to `q=t^2`.  The result is a single irreducible equation

```text
G(q,r) = 0,
```

with

```text
deg_q(G) = 18,
deg_r(G) = 5,
total_degree(G) = 23,
number_of_terms = 104.
```

The affine curve has geometric genus `16`.  Its affine singular points are

```text
(q,r) = (-1,0), (0,0), (0,1), (1,0), (1,1).
```

These are boundary points for the original family or not rational-square `q`
values.  Thus rational `2`-torsion on this `[23]` family is not a low-degree
auxiliary condition; it is a high-genus plane curve plus the extra condition
that `q` is a rational square.

The denominator in the original `t` parameter factors as

```text
(t - 1)^4 * t^4 * (t + 1)^4 * (t^2 + 1)^6
* (t^8 + 2*t^4 + 1/5)^2.
```

Over `Q`, the only rational boundary parameters are `t=0,+/-1`; the other
factors matter in finite-prime reductions.


Follow-up quotient and boundary analysis is in

```text
notes/order23_extra2_quotients.md
```

Main derivation script:

```text
sage code/order23_extra2_conditions.sage
```

A full text output, including the polynomial factorization and finite-prime
filters, is in

```text
data/order23_extra2_conditions.txt
```

## Finite-prime filters

The good affine chart has no square-`q` points at `p=3` or `p=5`; any rational
point must reduce to a boundary/bad residue there.  For larger primes the
condition gives small allowed sets of `q=t^2` residues.  For example:

```text
p=7:  q in {4}
p=11: q in {3,5}
p=13: q in {4,9,10}
p=17: q in {2,8,13,15}
p=23: q in {4,8,12}
```

The search script uses these residue sets through `p=97`, allowing boundary
residues separately.

## Search

The search is not a blind `(t,r)` box.  It enumerates rational `t` of bounded
height, applies the finite-prime residue filter coming from `G(q,r)=0`, and
only then factors the exact quintic `f_t(x)` over `Q`.

Run:

```text
sage code/order23_extra2_search.sage 2000
```

Output summary:

```text
height 2000
checked 4866348
survivors 3704
factored 3704
singular 0
hits 0
```

So there is no rational finite branch point for `t` of height up to `2000`
after the algebraic filter.

## Interpretation

The `[23]` family itself is very good: it has exact `[23]` examples and simple
Jacobians.  But adding rational `2`-torsion looks arithmetically restrictive.
The obstruction is not merely search size: the rational-root condition is an
irreducible degree `(18,5)` curve of geometric genus `16`, and small primes
force many reductions onto bad or narrow residue classes.

The next serious step, before any larger search, would be to study the genus-16
curve `G(q,r)=0`: look for maps or low-genus quotients, and try to prove
whether it has any rational points with `q` a square away from the boundary.
