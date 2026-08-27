# Frobenius and endomorphisms of the genus-6 trigonal quotient

## Outcome

The Jacobian `J` of the genus-6 quotient is **absolutely simple**, and in
fact

```text
End_{Qbar}(J) = Z.
```

Both statements are rigorous certificates, not heuristics.  In particular,
the trigonal curve has no positive-genus lower-dimensional quotient over
`Qbar`, and its Jacobian has no elliptic, genus-2, genus-3, or other proper
isogeny factor over any number field.  There is therefore no hidden quotient
or correspondence that reduces the rational-point problem to smaller
Jacobians.

The two good reductions also give

```text
gcd(#J(F_13), #J(F_17)) = gcd(5716043,25773047) = 1,
```

so `J(Q)_tors` is trivial.

## Exact Frobenius polynomials

Write

```text
P_p(X) = det(X-Frob_p | H^1_et).
```

At `p=13`, normalization after reduction still has genus 6, and

```text
P_13 = X^12 + 2 X^11 + 2 X^10 + 15 X^9 + 302 X^8
       + 341 X^7 + 437 X^6 + 4433 X^5 + 51038 X^4
       + 32955 X^3 + 57122 X^2 + 742586 X + 4826809.
```

This polynomial is irreducible over `Q`, is ordinary, and has exact Galois
group

```text
12T293,  order 46080 = 2^6 * 6!,
```

the full hyperoctahedral group.  Here

```text
#D(F_13) = 16,       #J(F_13) = 5716043.
```

At `p=17`, the corresponding data are

```text
P_17 = X^12 + 2 X^11 - 14 X^10 - 18 X^9 + 200 X^8
       + 31 X^7 - 5037 X^6 + 527 X^5 + 57800 X^4
       - 88434 X^3 - 1169294 X^2 + 2839714 X + 24137569,

#D(F_17) = 20,       #J(F_17) = 25773047.
```

Again the polynomial is irreducible and ordinary, with Galois group
`12T293` of order `46080`.

The branch-divisor discriminant factors as

```text
2^302 * 3^172 * 5^15 * 7 * 11 * 269^3,
```

so neither 13 nor 17 is a branch-collision prime.  Direct normalization in
both characteristics gives the unchanged genus 6.

## Why absolute simplicity is proved

Let the roots of `P_p` be arranged as

```text
alpha_1, p/alpha_1, ..., alpha_6, p/alpha_6.
```

The full hyperoctahedral Galois group contains every independent swap
`alpha_i <-> p/alpha_i`.  If two distinct Frobenius roots acquired equal
`n`-th powers over some finite extension, applying one independent swap
would force

```text
(alpha_i^2/p)^n = 1.
```

That would give `alpha_i` slope `1/2`, contradicting ordinarity.  Hence the
powered roots remain distinct for every `n`; transitivity of the full group
then makes every extension Frobenius polynomial irreducible.  Thus the
reduction remains simple over every finite extension and is absolutely
simple.  Good-reduction specialization then proves that `J/Qbar` itself is
absolutely simple.

As a finite sanity check, the script also verifies irreducibility after
extensions of every degree `2 <= n <= 12` at both primes.

## Why the geometric endomorphism ring is Z

For each ordinary absolutely simple reduction, its geometric rational
endomorphism algebra is the CM Frobenius field `K_p=Q(alpha_p)`.  Exact
subfield computations show that `K_13` and `K_17` each have only two
nontrivial subfields listed by Magma: the degree-6 maximal real field and
the full degree-12 CM field.

The real Weil polynomials are

```text
p=13: y^6 + 2 y^5 - 76 y^4 - 115 y^3
      + 1719 y^2 + 1446 y - 11133,

p=17: y^6 + 2 y^5 - 116 y^4 - 188 y^3
      + 3753 y^2 + 3839 y - 29755.
```

Their maximal-order discriminants are respectively

```text
117696160151333626425
4659781303841945487531457,
```

so the real fields are not isomorphic.  The full CM-field discriminants
are also unequal.  Specialization embeds `End^0_{Qbar}(J)` into both CM
fields.  Since neither possible nontrivial subfield is common to them, this
endomorphism algebra is `Q`, and integrally `End_{Qbar}(J)=Z`.

## Canonical-model probe

The independent canonical computation finishes quickly.  The normalized
curve is nonhyperelliptic and its canonical model has

```text
ambient space:             P^5
canonical degree:          10
canonical ideal degrees:   six quadrics and three cubics.
```

The six quadrics cut out a nonsingular rational normal scroll of dimension
2 and degree 4.  `CliffordIndexOne` recovers a map to `P^1`, namely the
expected unique trigonal pencil.  No extra canonical component or second
low-degree pencil appears.  The scroll records trigonal geometry but does
not split the Jacobian.

## Conductor and analytic rank status

The displayed branch discriminant identifies the branch-collision primes
`{2,3,5,7,11,269}`, but it is not by itself a conductor calculation.
Determining conductor exponents and local root numbers requires stable
reduction at those primes.  No rigorous conductor, global root number, or
analytic-rank claim is made here.

The remaining plausible global route is therefore direct arithmetic on the
full trigonal Jacobian: establish `rank J(Q)` and use Chabauty if it is less
than 6, or use a trigonal descent.  Lower-dimensional quotient methods are
now ruled out.

## Reproduction

```bash
magma -b Prime:=13 ComputeGalois:=true \
  code/contact30_c3root_genus6_frobenius.m
magma -b Prime:=17 code/contact30_c3root_genus6_frobenius.m
magma -b code/contact30_c3root_genus6_endomorphism_certificate.m
magma -b code/contact30_c3root_genus6_canonical.m
```

The first two commands regenerate the expensive zeta functions.  The
two-prime endomorphism certificate runs from the recorded exact polynomials
in a few seconds.
