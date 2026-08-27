# Marked local sieve for cyclic `[49]` in the contact-7 family

This records a bounded-memory search in the two-parameter family

```text
h = 1 - (7/2)*x + a*x^2 + b*x^3,
f = (h^2 + (x-1)^7)/x^2,
D7 = [(x-1),h(1)].
```

The implementation and recorded run are:

```text
code/z49_local_marked_h50_sieve.m
data/z49_local_marked_h50_p43.txt
data/z49_local_marked_h50_p43_time.txt
```

## Correct finite condition

It is not enough that `49` divide `#J(F_p)`.  It is not even enough that
the exponent of `J(F_p)` be divisible by `49`: the resulting order-49
subgroup could have the wrong order-7 subgroup.

If the rational torsion group of a smooth member is cyclic `[49]`, then its
marked nonzero class `D7` generates the unique order-7 subgroup.  Therefore

```text
D7 is in 7*J(Q).
```

At a good prime `p != 7`, this relation reduces to

```text
D7_bar is in 7*J(F_p).
```

The code computes `G,phi := AbelianGroup(J)` and the invariant-factor
coordinates of `D7 @@ phi`.  In a coordinate `Z/nZ`, the equation
`7*q=d` is soluble exactly when `gcd(7,n)` divides `d`.  Thus the mask is an
exact marked-divisibility test, not an order surrogate.  Singular fibers and
parameters nonintegral in the displayed chart are retained as unresolved;
only good marked fibers failing this test are rejected.

The local counts demonstrate the distinction:

```text
p    49|#J    49|exponent    D7 in 7J    exponent-only false positives
29     224        107            99                    8
37     228        201           196                    5
43     463        240           218                   22
```

For `p=3`, five of the nine affine fibers are good and none passes marked
7-divisibility.  Hence any rational cyclic-49 specialization in this chart
must enter a singular or nonintegral 3-adic boundary residue.

## Height-50 search

The parameter set consists of every reduced rational `n/d` with

```text
|n| <= 50,  1 <= d <= 50.
```

There are `3095` values for each of `a,b`, hence `9,579,025` pairs.  The
search used the primes

```text
3,5,11,13,17,19,23,29,31,37,41,43.
```

Masks and rational residue tables are precomputed, but parameter pairs are
streamed and never stored.  The result was

```text
checked               9,579,025
raw mask survivors        1,010
smooth survivors               0
```

Thus there is no smooth contact-7 member of height at most `50` satisfying
these necessary marked local conditions, and in particular no cyclic `[49]`
example in the box.

## Classification of the raw survivors

The discriminant factors as

```text
256*Disc(f) = (2*a+2*b-5)^7 * Q5(a,b),
```

where `Q5` is written explicitly in the code.  All `1,010` raw survivors
are singular:

```text
only 2*a+2*b-5 = 0       1,006
only Q5(a,b) = 0             3
intersection of both             1
unclassified                     0
```

The dominant component is exactly `h(1)=0`, where the marked contact point
collides with the singular fiber.  The finite masks are deliberately
conservative and retain such residues, and the characteristic-zero check
then removes them.

## Resources and scope

The recorded run completed with exit status zero in `27.01` seconds and used
`26,880 KB` maximum resident memory, under a hard `200 MB` address-space cap.

This closes the stated rational-height box rigorously.  It is not a global
nonexistence result: a cyclic `[49]` specialization of larger height would
still have to negotiate the forced 3-adic boundary behavior.
