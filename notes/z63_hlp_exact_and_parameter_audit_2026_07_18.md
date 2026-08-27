# Exact split HLP `[63]` and parameter audit

Audit date: 2026-07-18.

## Bottom line

There was no targeted `[63]` construction or search in the repository before
this audit.  The known Howe--Leprévost--Poonen (HLP) example has now been
transcribed and verified exactly:

```text
C: y^2 = 897*x^6 - 197570*x^4 + 79136353*x^2 - 146398496,
J(C)(Q)_tors = Z/63Z.
```

It is geometrically split.  The even sextic has the non-hyperelliptic
involution `x -> -x`, and its two elliptic quotients have rational torsion
groups `[9]` and `[7]`.  Equivalently, the HLP construction gives a `(2,2)`
isogeny from the product of those elliptic curves to `J(C)`.

No new underlying 7-by-9 gluing pair, hence no new `[63]` curve, was found in
the exact rational-parameter search through height `200` described below.

## Exact marked divisor and full torsion

HLP print a conjugate pair of points with

```text
x_R = (-69+sqrt(4369))/2,
y_R = 4515015-68241*sqrt(4369).
```

The rational Mumford representation of their sum minus the two points at
infinity is therefore

```text
D = [x^2+69*x+98, -136482*x-193614].
```

The verifier checks directly that `v^2-f` is divisible by the support
polynomial and that

```text
63*D = 0,   21*D != 0,   9*D != 0.
```

Thus `D` has exact order `63`.  Independently, Magma computes the complete
rational torsion invariants as `[63]`.  It also reproduces HLP's sharp
good-reduction certificate

```text
#J(F_5) = 63.
```

The run took `0.19` wall seconds and `82,184 KB` peak RSS under a hard
`250,000,000`-byte address-space cap.

## Exact bounded parameter audit

HLP obtain the curve by gluing the universal elliptic curve with a marked
7-torsion point at `t=-16/3` to the universal curve with a marked 9-torsion
point at `u=4`.  A necessary condition for their 2-torsion cubic fields to be
isomorphic is equality of the rational squareclasses

```text
Delta_7(t) = t*(t-1)*(t^3-8*t^2+5*t+1),
Delta_9(u) = u*(u-1)*(u^2-u+1)*(u^3-6*u^2+3*u+1).
```

The exact sieve enumerated every reduced rational `a/b` with
`max(|a|,b) <= 200`, namely `48,927` parameters on each side.  It found

```text
24,462 nonzero squareclasses on the 7-side,
24,463 nonzero squareclasses on the 9-side,
2 common squareclasses,
18 marked parameter pairs.
```

The 18 pairs are two `3 x 3` blocks:

```text
t in {-16/3, 19/16, 3/19}, u in {4, -1/3, 3/4}, squareclass -741;
t in {-2, 3/2, 1/3},       u in {-1, 2, 1/2}, squareclass -6.
```

Exact `j`-invariant comparisons show that the three parameters in each set
are merely three markings of the same underlying elliptic curve.  All nine
pairs in the first block have isomorphic 2-torsion cubic fields, both of
field discriminant `-2964`; this is precisely the known HLP pair.  The cubic
field is non-Galois with `S3` root action, so the centralizer of that action
is trivial and the Galois-equivariant identification of the two root sets is
unique.  There are not three additional gluings hidden at this seed.

The second block is eliminated by the sufficient cubic-field test: the two
fields have discriminants `-1176` and `-216` and are not isomorphic.  Thus the
height-200 audit produces no new underlying gluing pair.

The sieve used `70,692 KB` peak RSS and `5.20` wall seconds.  The cubic-field
test used `26,624 KB` and `0.04` wall seconds, both under the same 250 MB cap.

## Is a family realistic?

Not by this HLP 7-by-9 route.  HLP Table 1 labels `[63]` by `P^0`, meaning a
single curve rather than a positive-dimensional family, and Section 3.6 says
that only finitely many rational `[63]` curves seem likely and that their
curve might even be unique.  This is not a theorem, but the height-200 exact
audit supports their assessment: after the cheap necessary squareclass test,
only the known modular orbit survives the actual cubic-field test.

A larger bounded search is computationally realistic with this sieve, but it
is an isolated-example search with rapidly growing rational height, not a
plausible path to an infinite family.  A geometrically simple `[63]` would
also require leaving the HLP construction entirely, because every HLP output
here is `(2,2)`-split.

## Reproduction

```bash
prlimit --as=250000000 magma -b MemMB:=220 code/z63_hlp_exact_verify.m
prlimit --as=250000000 python3 code/z63_hlp_squareclass_sieve.py --height 200 --show 0
prlimit --as=250000000 magma -b MemMB:=180 code/z63_hlp_cubic_field_test.m
```

Machine-readable concise outputs are in
`data/z63_hlp_exact_verify_2026_07_18.txt` and
`data/z63_hlp_parameter_sieve_h200_2026_07_18.txt`.

Primary source: Everett W. Howe, Franck Leprévost, Bjorn Poonen,
*Large torsion subgroups of split Jacobians of curves of genus two or three*,
Forum Math. 12 (2000), Section 3.6 and equation (4):
<https://math.mit.edu/~poonen/papers/large.pdf>.
