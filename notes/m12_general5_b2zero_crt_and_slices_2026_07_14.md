# Generic linear-`B` `M(12)+5`: CRT reconstruction and exact slices

## Scope and resource limit

This continues the generic squarefree/coprime `b2=0` sign quotient from
[`m12_general5_b2zero_rootquotient.md`](m12_general5_b2zero_rootquotient.md).
It does not treat the repeated-support or shared-Weierstrass collision
strata. Every Magma run used `MemGB=2` and an external virtual-memory cap of
2.8 GB; the Python reconstruction used a 2.2 GB cap. No run exceeded its cap.

Artifacts:

```text
code/m12_general5_b2zero_crt_reconstruct.py
code/m12_general5_b2zero_slice_probe.m
data/m12_general5_b2zero_crt_reconstruct_2026_07_14.txt
data/m12_general5_b2zero_slice_c1_2026_07_14.txt
data/m12_general5_b2zero_slice_c0_2026_07_14.txt
data/m12_general5_b2zero_slice_d1_2026_07_14.txt
```

## Ten-prime CRT reconstruction

The Python script independently enumerates all fully open quotient points
for which the recovered `tau` is a nonzero square. The counts for

```text
p = 7,11,13,17,19,23,29,31,37,41
```

are

```text
1,1,1,2,5,5,9,16,9,11.
```

Thus there are only `712800` ways to select one local quotient orbit at all
ten primes. Coordinate-wise CRT followed by canonical balanced rational
reconstruction tested every selection. A reconstructed tuple was accepted
only after exact substitution in all four characteristic-zero quotient
equations and an exact rational-square check for `tau`.

The ten-prime modulus is

```text
10141675450907
```

and its canonical reconstruction bound is `2251852` for the numerator and
denominator of each of `(b,w,c,d,e,tau)`. There were `1242` simultaneous
coordinate reconstructions and **no exact point**.

Repeating with every nine-prime subset permits one of the ten primes to be a
bad/boundary prime. In total this tested `3760470` CRT selections and
produced `10060` simultaneous coordinate reconstructions, again with **no
exact point**. The smallest nine-prime balanced bound is `351680`; the
largest is `851120`.

This is a precise bounded statement about the quotient coordinates, not a
height bound solely in the base pair `(b,w)`. It also does not cover a point
having bad/boundary reduction at two or more of the ten primes.

## Exact `c=1` slice

The condition `c=1` puts the root `x=-c` of `B` at the visible marked point
`x=-1` of the compact `M(12)` presentation, making it the most structurally
motivated generic hyperplane.

After specializing before elimination and saturating only by the cheap chart
factors, the slice is zero-dimensional. Its lexicographic eliminant in `e`
has degree `242` and factors over `Q` as two irreducible factors of degrees

```text
46 + 196.
```

There is no linear factor, hence no rational point on this cheap-open slice.
This factorization is a fiber statement; by itself it does not prove that the
global quotient curve has two components.

## `c=0` and `d=1` slices

Both cheap-open saturated slices were certified zero-dimensional:

```text
c=0: final grevlex basis length 157
d=1: final grevlex basis length 182
```

Their rational lex conversions did not finish within the 300-second wall cap
and were stopped. In particular, no rational-point conclusion is drawn from
those two slices. The partial transcripts record the completed exact
saturations.

The finite-field census gives an additional warning about `d=1`: it has no
fully open signed reduction at `p=11,13,17,23,31,41`. Therefore a rational
point on that exact slice would have to enter a denominator or excluded
boundary chart at each of those primes; this is evidence against the slice,
not a global obstruction.

## Decision

The generic squarefree/coprime linear-`B` curve is still locally viable, but
it now looks arithmetically cold:

1. no balanced ten-prime reconstruction through coordinate height
   `2251852`;
2. no reconstruction through the nine-prime tests allowing one bad prime;
3. no rational point on the marked-root slice `c=1`;
4. high (`242`) degree already on a simple hyperplane fiber.

This does not prove that the generic curve has no rational point. It does
make further blind height search or monolithic characteristic-zero
decomposition a poor use of memory. For actually obtaining cyclic `[60]`,
the previously uncomputed shared-Weierstrass and repeated-support collision
strata should rank ahead of this generic open. If the generic lane is
revisited, the next credible step is modular normalization/component genus,
not a larger CRT search.
