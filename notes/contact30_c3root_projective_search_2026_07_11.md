# Contact-30 cubic-root cover: projective search

## Result

The rational-root gate for the distinguished Richelot route to `[2,60]` has
no smooth point on either order-30 branch through parameter height `25000`.
The exhaustive box was

```text
R=a/b,  b>0,  gcd(a,b)=1,  |a|<=25000,  b<=25000.
```

It contains `759,920,359` rational parameters, or `1,519,840,718` signed
branch tests.  A 61-prime projective sieve left four parameters on each sign:

```text
R = 1, 2, 3, 7/3.
```

All four are exact global boundary points.  After removing them there are
zero open survivors and hence no rational-root candidate requiring a
Jacobian computation.

## What was sieved

For each parameter and sign, write

```text
C3(x) = 2*x^3 + (A-3)*x^2 + (B+3)*x + (C-1).
```

At every selected odd prime, every point of `P^1(F_p)`, including infinity,
was classified as:

1. a good smooth fiber where `C3` has an `F_p` root;
2. a good smooth fiber where it has no root; or
3. a boundary fiber (`R`-pole, `u=0`, `c=0`, `q`-pole, wrong degree, or
   singular genus-2 model).

Classes (1) and (3) were retained; only class (2) was rejected.  This matters:
a rational parameter can reduce into a boundary disk at some prime, and such
a disk cannot safely be treated as a local obstruction.  The program prints
the boundary disks separately at every prime.  Consequently the final empty
open intersection is an exhaustive bounded exclusion, not a good-reduction-
only heuristic.

At a good residue, the implication used by the sieve is exact.  Since the
leading coefficient of `C3` is `2`, a rational root is integral at every odd
good prime and reduces to an affine `F_p` root.

As an independent implementation check, the existing Magma Richelot code was
rerun at `p=11,13,17,23`.  Its numbers of affine good parameters with a `C3`
root agree with the new C++ masks; after adding the projective infinity fiber,
the counts agree exactly on both signs.

## Exact boundary verification

The four survivors arise from

```text
c = (u^2-1)/(2*u) = 0.
```

Over `Q`, this forces `t=+/-1`, and direct factorization gives

```text
t=+1:  (R-2)*(R-3)=0,
t=-1:  (R-1)*(R-7/3)=0.
```

At `R=1,7/3`, the recovery denominator for `q` vanishes.  At `R=2,3`, the
formal cubic factors as `x^2*(x+1)`, but the source curve is degenerate
because `c=0`.  Thus even the apparent rational roots at `R=2,3` do not give
points of the open trigonal cover.

The exact assertions and factorizations are in
`code/contact30_c3root_boundary_verify.m`.

## Interpretation

A rational root of `C3` is necessary before constructing the distinguished
Richelot neighbor and testing its dual 2-classes.  Therefore there is no
`[2,60]` candidate from this route in the stated height box.  This does not
prove that the trigonal cover has no rational points outside the boundary;
it is a bounded result.  Enlarging the box is cheap compared with Chabauty,
but the complete absence of an open survivor after more than 1.5 billion
signed tests is substantial negative evidence.

## Reproduction

```bash
c++ -O3 -std=c++17 code/contact30_c3root_projective_sieve.cpp \
    -o /tmp/contact30_c3root_projective_sieve
/tmp/contact30_c3root_projective_sieve 25000 quiet
magma -b code/contact30_c3root_boundary_verify.m
```

Files:

- `code/contact30_c3root_projective_sieve.cpp`
- `code/contact30_c3root_boundary_verify.m`
- `data/contact30_c3root_projective_sieve_h25000.txt`
