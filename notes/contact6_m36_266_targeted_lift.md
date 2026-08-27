# Targeted non-automorphic `[2,6,6]` lifting diagnostic

Date: 2026-07-10

This is a bounded follow-up to the extra-root/cubic-contact search.  It does
not extend the rational height box.  It tests two fixed `(eps,r,b)` fibers of
the three contact equations in `(M,U,v)`, where `M=L^2`:

```text
F1(M,U,v) = F2(M,U,v) = F3(M,U,v) = 0.
```

The reproducible implementations and output are:

```text
code/contact6_m36_266_targeted_lift_sage.py
code/contact6_m36_266_targeted_lift.m
data/contact6_m36_266_targeted_lift.txt
```

The completed run was:

```text
sage code/contact6_m36_266_targeted_lift_sage.py \
  --skip-exact --lift-precision 6
```

## Correction to the `p=19`, `r=2`, `b=3` interpretation

The point with `(eps,r,b)=(1,2,3)` modulo `19` has an irreducible Frobenius
quartic, but it is **not** a non-automorphic finite point.  At `r=2,b=3`, both
explicit Mobius-pairing equation sets vanish modulo `19`:

```text
AutoA = [4256,-6992,608]  = 0 mod 19,
AutoB = [12768,-6992,608] = 0 mod 19.
```

Thus the irreducible quartic is compatible with elliptic factors that become
conjugate over the finite field.  It cannot be used as evidence that this
residue avoids the Humbert/extra-involution locus.  Over `Q` the same fixed
base `(r,b)=(2,3)` is not on either explicit automorphism locus, since the
displayed integers are nonzero.

For a genuinely non-automorphic `p=19` test, the run used the existing finite
sample

```text
(eps,r,a,b,L,U,v) = (1,4,18,5,5,9,3) mod 19.
```

The exact base is `(r,a,b)=(4,-59/2,5)`.  Neither automorphism equation set
vanishes over `Q` or modulo `19`.

## Smooth roots and Hensel lifts

For `(r,b)=(2,3)`, all good smooth roots in the selected fibers are:

```text
p=13: (M,U,v) = (4,0,5), (9,0,8)
p=19: (M,U,v) = (7,3,14), (7,16,5)
```

For the non-automorphic base `(r,b)=(4,5)` they are:

```text
p=13: (M,U,v) = (1,4,6), (1,9,7)
p=19: (M,U,v) = (6,9,3), (17,9,16)
```

Every listed root has invertible `3 x 3` Jacobian and lifted uniquely through
`p^6`.  The wider prime scan for `(4,5)` found no smooth good roots at
`p=11,17,29,31`; it found one at `p=23`, but that reduction is on the explicit
automorphism locus.

## CRT result and scope

All four pairings of the `p=13` and `p=19` lifts were combined.  The modulus
and balanced rational-reconstruction bound were

```text
13^6 * 19^6 = 227081481823729,
|numerator|, denominator <= 10655549.
```

Results:

```text
(r,b)=(2,3): 4 CRT pairs, 1 coordinate-wise reconstruction, 0 exact points
(r,b)=(4,5): 4 CRT pairs, 0 coordinate-wise reconstructions, 0 exact points
```

Consequently, neither fixed fiber contains a rational `(M,U,v)` point within
that balanced bound whose denominators are prime to `13*19` and whose
reductions are among these smooth good roots.  This is a bounded negative
diagnostic, not a proof that either fiber has no rational point.

## Exact-component and verification limits

An exact saturated Groebner calculation over `Q` for the first fixed fiber did
not finish within the 90-second cap and was interrupted.  Primary decomposition
was not attempted after that stop condition.  Magma was not installed on this system.  Since there was no
exact CRT hit, no exact `TorsionSubgroup` or Sage/Lombardo simplicity check was
triggered.

The next justified step is not a deeper lift of these same isolated fibers.
It is to saturate the full `(r,b,M,U,v)` contact surface by the two explicit
automorphism ideals and boundary factors, decompose the surviving locus modulo
`13` and `19`, and lift a positive-dimensional component with compatible
tangent directions.  Fixed-fiber CRT should resume only after such a component
supplies a reason to pair residues across primes.
