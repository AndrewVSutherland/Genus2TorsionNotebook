# W1: component-complete local obstruction

Date: 2026-07-19

## Result

The entire repeated fiber

\[
  W1=(-2254,-2162,2303,4900T^2)
\]

is excluded from the search for rational torsion `[2,2,2,24]`.  This is an
intrinsic local obstruction at (p=13), not a finite Mordell--Weil-bank or
sampled tangent-bank result.  There are **no surviving coefficient classes**
in the rank-three quotient, at any height.

## Complete finite chart

For a signed tuple ((a,b,c,d)), the full `A(2,2,2,8)` cover requires the
four quantities

\[
\begin{aligned}
R_0&=abcd,\\
R_1&=a(a+b)(a+c)(a+d),\\
R_2&=b(b+a)(b+c)(b+d),\\
R_3&=c(c+a)(c+b)(c+d)
\end{aligned}
\]

to be squares.  Exhausting all (T\bmod 13) leaves exactly

\[
  T\equiv 1,4,9,12=\pm1,\pm4\pmod {13}.
\]

All four radicands are nonzero on these disks, and the four squared branch
parameters are distinct, so the genus-two curve has good reduction.  There
are only two reduced curves, according to (T^2\), and both have

\[
  \#J(\mathbf F_{13})=128\equiv2\pmod3.
\]

Because the kernel of good reduction is pro-(13), it has no (3)-torsion.
Consequently (J(\mathbf Q_{13})[3]=0) on every finite full-cover disk.
Direct enumeration of the corrected cubic-contact equations independently
gives zero contact presentations on all four disks.

The residue details are:

| (T\bmod13) | cover radicands | (#J(\mathbf F_{13})) | contact points |
|---:|---|---:|---:|
| 1, 12 | `(12,4,9,12)` | 128 | 0 |
| 4, 9 | `(10,1,10,1)` | 128 | 0 |

## Complete infinity chart

If (v_{13}(T)<0), put (z=1/T) and projectively rescale the signed tuple to

\[
  (-2254z^2,-2162z^2,2303z^2,4900).
\]

After removing the common even power (z^6), the four limiting cover units
modulo (13) are

\[
  (12,5,7,1).
\]

The middle two units are nonsquares.  Thus the full cover has no point in
the nonintegral-(T) chart.  Together with the finite calculation, this
exhausts (mathbf P^1(\mathbf Q_{13})).

## Mordell--Weil congruence mask and periodic bank

The pair `<1,3>` quotient is

\[
 E:y^2=x^3-x^2+34591x+34593,
\]

of rank three with torsion `[2]`.  Although this quotient has awkward
reduction at (13), no quotient-map calculation is needed: every possible
value of (T) is already excluded on the target genus-two fiber.  Hence the
complete congruence condition on its coefficients is simply

\[
  \boxed{\text{FALSE for every }(m,n,k,ti)}.
\]

The existing periodic (N=1000) computation represents exactly

\[
  2(2\cdot1000+1)^3=16,024,012,002
\]

raw coefficient/torsion classes.  Earlier profiles left 1,697 periodic rows,
and the extension through (p=167) left the single row
`(2,-1,1,ti=2)`.  The complete (p=13) obstruction sends both banks to the
empty set.  The exact reconstruction

\[
T=-6015613/4711271
\]

and its nonsquare middle cover ratio are therefore consistent but redundant.

This is stronger than any enlargement of the periodic box: W1 has no
rational specialization yielding the target, at any Mordell--Weil height.

## Files

```text
code/target_22224_W1_p13_complete_obstruction.m
results/target_22224_W1_p13_complete_obstruction.log
results/target_22224_W1_p13_complete_obstruction.tsv
results/target_22224_W1_p13_coefficient_mask.tsv
results/target_22224_W1_rank3_periodic_N1000_p13complete.tsv
```

