# S1: component-complete obstruction

Date: 2026-07-19

## Result

The entire repeated fiber

\[
  S1=(-4,9,30,(-2166/245)T^2)
\]

is excluded from the search for rational torsion `[2,2,2,24]`.  The complete
obstruction occurs at (p=71), where every finite full-cover residue has
good reduction and Jacobian order prime to (3), while the nonintegral-(T)
chart has no full-cover point.

Thus there are **no S1 Mordell--Weil coefficient classes at any height**.
This supersedes the earlier conclusion relative to the sampled 144-key
(13^5) tangent bank.

## Integral model and full-cover equations

Multiplying all four signed branch parameters by (245) gives the
projectively equivalent integral tuple

\[
  (a,b,c,d)=(-980,2205,7350,-2166T^2).
\]

This common scaling does not alter the full-cover square classes or the
genus-two Jacobian up to rational isomorphism.  It is a unit at (71).

The full `A(2,2,2,8)` cover requires

\[
\begin{aligned}
R_0&=abcd,\\
R_1&=a(a+b)(a+c)(a+d),\\
R_2&=b(b+a)(b+c)(b+d),\\
R_3&=c(c+a)(c+b)(c+d)
\end{aligned}
\]

to be squares.  Exhaustion of all (T\bmod71) leaves exactly twelve
residues:

\[
  T\equiv\pm1,\ \pm4,\ \pm9,\ \pm20,\ \pm33,\ \pm34\pmod {71}.
\]

Every cover radicand is nonzero on these disks, and the squared branch
parameters are pairwise distinct.  Hence every disk has good reduction.

## Jacobian orders and contact components

The reduced Jacobian orders, grouped by the sign symmetry (T\mapsto-T),
are:

| (T\bmod71) | (#J(\mathbf F_{71})) | modulo (3) |
|---:|---:|---:|
| (pm1) | 5632 | 1 |
| (pm4) | 4096 | 1 |
| (pm9) | 5312 | 2 |
| (pm20) | 4480 | 1 |
| (pm33) | 4096 | 1 |
| (pm34) | 3584 | 2 |

For good reduction at (71\ne3), the reduction kernel is pro-(71).
Therefore it contains no (3)-torsion, and the displayed orders prove

\[
  J(\mathbf Q_{71})[3]=0
\]

on every finite full-cover disk.  In particular, there is no cubic-contact
component over any of them.  As an independent formula audit, the code also
enumerates the corrected (M=L^2\), (L\ne0) cubic-contact equations over
(\mathbf F_{71}); their count is zero on all twelve residues.

The two known rational full-cover parameters occur among these disks:

\[
  T=1\longmapsto1,\qquad T=7/19\longmapsto34\pmod {71},
\]

and are killed by orders (5632) and (3584), respectively.

## Infinity chart

For (v_{71}(T)<0), put (z=1/T) and rescale the signed tuple to

\[
  (-980z^2,2205z^2,7350z^2,-2166).
\]

After removing the common even power (z^6), the limiting full-cover units
modulo (71) are

\[
  (29,35,15,47).
\]

Their square flags are `(true,false,true,false)`.  Thus the full cover has
no point in the nonintegral-(T) chart.  This and the finite calculation
exhaust (\mathbf P^1(\mathbf Q_{71})).

## Mordell--Weil mask and (10^8) bank

The pair `<1,2>` quotient is the rank-one curve

\[
 E:y^2=x^3+x^2+3392x+62288,
 \qquad E(\mathbf Q)=\langle(152,-2028)\rangle\oplus\langle(-17,0)\rangle.
\]

The obstruction occurs on the target genus-two curve, below the quotient
map.  Consequently the complete coefficient condition is simply

\[
  \boxed{\text{FALSE for every }(m,ti)}.
\]

The existing (|m|\le10^8) run represents

\[
  2(2\cdot10^8+1)=400,000,002
\]

raw coefficient/torsion classes.  Its finite-prime periodic sieve enumerated
2,051,283 classes and retained 3,412 rows.  The intrinsic (p=71) mask sends
all 3,412 to the empty set.  There are no classes requiring auxiliary-prime
or exact reconstruction, and increasing the Mordell--Weil bound cannot
reopen S1.

## Files

```text
code/target_22224_S1_auxprime_scan.m
code/target_22224_S1_p71_complete_obstruction.m
results/target_22224_S1_auxprime_scan.log
results/target_22224_S1_auxprime_scan.tsv
results/target_22224_S1_p71_complete_obstruction.log
results/target_22224_S1_p71_complete_obstruction.tsv
results/target_22224_S1_p71_coefficient_mask.tsv
results/target_22224_S1_rank1_N100000000_p71complete.tsv
```

