# M(2,2,2,6): trying to halve a 2-torsion point

This is the first route: impose that one of the rational 2-torsion classes is divisible by 2.

Use the same odd quintic model as in `m2226_order6_doubling.md`:

```text
Y^2 = L1 L2 L3 L4 L5,
P = (0,2),
```

with branch points on `P^1`

```text
b0 = infinity = (1:0),
b_i = (-A_i : B_i),
(A_1,A_2,A_3,A_4,A_5) = (1,1,1,2,2),
B_1 = 2s^2 - sn,
B_2 = 2s^2 + sm - 2sn - mn,
B_3 = 2s^2 + sm - sn - mn,
B_4 = -mn,
B_5 = 4s^2 - 4sn - mn.
```

For an ordered pair `(b_inf,b_zero)`, put

```text
lambda_k = det(b_k,b_zero)/det(b_k,b_inf)
```

for the four remaining branch points.  The 2-torsion class `b_zero - b_inf` is divisible by 2 if and only if, after choosing one remaining point as reference, the other three ratios `lambda_k/lambda_ref` are squares.  This is just the standard full-rational-2-torsion criterion after moving `b_inf` to infinity and `b_zero` to zero.

## Search

`code/m2226_order2_halving_search.m` tests all 15 unordered 2-torsion classes.

A first pass found many apparent points, but they were all boundary points where one transformed finite branch point moves to infinity.  Concretely the visible families were only the classes `[1,4]` and `[1,6]`, and they forced `B_3=0` or `B_5=0`, respectively, so the original sextic is singular.

After adding the required nonsingularity filter `B_i != 0` for all `i`, plus branch distinctness, the search found no nondegenerate hits up to height 50.

## Local good-reduction sieve

`code/m2226_order2_halving_local_sieve.m` checks the same squareclass conditions over finite fields on the good-reduction locus.

For every one of the 15 classes, there is no good `F_7` point satisfying the halving conditions.  There are 10 good parameter points in `P^2(F_7)`, and all fail the squareclass test for every class.

This is a good-reduction obstruction modulo 7, not yet a complete proof over `Q`: a rational point could reduce into the boundary at 7, so a fully rigorous exclusion would need a small p-adic/blowup analysis of the boundary residue classes.  Empirically, the direct height-50 search and the boundary behavior both point to the same conclusion: this first route does not produce a nondegenerate rational curve in the family.
