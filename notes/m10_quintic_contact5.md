# Quintic-contact 5-torsion family

This records the first usable odd-prime torsion base for the larger-torsion
search.  The starting point is the elementary quintic-contact construction

```
h(x) = 1 + a*x + b*x^2,
f(x) = h(x)^2 - (1 + a + b)^2*x^5.
```

On `C: y^2 = f(x)`, the point `P=(0,1)` satisfies

```
div(y - h(x)) = 5(P) - 5(infinity),
```

provided `f` is squarefree of degree `5`.  Thus `[P - infinity]` is a
rational `5`-torsion point.  The normalization forces `x=1` to be a root of
`f`; this loses no generality after choosing one rational root away from `P`.

If `f` splits completely over `Q`, then the Jacobian also has full rational
`2`-torsion.  This gives at least

```
(Z/2Z)^3 x Z/10Z.
```

The script implementing the test is

```
code/m10_quintic_contact5_search.m
```

Typical runs:

```
magma -b mode:="finite" code/m10_quintic_contact5_search.m
magma -b mode:="search" height:=30 code/m10_quintic_contact5_search.m
```

## Finite-field split-density check

The finite check asks for smooth members of the family for which the quintic
splits completely modulo `p`.

```
p 3  total 9    smooth 3    full_split 0   five_ok 0
p 5  total 25   smooth 15   full_split 0   five_ok 0
p 7  total 49   smooth 32   full_split 0   five_ok 0
p 11 total 121  smooth 96   full_split 1   five_ok 1
p 13 total 169  smooth 134  full_split 0   five_ok 0
p 17 total 289  smooth 242  full_split 0   five_ok 0
p 19 total 361  smooth 312  full_split 5   five_ok 5
p 23 total 529  smooth 464  full_split 10  five_ok 10
p 29 total 841  smooth 762  full_split 20  five_ok 20
p 31 total 961  smooth 876  full_split 21  five_ok 21
p 37 total 1369 smooth 1262 full_split 30  five_ok 30
p 41 total 1681 smooth 1566 full_split 46  five_ok 46
p 43 total 1849 smooth 1724 full_split 40  five_ok 40
```

The important local consequence is that this normalized full-split family has
no good full-split points modulo `3`, `5`, `7`, `13`, or `17`.  Therefore any
rational full-split specialization in this chart must reduce to the bad or
boundary locus at those primes.  This is a strong congruence constraint, not
just a height-search artifact.

## Rational height search

The first exact rational search used the good-prime split filter at
`3,5,7,11,13,17,19`, then factored the surviving quintics over `Q`.

Run:

```
magma -b mode:="search" height:=25 max_hits:=5 code/m10_quintic_contact5_search.m
```

Output:

```
Quintic-contact 5-torsion rational split search
height 25 parameters 799
filter_primes [ 3, 5, 7, 11, 13, 17, 19 ]
DONE height 25
checked 637802 smooth 637564 filter_survivors 387 full_split 0 torsion_tests 0 hits 0
```

So the family is mathematically sound as a `5`-torsion base, but the naive
two-parameter rational split search found no full rational `2`-torsion
specialization through height `25`.

## Conclusion

This route is more reliable than the provisional `M_1(8)` plus `5` scripts,
because the divisor identity gives a genuine rational `5`-torsion point before
any search begins.  However, asking for full rational `2`-torsion in this
normalized chart is locally forced onto boundary/bad residues at several small
primes.

The next logical computation is not another blind height extension.  It is to
parameterize or explicitly enumerate the boundary charts where the reductions
are bad at the obstructing primes, then test those charts exactly.  In
practice this means using the congruence constraints at `3`, `5`, `7`, `13`,
and `17` as an input to the search rather than as a final filter.


## Boundary-focused search

The `boundary` mode converts the small-prime obstruction into a direct search.
For the primes where there are no good full-split residues, namely
`3,5,7,13,17`, it only keeps rational parameters whose reduction is bad or
outside the affine chart.  At the additional primes
`11,19,23,29,31,37,41,43`, it keeps parameters that are either bad/boundary or
full-split modulo that prime.  Only then does it factor `f` over `Q`.

Run:

```
magma -b mode:="boundary" height:=80 max_hits:=5 code/m10_quintic_contact5_search.m
```

Output summary:

```
height 80 parameters 7863
checked 61820872
boundary_survivors 958333
split_survivors 2327
smooth 0
full_split 0
torsion_tests 0
hits 0
```

The same phenomenon already appeared at lower height:

```
height 10: split_survivors 36,   smooth 0
height 40: split_survivors 581,  smooth 0
height 80: split_survivors 2327, smooth 0
```

Thus, in this normalized quintic-contact chart, the simultaneous boundary
conditions forced by `3,5,7,13,17` appear to collapse the surviving rational
search space onto the global discriminant locus, at least through height `80`.
This is stronger evidence against finding the desired full-split example here
than a blind height search would provide.
