# Contact-9 family

This records the direct contact-9 scaffold tried after the contact-7
`[28] -> [56]` route was largely closed off.

## Construction

On an odd genus-2 model `C: y^2 = f(x)`, the function

```text
x^2*y - h(x)
```

has pole order `9` at infinity.  Force

```text
h(x)^2 - x^4*f(x) = -(x - 1)^9.
```

Taking `h` to match `(1-x)^(9/2)` modulo `x^4` gives

```text
h = 1 - (9/2)*x + (63/8)*x^2 - (105/16)*x^3 + a*x^4
f = (h^2 + (x - 1)^9)/x^4.
```

Then `P=(1,h(1))` gives a rational divisor class of order `9`, provided the
curve is smooth and `h(1) != 0`.

Code:

```text
code/contact9_family_search.m
```

## Verification

Run:

```text
magma -b mode:=verify code/contact9_family_search.m > data/contact9_verify.txt
```

For the sample values

```text
a = 0, 1, -1, 2, -3/2, 5/4
```

Magma confirmed

```text
marked_order 9
torsion [9]
```

and each sample had a simple certificate at a small prime.

## Direct family search

Run:

```text
magma -b mode:=family height:=30 max_exact:=2000 \
    code/contact9_family_search.m > data/contact9_family_h30.txt
```

Result:

```text
checked 1111
smooth 1110
marked9 1110
exact_tests 1110
hits_ge18 1
TORSION_COUNTS:
  [18] 1
  [9] 1109
```

The one simple `[18]` hit was

```text
a = -7/4
torsion [18]
simple true, pcert 17
```

with integral model

```text
y^2 = 65536*x^5 - 389120*x^4 + 3864576*x^3
      - 4488960*x^2 + 2515968*x - 551936.
```

## Rational-root subfamily

Set

```text
r = 1 - s^2
h(r) = eps*s^9,  eps = +/-1.
```

This determines

```text
a = (eps*s^9 - (1 - (9/2)*r + (63/8)*r^2 - (105/16)*r^3))/r^4
```

and forces `f(r)=0`, giving rational `2`-torsion together with the marked
`9`-torsion.  Thus the root subfamily generically has torsion containing
order `18`.

Finite diagnostics:

```text
magma -b mode:=root_finite target:=36 code/contact9_family_search.m \
    > data/contact9_root_finite_target36.txt
magma -b mode:=root_finite target:=72 code/contact9_family_search.m \
    > data/contact9_root_finite_target72.txt
```

The root chart has no good open points at `p=3` or `p=5`, so those are boundary
primes for this parametrization.  For primes `p >= 7`, there is local room for
both targets.  For example, target `72` has:

```text
p 7  good 2   pass_target 2
p 11 good 10  pass_target 2
p 13 good 16  pass_target 12
p 17 good 24  pass_target 12
p 19 good 24  pass_target 12
```

## Root searches for 36 and 72

Target `36`, height `80`, with all local survivors exact-tested:

```text
magma -b mode:=root height:=80 target:=36 progress_interval:=4000 \
    max_exact:=2000 code/contact9_family_search.m \
    > data/contact9_root_target36_h80_all_exact.txt
```

Result:

```text
checked 15726
smooth 15718
marked9 15718
target_survivors 1212
exact_tests 1212
hits 0
TORSION_COUNTS:
  [18] 1212
```

Target `72`, height `80`:

```text
magma -b mode:=root height:=80 target:=72 progress_interval:=2000 \
    max_exact:=1000 code/contact9_family_search.m \
    > data/contact9_root_target72_h80.txt
```

Result:

```text
checked 15726
smooth 15718
marked9 15718
target_survivors 4
exact_tests 4
hits 0
TORSION_COUNTS:
  [18] 4
```

The target-72 survivors were only the sign pairs

```text
s = +/-30/19
s = +/-23/76
```

with exact torsion `[18]`.

## Conclusion

The contact-9 scaffold is valid and quickly gives simple Jacobians with
`[9]` and `[18]`, but the rational-root subfamily did not produce `[36]` or
`[72]` up to height `80`.  The next logical contact-9 move would be local
analysis at the bad root-chart primes `p=3` and `p=5`, analogous to what was
done for contact-7.  However, based on the height-80 exact results, the
generic root subfamily appears to stop at `[18]`.


## Direct contact-9 plus 5 attempt for `[45]`

I also tested the cleaner coprime odd target

```text
9 + 5 -> 45.
```

This avoids the repeated second-halving obstruction seen in the `[56]` and
`[72]` attempts.  The search uses the direct one-parameter contact-9 family
and filters for the prime-to-`p` part of `45` in `#J(F_p)`.

Code:

```text
code/contact9_plus5_search.m
```

Finite diagnostic:

```text
magma -b mode:=finite code/contact9_plus5_search.m \
    > data/contact9_plus5_finite.txt
```

The direct affine chart is already obstructed at `p=3`:

```text
p 3 total 3 good 1 pass_target 0
```

For `p >= 5` there is local room, for example:

```text
p 5  good 2   pass_target 2
p 7  good 4   pass_target 1
p 11 good 9   pass_target 1
p 13 good 11  pass_target 2
p 17 good 15  pass_target 1
```

Thus any rational `[45]` point in this direct family would have to reduce to
the bad or nonintegral boundary at `p=3`.

Rational height `200`:

```text
magma -b mode:=search height:=200 progress_interval:=5000 max_exact:=1000 \
    code/contact9_plus5_search.m > data/contact9_plus5_h200.txt
```

Result:

```text
checked 48927
smooth 48925
marked9 48925
target_survivors 0
exact_tests 0
hits 0
```

First-kill distribution:

```text
p 3:  12319
p 7:  13781
p 11: 15239
p 13: 4871
p 17: 2149
p 19: 355
p 23: 137
p 29: 57
p 31: 13
p 37: 1
p 41: 3
```

Conclusion: the direct contact-9 plus 5 route is not promising in the affine
height-200 search.  The only plausible continuation would be a focused
`p=3` boundary analysis, but this is now another small-prime boundary problem
rather than a clean high-torsion construction.
