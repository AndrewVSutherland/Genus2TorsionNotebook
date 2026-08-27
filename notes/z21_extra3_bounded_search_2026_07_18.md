# Bounded 3-primary lifts of Lepr\'evost's `[21]` family

## Result

No `[3,21]` or `[63]` specialization survived a rational parameter search to
height `5000` using all usable good reductions through `97`.

The exact box was

```text
t=a/b, gcd(a,b)=1, |a|<=5000, 1<=b<=5000, t not in {0,-1}.
```

It contains `30,401,829` smooth parameters.  There were zero survivors for
both targets, so the strengthened run did not need an expensive rational
`TorsionSubgroup` computation.  Peak RSS was only `25,856 KB`, well below the
hard `300 MB` Magma limit.

Implementation and certified numerical summary:

```text
code/z21_extra3_residue_search.m
data/z21_extra3_residue_h5000_p97.txt
```

## Why the sieve is rigorous inside the box

Let `J_t` be the Jacobian of Lepr\'evost's direct family.  It already has a
marked rational point of order `21`, and hence one rational 3-direction and a
rational point of order `7`.

At every prime `p != 3` where the displayed model has good reduction, rational
3-primary torsion injects into `J_t(F_p)`.  Consequently:

- `[3,21]` requires at least two invariant factors of `J_t(F_p)` divisible by
  `3`, i.e. `dim_F3 J_t(F_p)[3] >= 2`;
- `[63]` requires rational 9-torsion, so at least one finite invariant factor
  must be divisible by `9`.  Combining such a rational 9-class with the
  already marked rational 7-class would give an element of order `63`.

For every prime `5 <= p <= 97`, `p != 3`, the script precomputes the good
residue classes satisfying each necessary condition.  A prime is skipped for
a rational parameter only if its denominator is divisible by `p` or the
displayed reduction is singular; it is never used to reject in those cases.
Thus the sieve has no false negatives in the stated parameter box.

Two useful local snapshots explain its strength:

- modulo `5`, four of five affine parameter residues give singular displayed
  models, while the sole good residue has 3-rank one and no 9-torsion;
- modulo `47`, only `t=0,-1` give singular displayed models, and all 45 good
  residues have 3-rank one.  Fifteen of those residues do have 3-primary
  exponent at least `9`, so the cyclic-`63` test genuinely differs from the
  `[3,21]` rank test.

The full rank and exponent residue sets can be regenerated with

```text
magma -b mode:="finite" prime_bound:=97 MemMB:=300 \
  code/z21_extra3_residue_search.m
```

The bounded search is

```text
magma -b mode:="search" height:=5000 prime_bound:=97 \
  max_exact:=20 MemMB:=300 code/z21_extra3_residue_search.m
```

## Interpretation

This rules out small-height lifts on the direct, generically geometrically
simple Lepr\'evost family.  It does not prove that no exceptional
specialization at larger height exists, and it does not search the entire
moduli locus with rational 21-torsion.

An earlier weaker pass through prime `43` left `342` residue survivors at
height `1000`; exact `TorsionSubgroup` on its first 20 survivors returned
`[21]` every time.  Extending the local tables through `97` removed every
survivor already at height `1000`, and the final height-`5000` pass again left
none.
