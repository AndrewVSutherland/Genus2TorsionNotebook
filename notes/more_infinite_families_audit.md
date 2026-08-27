# More Infinite Families Audit

Goal: after reconstructing the Elkies `[32]` family, check whether the current
repository already implies additional infinite simple torsion families that
were missing from the inventory.  I did not count aliases written in
non-invariant notation; all groups below are normalized to invariant factors.

## Method

For each established simple-backed source group, I enumerated subgroup types in
Sage:

```text
[2,12], [2,20], [2,2,2,10], [4,8], [2,4,8],
[30], [32], [22], [23], [20], [18], [14].
```

Sage reports elementary factors, so I normalized them to invariant factors by
combining coprime primary parts.  Thus, for example:

```text
[2,3]       -> [6]
[4,5]       -> [20]
[2,2,5]     -> [2,10]
[2,3,5]     -> [30]
```

## Result

After normalization, every inherited subgroup row was already present in
`notes/infinite_families_inventory.md` except:

```text
[16]
```

This comes from the reconstructed Elkies `[32]` family.  In the contact model

```text
y^2 = (a*x^3 + b*x^2 + c*x + 1)^2 - a^2*x^5*(x+1),
```

the marked class `D=(0,1)-infinity` has order `32` generically, so `2D` gives a
cyclic order-`16` subgroup on the same one-parameter component.

The `[32]` row itself is now local, not literature-only.  The reconstructed
component is the genus-0 plane curve in
`data/elkies32_reconstruct_conditions.txt`, and the printed Elkies member has
a local `Q`-simplicity certificate:

```text
L_7(T) = 49*T^4 + 7*T^3 + 6*T^2 + T + 1,
```

which is irreducible.

## Conclusion

The only new inventory row produced by this pass is the inherited `[16]`
family.  This is real and useful for completeness, but it is not a new
large-torsion construction.  The meaningful constructive work remains:

- find a smaller rational parametrization of the reconstructed `[32]` genus-0
  component;
- or build a genuinely new native family, rather than adding more inherited
  subgroup rows.
