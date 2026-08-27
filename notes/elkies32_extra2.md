# Trying for `[2,32]` in the Reconstructed Elkies `[32]` Family

Goal: start with the reconstructed Elkies `[32]` family and force one more
independent rational `2`-torsion point.  The built-in rational root `r` gives

```text
(r,0) - infinity = 16*((0,1)-infinity),
```

so it is not independent.  A distinct rational finite Weierstrass root `s != r`
would give the desired extra independent `2`-torsion, hence a subgroup
`Z/2 x Z/32` when the marked class has exact order `32`.

## Corrected Scale

While setting this up, I found that the factorized `p=Pnum/Pden` formula copied
from `data/elkies32_reconstruct_conditions.txt` was missing the rational unit
from Sage's factorization output.  At the printed point the raw ratio gives
`-360/11`, while the correct value is `-1440/11`.  The correct formula is

```text
p = 4*Pnum(z,r)/Pden(z,r).
```

I corrected the `[2,32]` scripts and also regenerated the order-64 condition and
search data using this scale.

## Algebraic Condition

Use the reconstructed contact model

```text
y^2 = (a*x^3 + b*x^2 + c*x + 1)^2 - a^2*x^5*(x+1),
```

with

```text
a = z*p,
c = (-z^4*p + 4*z^3*p - 8*z^2*p + 8*z^2 + 8*z*p - 16)/(4*z^2),
b = (-z^4*p + 8*z^3*p - 12*z^2*p + 4*z^2 + 8*z*p - 16)/(4*z^2),
p = 4*Pnum(z,r)/Pden(z,r).
```

The `[32]` base is still the genus-0 plane curve `C32(z,r)=0`.  Let
`Hextra(z,r,s)` be the numerator of `f(s)` after substituting these formulas.
Then the extra-root cover is

```text
C32(z,r) = 0,
Hextra(z,r,s) = 0,
s != r,
```

with the diagonal component `s=r` removed/saturated and with boundary
denominators excluded.

The derived polynomial has size

```text
Hextra degrees (z,r,s): (18,12,5)
Hextra terms: 1234
```

The diagonal check is now consistent: `Hextra(z,r,r)` reduces to `0` modulo
`C32(z,r)`.

## Local Information

Modulo `13`, the good parameter residues have no extra root.  More explicitly,
for the chosen degree-8 genus-0 parameter, the good residues are

```text
1, 4, 6, 12
```

and none have an extra root.  The only possible residues for a rational
candidate to hide modulo `13` are boundary/bad-reduction residues:

```text
W=0:     0, 2, 3, 5, 8, 9, 11
PDen=0:  7, 10, infinity
```

So there is a strong local obstruction away from the `p=13` boundary.

## Search

The search used the genus-0 parameter of the reconstructed `[32]` component.
It applied extra-root residue filters modulo

```text
19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73
```

while keeping `p=11` and `p=13` diagnostic-only.  The exact stage factored the
specialized quintic over `Q` and required a rational root distinct from `r`.

Corrected height-5000 run:

```text
total reduced parameters scanned: 30,401,831
passed residue sieve: 149
nonsingular exact checks: 149
hits: 0
```

All 149 exact survivors had no rational extra root.

## Files

- `code/elkies32_extra2_conditions.py` derives the extra-root cover.
- `data/elkies32_extra2_conditions.txt` records the cover size and diagonal check.
- `code/elkies32_extra2_param_search.m` performs the residue sieve and exact factor checks.
- `data/elkies32_extra2_param_search_B5000.txt` is the corrected height-5000 search log.

## Conclusion

No `[2,32]` specialization was found in this first targeted search.  The useful
new information is structural: good reduction at `p=13` blocks the extra-root
condition, so the next serious step is to analyze the `p=13` boundary residues
`W=0` and `PDen=0`, rather than simply pushing the height bound.
