# Orbit-12 Clebsch--Klein CRT/lattice pilot

Date: 2026-07-13.

## Scope

This is the bounded lattice pilot recommended after the local and height-box
search in

```text
notes/elkies22210_source_halving_local_and_search_2026_07_11.md.
```

It uses one selected product of the three displayed open Hensel disks at

```text
11^3 = 1331,  19^2 = 361,  23^2 = 529
```

for the fixed orbit-12 marked pair `{r_1^2,r_2^2}`.  It does **not** enumerate
all local cover states and it does not claim a global obstruction.

The complete rational Clebsch--Klein chart is

```text
R1 = 1+t(t+2)m,
R2 = tm(m-t-2),
R3 = -1+m+t(t+1)m^2,
R4 = 1+t-m-tm^2,
R5 = -(1+t)(1+tm^2).
```

Its inverse on the labelled smooth open simplifies to

```text
t = (r3+r4)/(r1+r2),
m = (r1+r2+r3)/r1.
```

The exact search is

```text
code/elkies22210_orbit12_ck_crt_lattice.py
```

and the output of the final bounded run is

```text
data/elkies22210_orbit12_ck_crt_lattice_r1000000.txt.
```

## Integral-chart correction

Using the original coordinate order gives the tempting combined residues

```text
t = 147365295,  m = 195384016  (mod 254179739).
```

Those are only chart residues, not a faithful product Hensel disk.  At the
`19^2` seed both the actual first homogeneous coordinate and the related
chart factor satisfy

```text
R1 = 1+t(t+2)m = 228 mod 361,  v_19(R1) = 1,
1+t*m^2 = 95 mod 361,
v_19(1+t*m^2) = 1.
```

Here `1+t*m^2` is the proportionality factor appearing in the inverse-chart
identity; it is not the coordinate used to normalize the program's CK
tuple.  Since the homogeneous tuple has a common factor `19` and its actual
normalizing coordinate `R1` also has valuation one, primitive normalization
divides by `19` and loses one digit.  In particular, substituting the two
residues modulo `361` does not recover the displayed normalized seed modulo
`361`.  Calling the resulting CRT vector a Hensel point would be incorrect.

Swapping `r4` and `r5`, i.e. using zero-based chart permutation

```text
(0,1,2,4,3),
```

keeps the marked pair in positions `1,2` and makes both `R1` and
`1+t*m^2` units at all three seeds.  The program asserts those two unit
conditions, exact projective recovery of every seed, and the four local
square conditions before searching.

The corrected local and combined parameter residues are:

| modulus | `t` | `m` | `R1 mod p` | chart factor `1+tm^2 mod p` |
|---:|---:|---:|---:|---:|
| `1331` | `837` | `1202` | `10 mod 11` | `10 mod 11` |
| `361` | `273` | `347` | `12 mod 19` | `5 mod 19` |
| `529` | `248` | `511` | `7 mod 23` | `14 mod 23` |
| CRT product `254179739` | `49651130` | `195384016` | -- | -- |

Thus a rational parameter `a/b` in the selected `t`-disk is a primitive
vector with positive unit denominator in

```text
a - 49651130*b = 0 mod 254179739,
```

and similarly `m=c/d` satisfies

```text
c - 195384016*d = 0 mod 254179739.
```

## Short lattices

The two determinant-`254179739` numerator--denominator lattices have
Gauss-reduced bases

```text
t: (11517,7766), (11032,-14631),
m: (-11435,2957), (4827,20980).
```

The reduced bases describe the congruence lattices; a basis vector need not
be primitive, so it is not automatically a reduced rational parameter in
the selected disk.  The program discards such vectors and enumerates every
primitive lattice vector `(numerator,denominator)` of Euclidean norm at most
the requested radius.  The shortest retained vectors are

```text
t = -11032/14631,  m = -11435/2957.
```

They are genuine rational chart parameters, but are not exact cover points.
Since the radius is below half the CRT modulus, a denominator has at most one
centered numerator.  This gives a streaming exact search with no large
lattice box in memory.

For each pair it:

1. constructs the primitive integral CK tuple from the homogeneous chart;
2. verifies `sum r_i = sum r_i^3 = 0` exactly;
3. removes the CK boundary;
4. evaluates the four exact Stoll--Zarhin orbit-12 radicands; and
5. applies exact integer-square tests, after harmless modular prefilters.

## Bounded result

The final command was run with hard wall-time and address-space caps:

```text
timeout 300s prlimit --as=8589934592 -- \
  python3 code/elkies22210_orbit12_ck_crt_lattice.py --radius 1000000
```

It completed within the cap and gave

```text
t short vectors                    3,152
m short vectors                    3,138
parameter pairs                9,890,976
smooth exact CK points         9,890,976
all-four-positive radicands    1,799,807
points with any exact square           0
exact four-square cover points         0
```

The `9,890,976` points are actual rational CK points.  They are not Hensel
approximations.  Conversely, the original coordinatewise CRT tuple and the
short lattice basis vectors by themselves are only seeds; none is an exact
halving point.

## Decision

**No-go for enlarging this same selected three-disk CRT coset by another
unstructured radius increase.**  Radius `10^6` already gives nearly ten
million smooth exact CK points, and none has even one of the four required
rational square radicands.  The complete source-halving route remains open,
because the three primes have many other compatible Hensel disks.  A future
restart should first organize those disks into finitely many integral chart
cosets or derive global component geometry; it should not extend this one
coset or rerun a blind projective height box.
