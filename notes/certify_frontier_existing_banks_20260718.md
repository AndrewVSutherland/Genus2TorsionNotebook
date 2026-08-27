# Strict simplicity certificates from the existing 2-primary banks

Date: 2026-07-18

The reusable scanner is `code/certify_frontier_existing_banks.m`; its final
test transcript is `results/certify_frontier_existing_banks_20260718.log`.
The run used

```text
magma -b max_bank_rows:=7 max_m244_curves:=5 witnesses_per_curve:=2 \
  code/certify_frontier_existing_banks.m
```

For each candidate, the scanner first calls Magma's exact
`TorsionSubgroup(Jacobian(C))`.  It then searches good primes.  At a successful
prime the characteristic polynomial `chi_p` is irreducible and, for a root
`pi`,

```text
Degree(MinimalPolynomial(pi^n)) = 4  for every n = 2,...,12.
```

The point-count reconstruction of `chi_p` is independently checked against
Magma's built-in `LPolynomial` before the root-power test.

## Results

All three groups have strict witnesses and therefore should no longer be
listed as groups not known to occur on a geometrically simple genus-2
Jacobian over `Q`.

### `[2,2,4,4]` (order 64)

The first row of `paper/scripts_and_data/tor2244.txt` suffices:

```text
(a,b,c,d) = (55944,64800,64935,65160),
C: y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2).
```

Expanded:

```text
y^2 = x^5 + 15791150961*x^4
      + 93063465104931331200*x^3
      + 242416649475931369876277760000*x^2
      + 235275596525928286349936035430400000000*x.
```

Magma gives exact torsion `[2,2,4,4]`.  Two root-power witnesses are

```text
p=67: chi = T^4 + 8*T^3 + 22*T^2 + 536*T + 4489,
p=83: chi = T^4 - 8*T^3 - 10*T^2 - 664*T + 6889.
```

For each polynomial, the degree list for `n=2,...,12` is
`[4,4,4,4,4,4,4,4,4,4,4]`.

A much smaller, publication-friendly representative occurs at row 26629:

```text
(a,b,c,d) = (36,57,64,132),
C: y^2 = x(x+1296)(x+3249)(x+4096)(x+17424)
         = x^5 + 26065*x^4 + 173387808*x^3
           + 414985109760*x^2 + 300512487407616*x.
```

It again has exact torsion `[2,2,4,4]`.  Its two witnesses are

```text
p=37: chi = T^4 - 4*T^3 - 2*T^2 - 148*T + 1369,
p=47: chi = T^4 + 8*T^3 + 30*T^2 + 376*T + 2209,
```

with degree list `[4,4,4,4,4,4,4,4,4,4,4]` in both cases.  This is the
preferred representative for exposition; row 1 remains an independent
baseline.

### `[2,2,2,8]` (order 64)

Row 3 of `paper/scripts_and_data/tor2228.txt` suffices:

```text
(a,b,c,d) = (1,55,99,125),
C: y^2 = x(x+1)(x+3025)(x+9801)(x+15625).
```

Expanded:

```text
y^2 = x^5 + 28452*x^4 + 230082726*x^3
      + 463480444900*x^2 + 463250390625*x.
```

Magma gives exact torsion `[2,2,2,8]`.  Two root-power witnesses are

```text
p=41: chi = T^4 + 4*T^3 + 6*T^2 + 164*T + 1681,
p=47: chi = T^4 + 4*T^3 + 30*T^2 + 188*T + 2209.
```

Both degree lists are `[4,4,4,4,4,4,4,4,4,4,4]`.  Rows 1 and 2 are
projectively the same tuple `(4,11,16,44)` and did not yield two witnesses in
the scanner's prime list.  The direct bielliptic diagnostics below explain
that failure.

In fact, row 3 is the smallest geometrically simple representative in this
bank when primitive tuples are ordered by `max(|a|,|b|,|c|,|d|)`.  The four
smaller primitive tuples have heights `44`, `46`, `117`, and `119`.  Direct
diagnostics on representatives at rows 1, 5, 6, and 14 give automorphism-group
order `4` and two rational degree-2 elliptic subcovers in every case, so all
four are bielliptic.  The next primitive height is `125`, attained by the
successful tuple `(1,55,99,125)` above.

### `[2,4,4]` (order 32)

The very first regenerated curve from the `M(2,4,4)` pilot has

```text
s=-4, t=-3, u=-4, v=-644/799,
C: y^2 = x(x+16)(x+(644/799)^2)(x^2+x+16).
```

Multiplying the right-hand side by the rational square `638401^2` gives the
integral `Q`-isomorphic model used by `TorsionSubgroup`:

```text
y^2 = 407555836801*x^5 + 7193217102753*x^4
      + 17542840688944*x^3 + 112806866289408*x^2
      + 67780576546816*x.
```

Magma gives exact torsion `[2,4,4]`.  Two root-power witnesses are

```text
p=97:  chi = T^4 + 8*T^3 + 46*T^2 + 776*T + 9409,
p=103: chi = T^4 - 4*T^3 + 46*T^2 - 412*T + 10609.
```

Both degree lists are `[4,4,4,4,4,4,4,4,4,4,4]`.  In addition, Magma finds
automorphism-group order `2` and zero rational degree-2 elliptic subcovers.
Most importantly, the root-power witnesses themselves rule out a geometric
elliptic splitting.  Thus the `M(2,4,4)` pilot family is not geometrically
split by construction.

## Consequence for the frontier list

`[2,2,4,4]`, `[2,2,2,8]`, and `[2,4,4]` are realized exactly on curves
passing the project's strict geometric-simplicity certificate.  Any proposed
top-ten list of groups still unknown in the geometrically simple case must
remove all three.
