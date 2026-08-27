# R = -25/4 S_B Descent Followup

Date: 2026-07-02

Superseded correction: this note missed the `609/256` scale factor in
the reduced fiber coordinate and therefore describes a spurious genus-5
pullback.  The corrected genus-3 `V_4` certificate is in
`agent_m18_416_R25_4_SB_v4_certificate.md`.

This follows the newest descent suggestion in
`agent_m18_416_p7_blowup_notes.md`: certify the obstruction on the rich
ELS fiber `R = -25/4`, or find a large point.

## Specialized S_B cover

The specialized reduced `S_B` condition is

```text
y^4 - 4 alpha_B(m) y^2 + d_B(m) = 0

alpha_B(m) =
  (29/100*m^4 - 16907/64*m^2 + 268888725/4096)/m^2

d_B(m) =
  (-841/100*m^4 + 2190805/256*m^2 - 7797773025/4096)/m^2
```

The helper identity is exact:

```text
alpha_B^2 - d_B/4 = G*h_B^2
h_B = (2/525*m^4 - 29/8*m^2 + 441525/512)/m^2
```

After clearing denominators with `Y = m*y`, the plane model is

```text
-861184*m^6 - 118784*m^4*Y^2 + 876322000*m^4
+ 108204800*m^2*Y^2 - 194944325625*m^2
+ 102400*Y^4 - 26888872500*Y^2 = 0.
```

Pulling this cover back to the elliptic fiber gives a projective curve of
degree 8 and genus 5.  Magma does not convert it to a genus-one model.

## Quotient check

The visible quotient by `m -> -m` and `Y -> -Y`, with `X = m^2` and
`Z = Y^2`, is

```text
-861184*X^3 - 118784*X^2*Z + 876322000*X^2
+ 108204800*X*Z - 194944325625*X
+ 102400*Z^2 - 26888872500*Z = 0.
```

Its discriminant in `Z` is

```text
1048576*X^4 - 1884160000*X^3 + 1318168320000*X^2
- 426513150000000*X + 53731529750390625
```

with monic factorization, suppressing the leading rational unit,

```text
(X - 15225/32)^2 * (X^2 - 13525/16*X + 231800625/1024)
```

Equivalently, `disc_Z/Gamma(X)` is a rational square.  Thus this easy
quotient only recovers the old fiber relation and does not supply the
missing descent certificate.

## Ordinary elliptic two-descent check

For the minimal elliptic fiber

```text
E: y^2 = x^3 - x^2 - 468441*x + 122191641
MW invariants: [2, 2, 0, 0, 0]
free generators:
  (545, -5336)
  (-489, -15300)
  (8289/25, 224112/125)
torsion generators:
  (361, 0), (429, 0)
```

Magma's ordinary elliptic two-descent, after removing torsion and all
three free generators, returns no residual covers:

```text
TwoDescent residual num_covers = 0
MordellWeilShaInformation = [3, 3]
```

This corrects the interpretation in the previous note: the obstruction
is global and descent-flavored, but not an unaccounted ordinary
`Sha(E)[2]` class for this elliptic fiber.  The actual remaining object
is the genus-5 square-condition cover above.

## Larger point search

The existing Mordell-Weil enumeration script was pushed from the saved
`N = 4` run to `N = 5`:

```sh
magma -b Rnum:=-25 Rden:=4 mex_num:=-25 mex_den:=2 N:=5 \
  code/agent_m18_416_els_mw_deep.m
```

Result:

```text
MW points visited = 5323
distinct m tested = 2534
skipped(height) = 610
passes = 0
```

No large point was found in this enlarged box.

## Current next target

The next certification step should treat the genus-5 cover directly:
compute a 2-cover/Jacobian or Prym-style descent for the square condition
`Z = Y^2` over the rank-3 fiber, rather than running ordinary
two-descent on `E` alone.
