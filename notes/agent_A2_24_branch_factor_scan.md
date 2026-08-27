# A(2,12) -> A(2,24) branch factor scan

Date: 2026-07-02.

This is the third-pass bounded scan for a better saturated halving branch than
the previous split fiber

```text
p=-5/3, z=1, r=2/3.
```

The accompanying script is

```text
code/agent_A2_24_branch_factor_scan.m
```

It reuses the A(12) chart and the square-quartic halving equations from the
earlier saturated-cover calculation, but enumerates small rational split
fibers instead of fixing one branch.

## Method

For each small rational triple `(p,z,r)` in the A(12) chart:

1. Build `R,Q,ell,F,f=R*F`.
2. Keep smooth fibers where the residual quartic `F` has a rational quadratic
   pairing.  In the completed height-4 run, every accepted `F` had factor
   degrees `[ <2,1>, <2,1> ]`.
3. Build the visible order-12 class `P12=P4+P6`.
4. Generate rational 2-torsion classes from the factorization of `f` and test
   each translated exact order-12 class `D=P12+T`.
5. For `D=[u,v]`, form

```text
ell_D = v + u*(M*x + N),
S     = (ell_D^2 - f)/u
      = s4*x^4 + s3*x^3 + s2*x^2 + s1*x + s0,
E1    = 8*s4^2*s1 - s3*(4*s4*s2 - s3^2),
E0    = 64*s4^3*s0 - (4*s4*s2 - s3^2)^2.
```

Then eliminate `N`, remove factors coming from `s4=0`, and record the
saturated affine factor degrees.

Run used for the table:

```text
magma -b Height:=4 MaxRows:=1000 Progress:=0 \
  code/agent_A2_24_branch_factor_scan.m
```

The script also includes the previous known fiber as a seed even though
`p=-5/3` is outside height 4.

Summary:

```text
checked=11132 split_fibers=29 order12_split_fibers=29
translated_order12_rows=116 low_rows=40 errors=0
```

In every row:

```text
raw resultant degree in M = 32,
s4-boundary contribution = [ <2,8> ],
saturated affine degree = 16,
gcd(E1,E0)=1.
```

So no positive-dimensional quotient component appeared in this box.  The
useful variation is the factorization of the saturated degree-16 fiber.

## Fiber table

In the table, `O/TR` records the two translations by `0` and by the cyclic
2-torsion class `TR=6P12`; `extra` records the two independent extra
2-torsion translations.  Each entry is the saturated affine factor degree list
after removing `s4=0`.

| source | `(p,z,r)` | `O/TR` | `extra` | best |
|---|---:|---:|---:|---:|
| previous seed | `(-5/3, 1, 2/3)` | `[16]` | `[16]` | `16` |
| height 4 | `(-1/2, -4, -4)` | `[16]` | `[8,8]` | `8` |
| height 4 | `(-1/2, -3/2, -3/4)` | `[16]` | `[16]` | `16` |
| height 4 | `(-1/2, -4/3, 2/3)` | `[16]` | `[8,8]` | `8` |
| height 4 | `(-1/2, -4/3, 4)` | `[16]` | `[8,8]` | `8` |
| height 4 | `(-1/2, -1, -1/4)` | `[16]` | `[16]` | `16` |
| height 4 | `(-1/2, -1/2, 2/3)` | `[16]` | `[8,8]` | `8` |
| height 4 | `(-1/2, 1/2, 2/3)` | `[16]` | `[8,8]` | `8` |
| height 4 | `(-1/2, 1, -1/4)` | `[16]` | `[16]` | `16` |
| height 4 | `(-1/2, 4/3, 2/3)` | `[16]` | `[8,8]` | `8` |
| height 4 | `(-1/2, 4/3, 4)` | `[16]` | `[8,8]` | `8` |
| height 4 | `(-1/2, 3/2, -3/4)` | `[16]` | `[16]` | `16` |
| height 4 | `(-1/2, 4, -4)` | `[16]` | `[8,8]` | `8` |
| height 4 | `(-1/3, -1, 4/3)` | `[16]` | `[4,4,8]` | `4` |
| height 4 | `(-1/3, 1, 4/3)` | `[16]` | `[4,4,8]` | `4` |
| height 4 | `(1/3, -1, -4/3)` | `[16]` | `[4,4,8]` | `4` |
| height 4 | `(1/3, 1, -4/3)` | `[16]` | `[4,4,8]` | `4` |
| height 4 | `(1/2, -4, 4)` | `[16]` | `[8,8]` | `8` |
| height 4 | `(1/2, -3/2, 3/4)` | `[16]` | `[16]` | `16` |
| height 4 | `(1/2, -4/3, -4)` | `[16]` | `[8,8]` | `8` |
| height 4 | `(1/2, -4/3, -2/3)` | `[16]` | `[8,8]` | `8` |
| height 4 | `(1/2, -1, 1/4)` | `[16]` | `[16]` | `16` |
| height 4 | `(1/2, -1/2, -2/3)` | `[16]` | `[8,8]` | `8` |
| height 4 | `(1/2, 1/2, -2/3)` | `[16]` | `[8,8]` | `8` |
| height 4 | `(1/2, 1, 1/4)` | `[16]` | `[16]` | `16` |
| height 4 | `(1/2, 4/3, -4)` | `[16]` | `[8,8]` | `8` |
| height 4 | `(1/2, 4/3, -2/3)` | `[16]` | `[8,8]` | `8` |
| height 4 | `(1/2, 3/2, 3/4)` | `[16]` | `[16]` | `16` |
| height 4 | `(1/2, 4, 4)` | `[16]` | `[8,8]` | `8` |

## Conclusion

The previous seed fiber is not representative of the small split fibers:
many extra 2-torsion translations split the saturated degree-16 halving fiber.

Best branches in this bounded scan:

```text
p=-1/3, z=-1, r=4/3,
p=-1/3, z= 1, r=4/3,
p= 1/3, z=-1, r=-4/3,
p= 1/3, z= 1, r=-4/3.
```

For each of these, the two extra translations have saturated affine factor
degrees

```text
[ <4,1>, <4,1>, <8,1> ].
```

No rational factor appeared: the minimum saturated factor degree in the height-4
box is `4`.  Still, these degree-4 branches are a better next target than the
irreducible degree-16 `p=-5/3,z=1,r=2/3` fiber.
