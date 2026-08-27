# Orbit-12 quotient-fiber search

Date: 2026-07-11.

This search imposes the first exact orbit-12 radicand before looking for
Clebsch--Klein points.  Normalize the marked pair to `r1=1`, `r2=q` and put

```text
s = r3*r4+r3*r5+r4*r5.
```

After quotienting by permutations of the three complementary roots, the
first-radicand cover is

```text
v^2 = q*(s-q)*(2-q^2+(q+2)*s).
```

For fixed `s` it has Weierstrass model

```text
Y^2 = X^3+(s^2-2*s-2)*X^2
          -4*s^2*(s+1)*X+4*s^2*(s+1)^2,
q = 2*s*(s+1)/X.
```

The quotient point lifts to rational labelled complementary roots exactly
when

```text
z^3+(1+q)*z^2+s*z-(1+q)*(q-s)
```

splits completely over `Q`.  The implementation

```text
code/elkies22210_orbit12_rank1_fiber_multiples.m
```

factors this cubic exactly and, for every split lift, rechecks both CK
identities, smoothness, and all four literal Stoll--Zarhin radicands.

## Positive control

The known first-radicand lift

```text
s = 59/49,
q = 8/7,
(r3,r4,r5) = (-9/7,-5/7,-1/7)
```

is built into the script as a regression test.  It maps to

```text
(X,Y) = (1593/343,38232/16807)
```

and its exact radicands have square mask

```text
[true,false,false,false].
```

Thus the test certifies the quotient map and split-cubic reconstruction,
while deliberately not being a full-cover hit.

## Completed search

The completed command was

```text
magma -b integer_bound:=30 s_height:=0 multiple_bound:=50 \
  code/elkies22210_orbit12_rank1_fiber_multiples.m
```

It used all nonsingular integer fibers `-30 <= s <= 30`, together with the
positive-control fiber `s=59/49`.  On every positive-rank fiber, for each
free generator `P_i`, it tested

```text
n*P_i+T,  1 <= n <= 50,
```

for every one of the four rational torsion points `T`.  On rank-two fibers
this is a search along both generator axes; mixed combinations
`n*P_1+m*P_2` are not included.

The exact results were

```text
fibers considered                         60
proven rank-zero fibers                   29
proven positive-rank fibers               31
ambiguous rank fibers                      0
distinct finite quotient q values      6,800
cubics with a rational root                3
cubics split completely                    1
smooth split lifts                         1
full four-radicand hits                    0
maximum exact square-mask weight           1
```

The sole completely split cubic was the positive control.  In particular,
the integer fibers contributed `6,400` quotient values and no completely
split cubic.  With `SetSeed(1)` fixing Magma's generator choices, the run
used `24.660` CPU seconds (`26.96` wall-clock seconds) in the final logged
run on the current workspace machine.

A second completed run searched genuinely rational fibers:

```text
magma -b integer_bound:=0 s_height:=8 multiple_bound:=20 \
  code/elkies22210_orbit12_rank1_fiber_multiples.m
```

Here `s_height=8` means all reduced nonintegral `s=a/b` with
`|a|,b <= 8`, again together with the positive control.  Its exact totals
were

```text
fibers considered                         71
proven rank-zero fibers                   42
proven positive-rank fibers               29
ambiguous rank fibers                      0
distinct finite quotient q values      2,640
cubics with a rational root                4
cubics split completely                    1
smooth split lifts                         1
full four-radicand hits                    0
maximum exact square-mask weight           1
CPU seconds                            3.030
```

Again the unique split lift was the positive control.  Thus none of the
height-eight rational `s` fibers produced a new labelled CK point even on
the first-radicand cover.

The reducible-but-not-split cubics give a short list of near misses:

```text
s=-17,  q=1:    (z+4)*(z^2-2*z-9)        [q=1 is boundary]
s=-8,   q=-7/2: (z+5/2)*(z^2-5*z+9/2)
s=-7/3, q=-8/3: (z+1/3)*(z^2-2*z-5/3)
s=-2/3, q=-1/6: (z+1/2)*(z^2+z/3-5/6)
s=7/2,  q=3/2:  (z+2)*(z^2+z/2+5/2)
```

In every nonboundary row the quadratic factor is irreducible over `Q`, so
none yields three rational complementary roots.

## Exceptional elliptic fiber

The branch degeneration `q=-2` is stronger than a negative bounded search.
Its first-radicand condition is

```text
y^2 = -x*(x-1)*(x+1)*(x-2),
```

whose elliptic model is

```text
y^2 = x^3+x^2-4*x-4.
```

Magma proves rank zero and torsion `[2,2]`; pulling the entire group back
gives only `x=-1,0,1,2`, all boundary points.  Hence this distinguished
fiber has no smooth point even before the remaining three radicands are
imposed.  The exact certification is in

```text
code/elkies22210_orbit12_first_radicand_geometry.m.
```

## Next fixed-ratio target

The shortest marked ratio in the selected compatible `11^3,19^3,23^2`
CRT disk is

```text
q = -38357/22386,
q = 242 mod 1331,
q = 248 mod 6859,
q = 392 mod 529.
```

This is the natural next fixed-`q` genus-two descent target.  The present
search did not test that fiber.
