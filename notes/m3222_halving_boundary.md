# M_1(8,2,2): halving the distinguished order-8 class

We worked in the odd model from `paper/NotesAndTodo.tex`, Section
`(16,2,2)`, with symmetric parameters

```text
s = u + v,    p = uv.
```

The odd quintic is

```text
f_{s,p}(X) =
((p-s+1)X^2 + (2-s)X + 1)((s+2)X + 1)
(-X^2 + (ps - s^2 + 2p - s - 2)X - (s^2 - p + s + 1)).
```

The marked point is `Q=(-1,p(s+1))`.  To test whether `Q` is divisible
by 2, use a monic quadratic `a=X^2+A X+B` and require

```text
f_{s,p}(X) - L (X+1) a(X)^2
```

to be a square of a quadratic, where

```text
L = coeff_X^5(f_{s,p}) = s^2 - sp + s - 2p - 2.
```

The boundary/discriminant factors in `(s,p)` are

```text
p,
s + 1,
s + 2,
p - s + 1,
s - p + 1,
s^2 - 4p,
2s^2 + 3s + p + 1,
s^3 - s^2p + s^2 - 4sp - 4p.
```

Finite-field result:

```text
mod 7:  no split open halving-square points
mod 11: no split open halving-square points
```

Thus a rational specialization with good parameter reduction at `7` or
`11` must reduce to one of the boundary residues below.

Split boundary residues with halving-square points mod 7:

```text
(s,p)=(0,0)   on p=0, delta=0, collision=0
(s,p)=(1,1)   on qdisc=0, collision=0
(s,p)=(1,2)   on s-p+1=0, delta=0
(s,p)=(2,6)   on qdisc=0, collision=0
(s,p)=(3,0)   on p=0, qdisc=0
(s,p)=(3,4)   on s-p+1=0, delta=0
(s,p)=(6,0)   on p=0, s+1=0, s-p+1=0, qdisc=0, collision=0
(s,p)=(6,2)   on s+1=0, delta=0
```

Split boundary residues with halving-square points mod 11:

```text
(s,p)=(0,0)   on p=0, delta=0, collision=0
(s,p)=(4,4)   on delta=0
(s,p)=(5,0)   on p=0, qdisc=0
(s,p)=(6,8)   on qdisc=0
(s,p)=(8,1)   on qdisc=0
(s,p)=(8,5)   on delta=0
(s,p)=(10,0)  on p=0, s+1=0, s-p+1=0, qdisc=0, collision=0
(s,p)=(10,3)  on s+1=0, delta=0
```

Rational searches performed:

```text
direct (u,v) search, height 20:
  256938 nonsingular pairs, 26 exact local survivors, 0 halves

symmetric (s,p,A,B) search, height 10:
  16129 base pairs checked
  627 split base pairs
  375 nonsingular curve base pairs
  6048375 (A,B)-tuples checked
  0 local square survivors
  0 exact checks
  0 hits
```

Conclusion: the simple `M_1(8,2,2)` halving route is strongly locally
constrained.  Any rational `[2,2,2,16]` example from this construction
must live on simultaneous `7`- and `11`-adic boundary charts, not on the
ordinary open part of the family.


## First Hensel Lifts

I then lifted the full local system in variables `(s,p,A,B,C,D,E,z)`, where

```text
f_{s,p}(X) - L(X+1)(X^2+AX+B)^2 = (CX^2+DX+E)^2,
z^2 = s^2 - 4p.
```

This keeps both the halving square root and the rational-splitting condition.

### Mod 7 to Mod 49

All `48` split boundary solutions mod `7` lift to mod `49`.

```text
level 1: 48 solutions on 8 base residues
level 2: 34104 solutions on the same 8 base residues
dead at first lift: 0
tangent dimensions: dim 2 for 24, dim 3 for 12, dim 4 for 12
```

Boundary-depth summary modulo `49`:

```text
still zero mod 49 on no boundary factor: 12936
still zero mod 49 on delta=0: 4116
still zero mod 49 on p=0,delta=0,collision=0: 4802
```

The only mod-`7` base residues with lifts that completely leave all boundary factors modulo `49` are

```text
(s,p)=(3,0): 7056 such lifts
(s,p)=(6,0): 5880 such lifts
```

The other six mod-`7` base residues stay on at least one boundary factor to order `7^2`.

### Mod 11 to Mod 121

All `72` split boundary solutions mod `11` lift to mod `121`.

```text
level 1: 72 solutions on 8 base residues
level 2: 226512 solutions on the same 8 base residues
dead at first lift: 0
tangent dimensions: dim 2 for 24, dim 3 for 36, dim 4 for 12
```

Boundary-depth summary modulo `121`:

```text
still zero mod 121 on no boundary factor: 121880
still zero mod 121 on delta=0: 43076
still zero mod 121 on p=0,delta=0,collision=0: 29282
```

The mod-`11` base residues with lifts that completely leave all boundary factors modulo `121` are

```text
(s,p)=(5,0): 57200 such lifts
(s,p)=(6,8): 9680 such lifts
(s,p)=(8,1): 11440 such lifts
(s,p)=(10,0): 43560 such lifts
```

So there is no immediate `7`-adic or `11`-adic obstruction on the boundary cover itself. The likely useful charts are the transverse boundary-exit charts listed above, especially the common `p=0` direction visible at both primes.


## The `p=0` Boundary Chart

On `p=0` the odd quintic factors as

```text
f_{s,0}(X) = (X+1)^2 (X+s^2+s+1)((s-1)X-1)((s+2)X+1).
```

Away from `s=1,-2`, the residual square condition forces `a(-1)=0`, so write

```text
a(X) = (X+1)(X+B),
root(X) = (X+1)(MX+N).
```

The reduced boundary-square equations are

```text
N^2 + s^2 B^2 + s B^2 + s^2 - 2B^2 + s + 1 = 0,
MN + (1/2)s^2B^2 + s^2B + (3/2)s^2 + (1/2)sB^2 + sB
   - B^2 - 2B + (3/2)s + 2 = 0,
M^2 - s^4 - 2s^3 + 2s^2B + s^2 + 2sB + 2s - 4B + 3 = 0.
```

Eliminating `M,N` gives a single plane curve in `(s,B)`.  The polynomial is printed by

```text
magma -b do_boundary_elimination:=true code/m3222_p0_chart_symbolic.m
```

For the transverse chart I used

```text
u=t,  v=s+tV,
A=B+1+tUA,  B=B+tUB,
C=M+tUC,    D=M+N+tUD,    E=N+tUE.
```

This is the corrected first-order chart; the earlier overly rigid choice with fixed `s,B,M,N` is too restrictive.

### CRT-Filtered Search on the Chart

I then searched the transverse chart using the local congruences

```text
p = 77q,
s mod 77 in {10,27,38,76},
A - B - 1 = 0 mod 77.
```

Run:

```text
magma -b height_s:=40 height_q:=15 height_b:=20 height_r:=5 \
  code/m3222_p0_chart_search.m
```

Result:

```text
base_checked 14606
base_split 73
base_curve_ok 64
tuples_checked 1275456
local_pass 0
quartic_square 0
exact_checked 0
hits 0

rejected_at 7  892632
rejected_at 11 119457
rejected_at 13 30718
rejected_at 17 2919
rejected_at 19 98
```

So the `p=0` chart has p-adic tangent directions, but the rational CRT-filtered search found no candidate even surviving the local square filters through prime `19`.


## Direct `[2,2,8] -> [2,2,16]` residue search

I then restarted the halving search with the explicit target of finding a
`[2,2,16]` example, not requiring any additional splitting beyond the
`M_1(8,2,2)` / `[2,2,8]` family.

The direct finite-open test was rerun:

```text
magma code/m3222_finite_halving_sieve.m   > data/m3222_finite_halving_sieve_rerun.txt
```

It confirms the open obstruction:

```text
p=7:  open 12, Q_order8 12, Q_divisible_by_2 0
p=11: open 24, Q_order8 24, Q_divisible_by_2 0
```

Thus any rational `[2,2,16]` example in this model must reduce to boundary at
both `7` and `11`.

I next ran the residue-filtered exact search:

```text
magma -b height:=30 prime_bound:=43   code/m3222_halving_residue_search.m   > data/m3222_halving_residue_h30_p43.txt
```

Summary:

```text
checked       1234321
curve_ok      1225214
residue_pass      318
exact_checked     318
hits                0
```

So through height `30`, every rational-open parameter satisfying the local
halving residue conditions through `p=43` exact-checks as nondivisible by `2`.

To push the residue stage more cheaply, I added a two-stage pipeline:

```text
code/m3222_halving_dump_residues.m
code/m3222_halving_residue_enum.py
code/m3222_halving_exact_file.m
```

The Magma step dumps open finite residue classes where the order-8 class is
divisible by `2`; the Python step scans rational parameters and keeps only
those that are either boundary or allowed at every prime; the final Magma step
exact-tests the survivors.  Validation at height `30`, `p<=43`, reproduced the
same `318` rational-open candidates and no hits.

Using stronger residues through `p=73`, the rational-open necessary-condition
scan at height `50` found no survivors at all:

```text
python3 code/m3222_halving_residue_enum.py   --height 50   --allowed data/m3222_halving_allowed_residues_p73.txt   --out data/m3222_halving_candidates_h50_p73_open_py.txt   --require-rational-open
```

Summary:

```text
height 50
params 3095
checked 9579025
survivors 0
require_rational_open True
```

A height-80 run with the same exact-rational open filter was started but
interrupted because the Python implementation was unoptimized for exact
`Fraction` arithmetic over `61,826,769` parameter pairs.  The height-50 result
is the clean completed bound.

Conclusion: the direct `[2,2,8] -> [2,2,16]` route remains locally possible
only on simultaneous `7`- and `11`-adic boundary, but the residue conditions
through `p=73` already kill every rational-open parameter through height `50`.
The next useful continuation would be an optimized integer/CRT enumerator for
the rational-open boundary conditions, or a component-wise local chart search
at the simultaneous `7`/`11` boundary.
