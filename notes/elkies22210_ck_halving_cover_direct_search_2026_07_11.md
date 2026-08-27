# Direct searches on the Elkies Clebsch--Klein halving covers

This note records searches on the covers themselves, rather than a further
enumeration of bounded integral points on the Clebsch--Klein base.  Throughout

```text
sum r_i = sum r_i^3 = 0,
a_i = r_i^2,
```

and all open conditions (`r_i != 0` and pairwise distinct `a_i`) are enforced.
All arithmetic in the search programs is exact integer/rational arithmetic.

## Orbit A: `{0,a_1}`

Normalize `r_1=1`.  The four square conditions are parametrized by

```text
r_j = (1-t_j^2)/(1+t_j^2),
1-r_j^2 = (2*t_j/(1+t_j^2))^2,       j=2,...,5.
```

If `x=r_2` and `y=r_3`, eliminate `u=r_4` and `v=r_5` from the two CK
equations.  They satisfy

```text
u+v = S = -1-x-y,
u*v = P = (x+y)*(1+x)*(1+y)/(1+x+y).
```

The exceptional denominator `1+x+y=0` gives no open solution: then the
linear equation gives `u+v=0`, while the cubic equation forces `x*y=0`.
Off that boundary, `u,v` are rational exactly when `D=S^2-4P` is a rational
square.  The program then tests `1-u^2` and `1-v^2` for rational squareness
and independently rechecks both CK equations.

This is implemented in

```text
code/elkies22210_ck_orbitA_direct.cpp
```

The height `B` means that `t_2=m_2/n_2` and `t_3=m_3/n_3` are reduced with
positive numerator and denominator at most `B`.  The recovered `t_4,t_5`
have **no height bound**.  Because the four unmarked coordinates are
symmetric, the run is complete for any cover point for which some two of the
four unit-circle parameters have height at most `B`.

The exact `B=150` result is

```text
circle_values                 13714
unordered (x,y) pairs         94030041
distinct input squares        94023184
positive discriminant         65774803
square discriminant           3
rational nonzero (u,v)        3
1-u^2 square                  0
1-v^2 square                  0
smooth cover points           0
```

The first two rational CK completions already occur by `B=50`; the third
appears by `B=150`:

```text
(t2,t3)=(41/2,16/9)
r=[1,-1677/1685,-175/337,2871/5729,384/28645]
1-u^2=24578800/32821441,  1-v^2=820388569/820536025

(t2,t3)=(11/8,17/9)
r=[1,-57/185,-104/185,252/185,-276/185]
1-u^2=-29279/34225,       1-v^2=-41951/34225

(t2,t3)=(31/24,27/121)
r=[1,-385/1537,6956/7685,-12584/26129,-153252/130645]
1-u^2=524367585/682724641,
1-v^2=-6418059479/17068116025
```

In the first row both remaining values are positive nonsquares.  In the
second they are negative.  In the third, one is a positive nonsquare and the
other is negative.  Thus none is a cover point.

There is also a useful symmetric form of the same cover equations.  Off the
boundary, the two CK equations and four circle parametrizations reduce to

```text
sum_j 1/(1+t_j^2) = 3/2,
prod_j (1-t_j^2) = 16,               j=2,...,5.
```

As an independent check, the pairwise meet-in-the-middle implementation

```text
code/elkies22210_orbit01_unit_circle_search.py
```

hashes the sum and product invariants.  At parameter height `50` it found

```text
1546 t-values,
1193512 admissible pairs,
0 complementary pair joins,
0 points.
```

That run bounds all four `t`-parameters by `50`; the elimination run above is
stronger at its stated bound because only two parameters are bounded and the
other two are recovered with no height cap.  Further local and
meet-in-the-middle details are cross-recorded in

```text
notes/elkies22210_source_halving_local_and_search_2026_07_11.md.
```

## Orbit B: `{a_1,a_2}`

For complementary indices `j=3,4,5`, division by `G_0` gives the exact
simplification

```text
G_j/G_0 = -(a_j-a_2)/(a_1-a_j).
```

Consequently the four specified quantities are squares if and only if
`G_0` is a square and all three displayed ratios are squares.  With
`r_1=1`, `r_2=q`, and ratio square root `w_j`, each remaining root lies on
the genus-one quartic

```text
w_j^2 = (q^2-r_j^2)/(1-r_j^2),
r_j^2 = (q^2-w_j^2)/(1-w_j^2).
```

This shows why Orbit B does not collapse to four independent conic
parametrizations: for a fixed marked root `q`, each complementary root lies
on a full-2-torsion elliptic quartic.

There is a sharp real sieve.  Every complementary `a_j` must lie outside the
interval with endpoints `a_1,a_2`, so the marked pair is adjacent after the
five positive `a_i` are sorted.  The sign of `G_0` then leaves only marked
positions `{1,2}` or `{3,4}`; positions `{2,3}` and `{4,5}` are impossible.
The direct code independently checks this ordering assertion on every row
where all four `G` values are positive.

For the bounded exact search, again normalize `r_1=1`, choose arbitrary
reduced rationals

```text
x=r_2,  y=r_3,
height(x),height(y) <= B,
```

and use the same formulas for `S,P,D` to solve `u=r_4,v=r_5` with no height
bound.  On every smooth rational CK completion the program evaluates the
four original expressions `G_0,G_3,G_4,G_5` literally, first applies the
real sign sieve, and then performs four exact rational-square tests.

The implementation is

```text
code/elkies22210_ck_orbitB_direct.cpp
```

The exact search has now reached `B=100`:

```text
rational values x or y        12172
ordered (x,y) pairs           148157584
smooth rational CK tuples     15004
all four G values positive    3682
individual square counts
  (G0,G3,G4,G5)               (3,1,1,2)
all four squares              0
```

The first new partial row beyond `B=50` is

```text
r=[1,2/15,74/45,-289/225,-112/75],
(G0,G3,G4,G5) square mask=(0,0,1,0).
```

No partial row has more than one of the four required square conditions.
As a check on the real sieve, at `B=20` the 147 positive rows split as 76
with marked ranks `{1,2}`, 71 with ranks `{3,4}`, and 0 in any other rank
pair.

## Full rational two-parameter chart of the CK surface

There is also a single rational chart covering the entire **labelled smooth
CK open**, not merely a dense subset.  Put

```text
r1=x+y,  r2=x-y,  r3=z+w,  r4=z-w,  r5=-2*(x+z).
```

After scaling `x=1` and writing `z=t`, the cubic equation is

```text
Y^2+t*W^2=(1+t)*(t^2+3*t+1).
```

The line

```text
W+t+2=m*(Y+1)
```

through the section `(-1,-t-2)` gives the projective point

```text
R1 = 1+t*(t+2)*m,
R2 = t*m*(m-t-2),
R3 = -1+m+t*(t+1)*m^2,
R4 = 1+t-m-t*m^2,
R5 = -(1+t)*(1+t*m^2).
```

Conversely, a labelled CK point has

```text
t = (r3+r4)/(r1+r2),
Y = (r1-r2)/(r1+r2),
W = (r3-r4)/(r1+r2),
m = (W+t+2)/(Y+1).
```

There is no missing smooth locus.  If `r1+r2=0`, then `r1^2=r2^2`
(and the CK equations also force one of `r3,r4,r5` to vanish).  Moreover
`Y+1=2*r1/(r1+r2)`, so its vanishing is `r1=0`.  Both cases lie on the
excluded discriminant boundary.  The universal Magma script verifies
symbolically both CK identities and that this inverse recovers the five
coordinates projectively modulo the conic equation.

For reduced `t=a/b` and `m=c/d`, multiplication by `b^2*d^2` gives exactly
the five integer expressions used by

```text
code/elkies22210_ck_rational_param_search.cpp.
```

Dividing their common gcd does not change the square tests: orbit-I
radicands are homogeneous of degree `2`, and orbit-II radicands of degree
`6`, so projective rescaling changes each by a rational square.  The code
performs the coordinate CK checks and orbit-I differences exactly in signed
`128`-bit arithmetic, and the three prime-power filters in modular
`128`-bit arithmetic.  It constructs GMP degree-six radicands only for the
final CRT survivors, where exact square decisions use GMP's perfect-square
predicate.  The modular tables are necessary filters only and retain zero
residues.  Since the tuple is already checked smooth, rejecting zero as well
as negative exact radicands is correct.
At the advertised maximum `B<=500`, each cleared coordinate has absolute
value at most `4*B^4 <= 2.5*10^11`; the coordinate and cubic CK checks are
therefore safely inside signed `128`-bit range before conversion to GMP.

At rational parameter height `B=100`, the completed full-chart search gives

```text
rational t- or m-values       12175
parameter pairs               148230625
smooth labelled CK tuples     148145150
orbit-I hits                  0
orbit-II hits                 0
```

The same run applies the compatible necessary square-residue filters to both
real-eligible orbit-II pairs on every tuple:

```text
real-eligible marked pairs    296290300
all four squares mod 11^3      11626790
also squares mod 19^2            543370
also squares mod 23^2             17209
exact four-square hits                 0
```

All `17,209` final congruence survivors have exact square mask `0`; the modular
tests therefore introduce no claimed false positive as an exact cover point.

The chart is globally complete on the labelled smooth open, while this
negative computation is bounded only in its two inverse parameters `t,m`;
the reported count is for labelled parameter tuples and does not quotient
relabeling symmetries of the resulting genus-2 curve.

## Conclusion

Neither direct chart produced a smooth rational cover point.  The negative
statements are bounded, not proofs of global insolubility.  Orbit A is very
sparse already at the rational-CK-completion stage.  Orbit B has many more
rational CK completions, but its three elliptic cross-ratio constraints are
severe: among 3,682 real-compatible rows at `B=100`, only seven individual
`G` square events occurred and none coincided.  The complete rational CK
chart adds more than 148 million smooth labelled parameter tuples without a
hit.
