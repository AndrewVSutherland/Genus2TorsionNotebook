# Rational Richelot component of the order-96 `[2,2,2,12]` record

Date: 2026-07-18

The reproducible driver is

```text
code/record_22212_richelot_bfs.m
```

Run from the repository root with

```text
magma -b max_depth:=3 \
  log_file:=results/record_22212_richelot_bfs_depth3.log \
  code/record_22212_richelot_bfs.m
```

The full depth-3 output, including every normalized curve and its absolute
`G2Invariants`, is in

```text
results/record_22212_richelot_bfs_depth3.log
```

The shorter default depth-2 run is in

```text
results/record_22212_richelot_bfs.log
```

## Source and factorization

For the collaborator's generalized model

\[
 y^2+(1+x^2)y=f_0(x),
\]

completing the square gives exactly the reported simplified sextic:

\[
\begin{aligned}
(1+x^2)^2+4f_0(x)={}&3027600x^6+2950382280x^5+602288814361x^4\\
&-63417934304484x^3-2122595910966478x^2\\
&+128056619498204124x+3322970988364151397.
\end{aligned}
\]

It factors over `Q` as

\[
\begin{aligned}
&(x-123109/1740)(x-59)(x+299/12)\\
&\qquad\cdot(x+39)(x+519)(x+75593/145),
\end{aligned}
\]

times the square leading coefficient `3027600 = 1740^2`.  In the ordering
used by the script, the six roots are

\[
(r_1,\ldots,r_6)=
(-75593/145,-519,-39,-299/12,59,123109/1740).
\]

Exact Magma computation gives

\[
J(\mathbf Q)_{\rm tors}\simeq
(\mathbf Z/2)^3\oplus\mathbf Z/12,
\qquad \#J(\mathbf Q)_{\rm tors}=96.
\]

## Source divisibility by 2

The four invariant-factor generators of orders `2,2,2,12` are individually
not divisible by 2 over `Q`; in particular, the displayed order-12 generator
does not lift to order 24.

Among the 15 nonzero rational classes in `J[2]`, exactly one is divisible by
2.  It is represented in Mumford form by

\[
u=(x+519)(x+39)=x^2+558x+20241,
\]

and `IsDivisibleBy` returns a rational half of exact order 4.  This agrees
with the 2-primary part `(Z/2)^3 x Z/4` of the exact torsion group.

## The 15 immediate Richelot codomains

For each partition of the roots into three pairs, write

\[
F=A B C,
\qquad
[A,B]=A'B-AB',
\]

and let `Delta` be the determinant of the coefficient matrix of `A,B,C`.
The script constructs

\[
g=\Delta[B,C][C,A][A,B].
\]

Every one of the 15 determinants and discriminants is nonzero.  Every
manual codomain is `Q`-isomorphic to exactly one of the 15 surfaces returned
by `RichelotIsogenousSurfaces`; the formula above has the correct twist sign
in all 15 cases.

| pairing | pairs of root indices | exact codomain torsion |
|---:|:---|:---|
| 1 | `12|34|56` | `[2,12]` |
| 2 | `12|35|46` | `[2,12]` |
| 3 | `12|36|45` | `[2,12]` |
| 4 | `13|24|56` | `[2,12]` |
| 5 | `13|25|46` | `[2,12]` |
| 6 | `13|26|45` | `[2,12]` |
| 7 | `14|23|56` | `[2,2,6]` |
| 8 | `14|25|36` | `[2,12]` |
| 9 | `14|26|35` | `[2,12]` |
| 10 | `15|23|46` | `[2,2,6]` |
| 11 | `15|24|36` | `[2,12]` |
| 12 | `15|26|34` | `[2,12]` |
| 13 | `16|23|45` | `[2,2,6]` |
| 14 | `16|24|35` | `[2,12]` |
| 15 | `16|25|34` | `[2,12]` |

Thus the immediate layer consists of twelve exact `[2,12]` Jacobians and
three exact `[2,2,6]` Jacobians, all of order 24.

## Complete rational Richelot component

Depth 2 adds exactly two new `Q`-isomorphism classes, both with exact torsion
`[2,2,6]`.  Expanding those two vertices at depth 3 adds no new vertex, so
the rational Richelot component is closed:

| torsion type | number of vertices |
|:---|---:|
| `[2,2,2,12]` | 1 |
| `[2,12]` | 12 |
| `[2,2,6]` | 5 |
| **total** | **18** |

There are 21 undirected Richelot edges.  The graph has a particularly simple
description:

- the source is joined to twelve degree-1 `[2,12]` leaves;
- the source, the two new depth-2 `[2,2,6]` vertices, and the three immediate
  `[2,2,6]` vertices form a `K_{3,3}`.

The outgoing-kernel counts are visible both from rational 2-torsion and from
the factorization types printed in the log:

- the source has `J[2](Q)` of dimension 4, hence all 15 maximal isotropic
  planes are rational; its sextic has factor type `[1,1,1,1,1,1]`;
- every `[2,12]` vertex has factor type `[2,2,2]`.  The three irreducible
  quadratic orbits give the dual pairing, and Magma finds no other
  Galois-stable pairing, so there is one rational Richelot edge.  Equivalently,
  its pointwise-rational 2-torsion has dimension 2 and is precisely the dual
  maximal isotropic plane;
- every `[2,2,6]` vertex has factor type `[1,1,1,1,2]`.  The conjugate roots
  of the irreducible quadratic must be paired together, while the four
  rational roots can be paired in exactly three ways.  Equivalently, its
  rational 2-torsion has dimension 3; the alternating pairing on this
  3-space has a one-dimensional radical and exactly three isotropic planes.

## Target outcome and simplicity

No vertex has any of

```text
[4,12], [2,4,12], [2,24], [2,2,24],
[2,2,2,24], [2,2,4,12].
```

Therefore the rational Richelot component of this particular record point
does not itself improve the order-96 record.

All 18 Jacobians are nevertheless geometrically simple: every vertex is
`Q`-isogenous to the source, and geometric simplicity is invariant under
isogeny.  The source certificate is

```text
code/verify_record_22212_order96.m
```

which gives independent root-power Frobenius witnesses at `p=37` and `p=73`.
