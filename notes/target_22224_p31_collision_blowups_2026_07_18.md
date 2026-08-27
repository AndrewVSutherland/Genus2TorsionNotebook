# The (p=31) collision boundary for ([2,2,2,24])

Date: 2026-07-18

## Result

The six branch-collision strata on the (q=(x+t)^2) slice can be
analyzed without choosing an affine ((M,N)) chart.  The key is that the
full marked order-12 divisibility condition collapses to an explicit local
2-descent square test in the four normalized branch parameters.  This test
agrees exactly with the independent finite-Jacobian computations at
(p=29,31,37,41).

The resulting (31)-adic picture is:

* All three internal collisions
  (A^2=B^2,A^2=C^2,B^2=C^2) are dead at **every smooth depth**.  The only
  12 residue sheets which were not killed by a nonsquare unit have two
  strict-transform constants of opposite squareclass.
* Each external collision
  (A^2=D^2,B^2=D^2,C^2=D^2), (D=sigma/ho), has 12 same-sign residue
  sheets which open into genuine target neighborhoods over
  (mathbf Q_{31}).  It also has 56 opposite-sign sheets which survive
  precisely at even collision depth and satisfy one explicit squareclass
  equation.
* Above (R=S=0), the first exceptional conic has 42 generic live sheets,
  186 killed sheets, and 12 sheets entering an explicitly described deeper
  odd-depth tower.  The two apparently deeper directions
  (a^2+ab+b^2=0) are actually killed because (R-S) has odd order 3.
* A rational search of the three balanced exceptional covers through
  height 300 tested 328,773 presentations.  It found no open point; its ten
  square solutions all lie on the deeper collision towers.
* Exact Magma computations show that both tractable exact collision curves
  have rank zero.  They have no nonzero rational torus points.  No genus-2
  candidate arose, so there was nothing smooth to certify for torsion or
  geometric simplicity.

Thus (p=31) does **not** obstruct the q-square slice globally.  It now gives
a much sharper search prescription: retain only the external same-sign
disks, the external even-depth covers, and the live (R=S=0) exceptional
disks.  All internal collision residues may be discarded rigorously.

## 1. The chart-free divisibility criterion

Put

\[
 (w_1,w_2,w_3,w_4)=(A,B,C,D),\qquad ABC=1,\qquad
 D=\frac{\sigma}{\rho}.
\]

On this slice (q=(x+t)^2).  Therefore the order-3 class

\[
 D_3=[q,h\bmod q]
\]

has trivial image under the local 2-descent map: evaluating (q) at every
Weierstrass point gives a square.  Since (D_{12}=D_3+D_4), divisibility of
(D_{12}) by 2 is equivalent to divisibility of (D_4).

For

\[
 g_4=(x-ab)(x-cd)-x(a+b)(c+d),
\]

exact factorization gives

\[
\begin{aligned}
g_4(0)&=abcd,\\
g_4(-a^2)&=a(a+b)(a+c)(a+d),\\
g_4(-b^2)&=b(a+b)(b+c)(b+d),\\
g_4(-c^2)&=c(a+c)(b+c)(c+d),\\
g_4(-d^2)&=d(a+d)(b+d)(c+d).
\end{aligned}
\]

After removing the common fourth power (s^4), the Kummer coordinates are

\[
 K_0=D,\qquad
 K_i=w_i\prod_{j\ne i}(w_i+w_j),\quad 1\leq i\leq4.
\]

Moreover

\[
 \prod_{i=1}^4K_i
 =D\prod_{i<j}(w_i+w_j)^2.
\]

Consequently the marked order-12 class is divisible by 2 over any field of
odd characteristic exactly when the four (K_i), (1\leq i\leq4), are
squares.  This includes halves on every projective halving chart; it is not
an affine ((M,N)) test.

As a computational check, exhaustive evaluation of this criterion gives:

| (p) | smooth slice sheets | Kummer target sheets | target presentations | previous group-theoretic count |
|---:|---:|---:|---:|---:|
| 29 | 384 | 24 | 72 | 72 |
| 31 | 192 | 0 | 0 | 0 |
| 37 | 768 | 24 | 72 | 72 |
| 41 | 768 | 72 | 216 | 216 |

The match is exact at all four primes.

## 2. Exact collision identities

Write

\[
 x=A^2,\quad y=B^2,\quad z=C^2=(xy)^{-1},
\]

and

\[
 R=x+y+z-3=\rho^2,\qquad
 S=xy+xz+yz-3=\sigma^2.
\]

Two identities drive all blow-ups:

\[
 R-S=(x-1)(y-1)(z-1),
\]

and, for any one of (x,y,z),

\[
 S-xR=-\frac{(x-1)^3}{x}.
\]

Thus an external collision (x=D^2=S/R) is not a reduced generic divisor:

\[
 x-D^2=\frac{(x-1)^3}{xR}.
\]

This cubic thickness is why the external opposite-sign components alternate
between dead and live as their collision depth changes parity.

## 3. All six non-(R=S=0) residue strata

There are 768 normalized branch sheets above the 192 collision bases modulo
31.  For each internal divisor the 32 sheets split as follows:

| linear collision | nonsquare-unit dead | initially zero | unit-square live |
|---|---:|---:|---:|
| same sign (w_i=w_j) | 16 | 0 | 0 |
| opposite sign (w_i=-w_j) | 12 | 4 | 0 |

For each external divisor the 224 sheets split as follows:

| linear collision | nonsquare-unit dead | initially zero | unit-square live |
|---|---:|---:|---:|
| same sign (w_i=D) | 100 | 0 | 12 |
| opposite sign (w_i=-D) | 56 | 56 | 0 |

The raw rows, including all Kummer values and their Legendre symbols, are in
`results/target_22224_p31_collision_blowups_residue.tsv`.

## 4. Internal collisions are completely dead

Consider (x=y).  The other two cases follow by permutation.  The two
q-square equations restrict on the exact collision divisor to

\[
 R=\frac{(x-1)^2(2x+1)}{x^2},\qquad
 S=\frac{(x-1)^2(x+2)}x.
\]

Hence this boundary requires

\[
 r^2=2x+1,\qquad k^2=x(x+2).
\]

Putting (Y=2k) and (x=(r^2-1)/2) gives

\[
 Y^2=r^4+2r^2-3.
\]

Magma sends this quartic to

\[
 E:y^2=x^3-16x^2+96x-192,
\]

with minimal model (y^2=x^3-x^2+x), rank 0, and torsion
(mathbf Z/4mathbf Z).  Its rational quartic points are only
(r=\pm1,Y=0) and the two points at infinity.  The affine points give
(x=0), outside the torus; the other two are boundary points.  Thus the
exact internal collision curve has no nontrivial rational q-square point.

The strict transform is stronger.  Suppose (w_i\equiv-w_j\pmod {31}),
and let

\[
 w_i^2-w_j^2=31^m n,qquad n\in\mathbf Z_{31}^{\times}.
\]

The common small sum satisfies

\[
 \frac{w_i+w_j}{31^m}=\frac{n}{w_i-w_j}.
\]

For every one of the 12 initially possible internal sheets, the two Kummer
coordinates containing this sum have normalized constants of Legendre
symbols

\[
 (+1,-1).
\]

If (m) is odd, their valuations are odd.  If (m) is even, multiplication
by the same normal unit (n) cannot make constants of opposite squareclass
both square.  Every smooth lift is therefore killed, at every depth.

In normalized lift counts this resolves all

\[
 96\cdot31^2=92,256
\]

internal lifts modulo (31^2) as dead; there is no residual (31^3) case.

## 5. External collisions: the live components

It is enough to treat (A^2=D^2); the other two are permutations.

### 5.1 Same-sign component

The only live orientation has

\[
 A=D=1,\qquad C=B^{-1}.
\]

On the boundary,

\[
\begin{aligned}
K_A=K_D&=\frac{2(B+1)^2}{B},\\
K_B&=(B+1)^2(B^2+1),\\
K_C&=\frac{(B+1)^2(B^2+1)}{B^4}.
\end{aligned}
\]

Thus its exact square cover is

\[
 e^2=\frac2B,\qquad f^2=B^2+1.
\]

Over (mathbf F_{31}), where 2 is a square, the six allowed (B)-residues
are

\[
 2,7,9,10,16,28.
\]

The two signs of (ho) give 12 sheets for each external divisor.  Every
one opens into a full (mathbf Q_{31})-target neighborhood because all
Kummer coordinates are units and their squareclasses cannot change.

Over (mathbf Q), put (B=2t^2).  The cover becomes

\[
 f^2=4t^4+1.
\]

Its elliptic model has minimal equation (y^2=x^3-x), rank 0, and torsion
((\mathbf Z/2\mathbf Z)^2).  The only affine rational points have (t=0),
so (B=0); the remaining points are at infinity.  Hence the exact singular
component contains no rational torus point, although its listed 31-adic
disks are live and must be retained in a global search.

### 5.2 Opposite-sign component

The only potentially live orientation has

\[
 A=-1,\qquad D=1,\qquad C=-B^{-1}.
\]

The two nonzero Kummer coordinates are already squares:

\[
 K_B=(B^2-1)^2,qquad
 K_C=\frac{(B^2-1)^2}{B^4}.
\]

Let

\[
 m=v_{31}(A^2-1),\qquad A^2-1=31^m u,quad
 u\in\mathbf Z_{31}^{\times}.
\]

Since (R) is a unit on this nonintersection stratum, the cubic identity
gives

\[
 v_{31}(A+D)=3m.
\]

After division by (31^{3m}), both vanishing Kummer coordinates have the
same squareclass

\[
 c(B)u^3,qquad c(B)=-\frac{B}{2(B^2-1)}.
\]

Therefore:

* odd (m) is killed by odd valuation;
* even (m) is live exactly on

  \[
  W^2=c(B)u.
  \]

Every even-depth exceptional component is rational; one parameterization is

\[
 u=-\frac{2(B^2-1)}B,W^2.
\]

Among the 168 opposite-sign residue sheets over all three external
divisors, the constant (c(B)) is a square on 84 and a nonsquare on 84.

### 5.3 Exact lift counts

The strict-transform count through (31^2) is

| external strata | resolved live | killed | deep mod (31^2) | total |
|---|---:|---:|---:|---:|
| all three divisors | 34,596 | 605,988 | 5,208 | 645,792 |

The 5,208 deep classes are exactly those with (m\ge2).  Lifting only
these to (31^3) gives

| status at (31^3) | count |
|---|---:|
| resolved live at (m=2) | 2,421,720 |
| killed at (m=2) | 2,421,720 |
| still deep, (m\ge3) | 161,448 |
| total | 5,004,888 |

The remaining depth is not conceptually unresolved: the odd/even rule and
(W^2=c(B)u) classify every finite smooth depth.

Combining all six collision divisors gives, modulo (31^2),

\[
 738,048=34,596\text{ live}+698,244\text{ killed}+5,208\text{ deep}.
\]

## 6. The intersection (R=S=0)

The four residue bases have (A^2=B^2=C^2=1).  Write

\[
 A=\epsilon_A(1+31a+\cdots),\qquad
 B=\epsilon_B(1+31b+\cdots).
\]

Then

\[
 R=4\cdot31^2(a^2+ab+b^2)+O(31^3),
\]

and the same formula holds for (S).  Put

\[
 q^2=a^2+ab+b^2.
\]

When (q\ne0), the first normalized branch deviations are

\[
 \ell=(\ell_A,\ell_B,\ell_C,\ell_D)
 =\left(a,b,-a-b,\frac{ab(a+b)}{q^2}\right).
\]

The exceptional conic is rational.  A convenient parameterization is

\[
 (a:b:q)=(2m-1:1-m^2:m^2-m+1).
\]

The sign constraint (ABC=1) permits zero or two negative signs among
(A,B,C).

* With all four normalized branches positive, every Kummer coordinate is
  (8) modulo 31, hence square.  All 15 q-square directions are live.
* With one or three negative signs, the four Kummer valuations are odd;
  these sheets are killed.
* With two negative signs and (D>0), the leading exceptional cover is

  \[
  Y_i^2=2\prod_{\epsilon_j=-\epsilon_i}(\ell_i-\ell_j),
  \qquad i=1,2,3,4.
  \]

  There are three sign patterns:
  ((+--+),(-+-+),(--++)).  For each pattern, among the 15 q-square
  directions, 2 are live, 11 are killed, and 2 meet a deeper external
  collision.

Including the two signs of the first square root gives 42 generic live
sheets, 186 killed sheets, and 12 deeper intersection-tower sheets.

For one of the 12 tower sheets, let (x_i-1=31^m u) with (m\ge2).
Here (v_{31}(R)=2), and the exact cubic identity gives

\[
 v_{31}(w_i+D)=3m-2,qquad v_{31}(K_i)=3m-1.
\]

Thus even (m) is killed.  At odd (m\ge3), the next strict-transform
equation is

\[
 W^2=c u,qquad c=\frac{\ell_i-\ell_k}{4},
\]

where (k) is the other opposite-sign vertex.  Three of the six projective
tower rows have square (c), and three have nonsquare (c); the two square
root sheets double these to 12.  This gives the requested next equation for
every deeper intersection component.

Finally, the two directions (q=0) are not live.  Indeed

\[
 R-S=-8\cdot31^3ab(a+b)+O(31^4),
\]

and the displayed coefficient is nonzero at both roots of
(a^2+ab+b^2=0) in (mathbf F_{31}).  Therefore (R) and (S) cannot
both jump to even valuation at least 4: one has valuation 3 and is not a
square.

## 7. Rational searches and exact verification

The balanced exceptional equations were searched after substituting the
conic parameterization above.  For every reduced

\[
 m=n/d,\qquad |n|,d\le300,
\]

and all three balanced sign patterns, the program tested whether all four
leading Kummer coordinates were rational squares.  It tested 328,773
presentations.

There were no open hits.  The ten square solutions occur only at

\[
 m\in\{-1,0,1,2,1/2\},
\]

and each has at least two zero Kummer coordinates; these are precisely the
deeper collision-tower points described above.

The exact Magma script independently verifies:

* the five factorizations defining the Kummer coordinates;
* the record point has four nonsquare Kummer coordinates (and (D) itself
  is nonsquare);
* the internal quartic has rank 0 and torsion (mathbf Z/4mathbf Z);
* the external quartic has rank 0 and torsion
  ((\mathbf Z/2mathbf Z)^2);
* all rational points on both quartics are the invalid zero/pole boundary
  points listed above.

No smooth rational curve was produced, so no alleged
([2,2,2,24]) hit remains unverified.

## 8. Files

* `code/target_22224_p31_collision_blowups.py`: finite controls, all residue
  sheets, strict-transform counts, (R=S=0) exceptional divisor, and the
  rational height-300 search.
* `code/target_22224_p31_collision_blowups.m`: exact Kummer factorization and
  rank-zero collision-curve certifications.
* `results/target_22224_p31_collision_blowups.log`: primary numerical log.
* `results/target_22224_p31_collision_blowups_exact.log`: Magma exact log.
* `results/target_22224_p31_collision_blowups_residue.tsv`: all 768
  nonintersection collision sheets.
* `results/target_22224_p31_collision_blowups_rs_intersection.tsv`: every
  projective direction and sign stratum above (R=S=0).
* `results/target_22224_p31_collision_blowups_balanced_rational.tsv`:
  rational exceptional-cover search output.

## 9. Recommended next computation

Replace the old undifferentiated “(p=31) boundary” flag in the rational
sieve by the following exact masks:

1. discard all internal collision residues;
2. on an external same-sign residue, retain only
   (A=D=1) (or its permutation) with
   (B\in\{2,7,9,10,16,28\});
3. on an external opposite-sign residue, enforce even collision depth and
   (W^2=c(B)u);
4. at (R=S=0), retain the 42 generic live sheets and the 12 odd-depth
   tower sheets with their displayed (W^2=cu) equation.

These conditions are cheap integer valuation and Legendre-symbol tests.
They should be combined with the positive masks at 29, 37, and 41 before
the next rational or lattice search.
