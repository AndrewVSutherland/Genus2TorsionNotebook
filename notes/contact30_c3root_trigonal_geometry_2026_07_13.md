# The contact-30 C3-root trigonal curve and its genus-6 quotient (2026-07-13)

## Result

The rational-root cover needed by the distinguished Richelot route to
\([2,60]\) is an irreducible degree-three cover of the original
genus-zero parameter line. Its normalization has genus \(12\). It has a
rational involution, and the exact quotient is a trigonal curve of genus
\(6\).

No contact-open rational point was found in either of two height-\(5000\)
searches. Every rational root found is on the previously identified
boundary. This is strong bounded negative evidence, not a determination of
the rational points.

The calculations are reproduced by
[the bounded geometry audit](../code/contact30_c3root_trigonal_bounded_audit.m)
and [the quotient/search companion](../code/contact30_c3root_trigonal_quotient_search.m).
The [compact numerical ledger](../data/contact30_c3root_trigonal_geometry_2026_07_13.txt)
records the commands and counts.

## Original trigonal curve

Let \(F(R,\rho)=0\) be the primitive cubic obtained by requiring the cubic
factor \(C_3\) of the simultaneous contact-\((5,6)\) curve to have root
\(\rho\). Direct calculation gives bidegree \((20,3)\), irreducibility over
\(\mathbf Q(R)\), genus \(12\), and different degree \(28\).

Matching the exact degree/ramification-index profile returned by
\(\operatorname{DifferentDivisor}\) to the factored discriminant identifies
the following branch ledger:

    R=3       e=2
    R=7/3     e=3
    R=2       e=2
    Q4(R)=0   e=2
    Q20(R)=0  e=2

The degrees contribute \(1+2+1+4+20=28\), agreeing with
Riemann--Hurwitz.

## Exact involution and quotient

The curve has the involution

\[
 \sigma(R)=\frac{5R-7}{3R-5},\qquad \sigma(\rho)=\rho.
\]

The companion verifies \(\sigma^2=1\) and equality of the monic cubics
\(F(\sigma(R),\rho)\) and \(F(R,\rho)\). An invariant coordinate is

\[
 v=\frac{3R^2-7}{3R-5}.
\]

Eliminating \(R\) from \(F(R,\rho)\) and
\(3R^2-3vR+5v-7=0\) gives a scalar times \(q(v,\rho)^2\).
The factor \(q\) is irreducible, has bidegree \((10,3)\), and has 44
terms. With Magma's monic normalization, \(3125q\) is primitive integral.

Its exact function field has genus \(6\) and different degree \(16\), with
support profile

    <1,2,1>, <1,3,2>, <1,2,1>, <2,2,1>, <10,2,1>.

Matching this profile to the factored quotient discriminant identifies the
quotient branch ledger:

    v=5        e=2
    v=14/3     e=3
    v=2        e=2
    Q2(v)=0    e=2
    Q10(v)=0   e=2

The displayed candidate base transformation
\(v\mapsto(106v-532)/(21v-106)\) fails to preserve both genuine
nonrational branch factors, so this candidate gives no quotient. This is not
a classification of all automorphisms; it only shows that there is no
obvious further base quotient of this form.

The scripts print the two ledger labels after these profile/discriminant
comparisons. They do not machine-assert the individual support polynomials,
so the matching step remains part of the documented audit.

The double cover back to the \(R\)-curve is controlled by

\[
 \Delta(v)=9v^2-60v+84.
\]

A rational quotient point lifts to rational \(R\) precisely when
\(\Delta(v)\) is a rational square.

## Conservative height searches

For each prime from \(7\) through \(83\), the search retains a
\(\mathbf P^1(\mathbf F_p)\) residue whenever the specialized cubic has a
root or drops degree. Keeping degree drops makes the filter conservative for
roots reducing to \(\rho=\infty\). Survivors of every mask are factored
exactly over \(\mathbf Q\).

For \(R=a/b\), with \(1\le b\le5000\), \(-5000\le a\le5000\), and
\(\gcd(|a|,b)=1\), the result was:

    tested rational parameters  30,401,831
    finite-sieve survivors           73,971
    rational-root hits
      (1,1), (2,-1), (2,0), (3,-1), (3,0), (7/3,1)

The cubic at \(R=\infty\) has no rational root. All six finite hits are
boundary points.

The identical search on the quotient \(v\)-line gave:

    tested rational parameters  30,401,831
    finite-sieve survivors           76,165
    rational-root hits
      (2,1), (5,-1), (5,0), (14/3,1), (32/7,1)
    corresponding Delta
      0, 9, 9, 0, -108/49

The points above \(v=5\) lift to \(R=2,3\), but both are contact boundary.
The zero-discriminant lifts are boundary as well, and the last value has
negative nonsquare discriminant. Thus there is no open lift among these
quotient points.

All three new modes were also run at height \(50\), under a 4 GiB memory
limit and a 120-second timeout. They reproduced the same rational-root lists
and no open lift.

## Consequence for the [2,60] lane

The tested local evidence leaves the rational-root cover locally viable,
while the bounded global searches show it to be thin in the tested boxes.
The genus-\(6\) quotient is materially better than the old
total-degree-\(32\) plane projection, but a larger undirected height search
is not justified. The next defensible step is arithmetic on this exact
curve: its Jacobian, ranks and possible maps to lower-genus curves, followed
by Chabauty or a Mordell--Weil sieve if the rank permits it.

Until such a rational-point argument is completed, the searches establish
only that this route produces no \([2,60]\) candidate in the stated boxes.
