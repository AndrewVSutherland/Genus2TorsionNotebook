# The symmetric \((2,2,2,12)\) curve

Date: 2026-07-20.

## Result

Let

\[
\mathcal S:\qquad
x^2+y^2+z^2=u^2+v^2,\qquad
x^4+y^4+z^4=u^4+v^4
\]

in \(\mathbf P^4\).  A rational point
\(P=[x:y:z:u:v]\) on the natural smooth open of this surface gives the
genus-two curve

\[
\boxed{\quad
C_P:\ W^2=T(T-x^2)(T-y^2)(T-z^2)(T-u^2)(T-v^2).
\quad}
\]

Its Jacobian contains

\[
(\mathbf Z/2\mathbf Z)^3\times\mathbf Z/12\mathbf Z.
\]

Thus the answer in the symmetric parameters is considerably simpler than
the original \((s,m,n)\)-equation: the five nonzero branch points are just
the five squared projective coordinates.

This statement concerns the marked subgroup contained in the Jacobian.  The
generic torsion group need not be computed in order to obtain the family.
For each of the three known rational points, Magma computes the full torsion
group to be exactly \([2,2,2,12]\).

## 1. The two surface identities

Put

\[
p=x^2,\quad q=y^2,\quad r=z^2,\quad b=u^2,\quad c=v^2.
\]

The power-sum equations defining \(\mathcal S\) give

\[
p+q+r=b+c,\qquad pq+pr+qr=bc.
\]

Consequently

\[
(T-p)(T-q)(T-r)=T(T-b)(T-c)-pqr.
\]

This single constant-difference identity explains both the simple curve
equation and its torsion.

## 2. Direct proof of the torsion subgroup

The six branch points \(0,p,q,r,b,c\) are rational, so the Jacobian has full
rational two-torsion:

\[
J[2](\mathbf Q)\simeq(\mathbf Z/2\mathbf Z)^4.
\]

For the order-three point, set

\[
G(T)=(T-p)(T-q)(T-r),\qquad H(T)=T(T-b)(T-c).
\]

Then \(G=H-pqr\) and \(W^2=GH\).  If
\(\phi=2W-G-H\), then

\[
\phi\,\iota(\phi)=(G-H)^2=(pqr)^2,
\]

where \(\iota(W)=-W\).  Comparison at the two rational points at infinity
gives, after labeling them,

\[
\operatorname{div}(\phi)=3\infty_+-3\infty_-.
\]

The class \(\infty_+-\infty_-\) therefore has exact order three.

For the order-four point, let \(W_b,W_c\) be the Weierstrass points over
\(b,c\).  The change

\[
R=\frac{T-b}{c-T}
\]

sends \(W_b\) to zero and \(W_c\) to infinity.  A remaining branch point
\(e\in\{0,p,q,r\}\) maps to

\[
a_e=\frac{e-b}{c-e},\qquad
-a_e=\frac{(b-e)(c-e)}{(c-e)^2}.
\]

The four numerators are squares:

\[
\begin{array}{c|c}
e & (b-e)(c-e)\\ \hline
0 & bc=(uv)^2\\
p & qr=(yz)^2\\
q & pr=(xz)^2\\
r & pq=(xy)^2.
\end{array}
\]

The scalar in the transformed odd-degree model is also the square
\((c-b)^2(vxyz)^2\).  The standard two-descent halving criterion now shows
that \(W_b-W_c\) has a rational half, necessarily of exact order four.
Full rational two-torsion, this order-four class, and the independent
order-three class generate

\[
(\mathbf Z/2\mathbf Z)^3\times\mathbf Z/12\mathbf Z.
\]

## 3. Smooth locus

For

\[
f(T)=T(T-x^2)(T-y^2)(T-z^2)(T-u^2)(T-v^2),
\]

one has

\[
\operatorname{Disc}(f)=
(xyzuv)^4
\prod_{\alpha<\beta}
(\alpha^2-\beta^2)^2,
\qquad
\alpha,\beta\in\{x,y,z,u,v\}.
\]

Hence \(C_P\) is smooth precisely when all five coordinates are nonzero and
their squares are pairwise distinct.

## 4. Birational map from the old square cover

Choose the chart

\[
(A,B,C)=(x,u,y)
\]

and put

\[
d=x^2+y^2-u^2=v^2-z^2,\qquad
e=x^2+2y^2-u^2.
\]

The two old sextics factor on \(\mathcal S\) as

\[
F=d^2z^2,\qquad G=d^2v^2.
\]

Thus their chosen square roots are \(Y_0=dz\) and \(Z_0=dv\).  The old
parameters become

\[
\boxed{\qquad
(s:m:n)=
\bigl(x^2e:2(y^2-u^2)e:2x^2d\bigr).
\qquad}
\]

Equivalently, in the affine normalization from the original formulas,

\[
s=\frac{x e}{2ud},\qquad
m=\frac{(y^2-u^2)e}{xud},\qquad
n=\frac{x}{u}.
\]

Here is an explicit check that the two curve equations agree without a
quadratic twist.  Put \(\Lambda=2x^2e\).  After substituting the displayed
projective \((s,m,n)\), scaling the old horizontal coordinate by
\(\Lambda\), and removing a square from the ordinate, its equation is

\[
\begin{aligned}
Y_1^2={}&R(R+pq)(R+(b-p)d)(R+bd)\\
&\mathrel{}\cdot(R+(b-q)d)(R+pq-d^2).
\end{aligned}
\]

On \(\mathcal S\),

\[
pq=(b-r)d,\qquad pq-d^2=(b-c)d.
\]

The substitution \(R=d(T-b)\), followed by division of the ordinate by
\(d^3\), gives exactly

\[
W^2=T(T-p)(T-q)(T-r)(T-b)(T-c).
\]

The other eleven old charts are obtained by permuting \(x,y,z\) and
interchanging \(u,v\).  The symmetric equation makes it immediate that all
twelve reconstructed curves are \(\mathbf Q\)-isomorphic.

Conversely, from a point \((A,B,C,Y_0,Z_0)\) on the old square cover with
\(D=A^2-B^2+C^2\ne0\), one obtains

\[
[x:y:z:u:v]=[AD:CD:Y_0:BD:Z_0]
\]

in homogeneous coordinates.

The extra old factor \(e\) is only a coordinate-chart artifact.  In fact,
\(e=0\) has no nonzero rational point on this surface: after clearing
denominators primitively it forces

\[
u^2=x^2+2y^2,\qquad z^2=2(x^2+y^2),
\]

so \(x,y\) must be odd and then \(u^2\equiv3\pmod 8\), a contradiction.
Thus the old and natural nondegenerate opens have the same rational points.

## 5. Known rational points

\[
\begin{array}{c|c}
(x,y,z,u,v)&(s,m,n)\\ \hline
(408,143,1015,437,1013)&(336396,-689185,-166464)\\
(120,143,266,218,241)&(2208,-8303,-7200)\\
(16660,78793,21456,78644,27593)&
(1500518600,253638081,138777800)
\end{array}
\]

The canonical unordered keys are

\[
[143,408,1015]\mid[437,1013],
\]

\[
[120,143,266]\mid[218,241],
\]

and

\[
[16660,21456,78793]\mid[27593,78644].
\]

The exact symbolic identities and all three torsion computations are checked
by \(\texttt{code/verify\_22212\_symmetric\_model.m}\).

