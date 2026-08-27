# Doubled-theta loci on the Clebsch--Klein surface

Date: 2026-07-14

## Objective

Test the proposed shortcut

\[
2P-K\in J[10]
\]

on the Clebsch--Klein family with rational torsion
\([2,2,2,10]\).  In the odd model

\[
 Z^2=g(X)=\prod_{i=1}^5(1-r_i^2X),\qquad
 T=(0,1)-\infty,\qquad \operatorname{ord}(T)=5,
\]

the CK equations are

\[
\sum r_i=\sum r_i^3=0.
\]

Up to sign and the \(S_5\)-symmetry, the nontrivial order-10
classes that can be of doubled-theta form reduce to

\[
A=T+S_{12},\qquad B=2T+S_{01},\qquad C=2T+S_{12}.
\]

The executable certificate is
[code/agent_ck_doubled_theta_geometry.m](../code/agent_ck_doubled_theta_geometry.m);
its saved output is
[data/agent_ck_doubled_theta_geometry_2026_07_14.txt](../data/agent_ck_doubled_theta_geometry_2026_07_14.txt).

## Direct CK-chart calculation

Interpolation constructs the Mumford quadratic for each class and
factors its discriminant on the rational two-parameter CK chart.
In every case there is exactly one odd irreducible factor.  The three
odd cores have the following chart data:

| class | bidegree | total degree | normalization genus | projective points of height at most 1000 |
|---|---:|---:|---:|---|
| \(A\) | \((8,8)\) | 16 | 13 | \((0:1:0),(1:0:0),(-1:1:1)\) |
| \(B\) | \((8,8)\) | 16 | 7 | \((0:1:0),(-1:0:1),(0:1:1),(1:0:0),(-1:1:1)\) |
| \(C\) | \((8,8)\) | 16 | 4 | the same five points as \(B\) |

The height search is only diagnostic.  The conclusions below use
exact quotient or positivity arguments.

## Symmetric discriminant identities

For \(A\) and \(C\), put

\[
s=r_1+r_2,\quad p=r_1r_2,\quad
q=r_3r_4+r_3r_5+r_4r_5,\quad u=r_3r_4r_5.
\]

The CK equations give \(u=s(p-q)\).  Up to nonzero square factors,
the relevant Mumford quadratics are

\[
U_A=(q-p)\left(2+(p+q-2s^2)Y-s^2(q-p)Y^2\right),
\]

\[
U_C=p(q-p)\left(2+2qY+p(q-p)Y^2\right).
\]

Their discriminants are

\[
\operatorname{disc}(U_A)
=(q-p)^2\left((p+q-2s^2)^2+8s^2(q-p)\right),
\]

\[
\operatorname{disc}(U_C)
=4p^2(q-p)^2\left(q^2-2p(q-p)\right),
\]

and the odd core in the latter identity is

\[
q^2-2p(q-p)=(q-p)^2+p^2.
\]

For \(B\), single out \(r_1\), and write the elementary symmetric
functions of the other four roots as

\[
e_1=-r_1,\quad e_2=q,\quad e_3=-r_1q,\quad e_4=w.
\]

After normalizing \(r_1=1\), the quadratic and its discriminant are

\[
U_B=w(wZ^2+2qZ+2),\qquad
\operatorname{disc}(U_B)=4w^2(q^2-2w).
\]

## Case A: exact genus-2 Chabauty quotient

On the smooth CK open, \(s\ne0\) and \(q-p\ne0\), so normalize
\(s=1\).  The doubled-theta core is parameterized by

\[
p=(v+1)^2,\qquad q=1+2v-v^2,\qquad u=2v^2.
\]

Rational splitting of the \(r_1,r_2\) pair gives

\[
v=-\frac{k^2-k+1}{k^2+1},\qquad
r_1=\frac1{k^2+1},\qquad r_2=\frac{k^2}{k^2+1}.
\]

The complementary cubic splits over \(\mathbf Q\) only if

\[
\begin{aligned}
D_{12}(k)={}&8k^{12}-80k^{11}+256k^{10}-640k^9+1097k^8
-1568k^7\\
&+1702k^6-1568k^5+1097k^4-640k^3
+256k^2-80k+8
\end{aligned}
\]

is a square.  This polynomial is irreducible and defines a genus-5
hyperelliptic curve.  In fact, for

\[
v=-\frac{k^2-k+1}{k^2+1},\quad
q=1+2v-v^2,\quad u=2v^2,
\]

the complementary cubic is \(z^3+z^2+qz-u\), and the script verifies
the exact identity

\[
\operatorname{disc}(z^3+z^2+qz-u)
=\frac{D_{12}(k)}{(k^2+1)^6}.
\]

Reciprocity gives the genus-2 quotient
\(h=k+k^{-1}\):

\[
y^2=8h^6-80h^5+208h^4-240h^3+145h^2-48h+4.
\]

Magma returns Jacobian rank bounds \([0,1]\), trivial torsion, and the
divisor

\[
D_A=[h^2,\,2-12h]
\]

has infinite order.  Hence the rank is exactly one.  Chabauty with
this generator proves that the only rational quotient points are
\((h,y)=(0,\pm2)\).  But \(h=0\) would require \(k^2+1=0\), so neither
point lifts to rational CK data.

## Case B: complete rank-zero elliptic enumeration

The complementary four roots must split in

\[
z^4+z^3+qz^2+qz+\frac{q^2}{2}.
\]

Its discriminant is

\[
\frac{q^3}{4}\left(32q^3+56q^2-35q-16\right).
\]

Thus complete splitting forces a rational point on the elliptic
curve

\[
y^2=q(32q^3+56q^2-35q-16).
\]

With \((q,y)=(0,0)\) as origin, a minimal model is

\[
E:\quad y^2+xy+y=x^3-82x-92.
\]

The exact rank bounds are \([0,0]\) and
\(E(\mathbf Q)_{\rm tors}\simeq\mathbf Z/3\mathbf Z\).
All three torsion points were mapped back to the quartic model; this
gives exactly

\[
(q,y)=(0,0),\qquad
\left(-\frac4{17},\pm\frac{316}{289}\right).
\]

The first is CK boundary.  At \(q=-4/17\), the required splitting
polynomial, after clearing denominators, is

\[
289z^4+289z^3-68z^2-68z+8,
\]

which is irreducible over \(\mathbf Q\).  Hence no open specialization
has four rational complementary roots.

## Case C: real obstruction

The omitted factors \(p=0\) and \(q=p\) are CK boundary because they
force respectively \(r_1r_2=0\) and
\(r_3r_4r_5=s(p-q)=0\).  On the smooth open the residual equation is

\[
(q-p)^2+p^2=0.
\]

It has only \(p=q=0\) over \(\mathbf R\), again on the boundary.

## Conclusion

None of the three symmetry types has a smooth rational CK point.
Consequently, the most direct doubled-theta attempt does **not**
upgrade this Clebsch--Klein \([2,2,2,10]\) family to
\([2,2,2,20]\).

This is a useful negative result rather than a rejection of the
\(2P-K\) strategy in general: it exhausts this specific CK family
because its symmetry reduces the order-10 translates to the three
cases above.

## Reproduction

The successful command was

    ulimit -v 4194304
    timeout 300s magma code/agent_ck_doubled_theta_geometry.m

Magma V2.29-4 completed in 2.200 seconds with peak reported usage
129.59 MB.  The script also sets its internal memory limit to 4 GB.
