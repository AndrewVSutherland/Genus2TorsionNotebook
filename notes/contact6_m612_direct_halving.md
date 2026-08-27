# Direct contact-6 attack on [6,12]

## Universal and compact halving covers

The contact-6 polynomial factors as

\[
f=xQ_1Q_2,\quad
Q_1=(b+3)x^2+(a-3)x+2,\quad
Q_2=2x^2+(b-3)x+(a+3).
\]

Every smooth fiber has rational 2-rank at least two. For the marked
order-six contact class \(D\), exact arithmetic gives \(3D=[Q_1,0]\).
A rational half of a nonzero 2-class raises the 2-primary subgroup from
[2,2] to [2,4]; with rational 3-primary [3,3], this gives [6,12].

All six universal square-quartic covariants are irreducible over
\(\mathbf Q\), and each pair has gcd one. For \(T_B=[Q_1,0]\), quadratic
recovery in the auxiliary \(N\) gives

\[
M^8(b+3)^4R_{26}(a,b,M),
\]

where \(R_{26}\) is irreducible with multidegree \((10,9,16)\) and 276
terms. On all four recorded exact [6,6] curves, its only rational linear
factor is the boundary \(M=0\). See
code/contact6_m612_halving_equations.m and
code/contact6_m612_tb_resultant.m.

For a compact chart, compare the leading and constant coefficients of

\[
S=Q_1(Mx+N)^2-xQ_2.
\]

On a smooth affine half, this forces \(b=2s^2-3\). Put
\(M=1/m,\ N=r/m,\ K=m^2\), and \(A_3=a-3+4s^2r-2K\). The complete square
condition becomes

\[
\begin{aligned}
H_1={}&s((a-3)r^2+4r-K(a+3))-rA_3,\\
H_2={}&8s^2(2s^2r^2+2(a-3)r+2-K(2s^2-6))
       -A_3^2-32s^3r.
\end{aligned}
\]

\(H_1\) is linear in \(a\), with denominator
\(a_{\rm den}=s(r^2-K)-r\). After recovering \(a\) and removing \(s^2K\),
\(H_2\) leaves one irreducible equation \(H(s,r,K)=0\), of total degree
eight, degree three in \(K\), and 31 terms. The exact half is

\[
M_{\rm tool}=\frac{2s^2}{m},\quad
N_{\rm tool}=\frac{2s^2r}{m},\quad
G=x^2+\frac{A_3}{4s^2}x+\frac rs.
\]

The generic search checked 1,232,100 rational \((s,r)\) pairs at height 30.
It found 146 square-\(K\) lifts and 81 smooth curves; every smooth curve was
killed by an exact finite [6,12] test.

## Genuine recovery-boundary family

The locus \(a_{\rm den}=a_{\rm num}=0\), with \(sr\ne0\), splits as

\[
\begin{array}{c|c}
r&K\\ \hline
1/s&0,\\
s(2s+3)&(s+1)^2(2s-1)(2s+3).
\end{array}
\]

The first component has \(a=4s+3\) and \(Q_1=2(sx+1)^2\), so it is
singular. The second is genuine. Its discriminant in \(a\) is

\[
1024s^2(s+1)^5(s-\tfrac12).
\]

Requiring both \(K\) and this discriminant to be squares gives a genus-zero
fiber product. Put

\[
\begin{aligned}
D_u&=(u^2-1)(u^2-9),\\
s&=-\frac{3u^4-14u^2+27}{2D_u},\\
y&=\frac{u^4-9}{D_u},\\
w&=-\frac{8u(u^2-3)}{D_u}.
\end{aligned}
\]

Then \(y^2=(s+1)(s-\tfrac12)\), \(w^2=(2s-1)(2s+3)\), and

\[
r=s(2s+3),\quad m=(s+1)w,\quad K=m^2,\quad b=2s^2-3.
\]

The two exact branches are

\[
a_-=
\frac{9u^{12}+190u^{10}-565u^8+2484u^6-10449u^4+21870u^2-19683}
{u^{12}-38u^{10}+559u^8-3924u^6+12879u^4-16038u^2+6561},
\]

and

\[
a_+=
\frac{-3u^{12}+30u^{10}-129u^8+276u^6-565u^4+1710u^2+729}
{u^{12}-22u^{10}+159u^8-436u^6+559u^4-342u^2+81}.
\]

Both depend only on \(u^2\). Specializing the independent cubic-contact
equations to either branch gives a cover in \((u,M,U,v)\):

| equation | total degree | \(u\)-degree | terms |
|---|---:|---:|---:|
| \(F_1\) | 21 | 16 | 62 |
| \(F_2\) | 36 | 32 | 289 |
| \(F_3\) | 27 | 24 | 130 |

Their common gcd is one.

This family is geometrically simple. At \(u=2\), \(s=19/30\) and
\(b=-989/450\). The two fibers are:

* \(a=1601/2025\): torsion [2,12], D4 certificate at \(p=29\), with
  \(\Phi=T^4+6T^3+34T^2+174T+841\).
* \(a=48791/1875\): torsion [2,12], D4 certificate at \(p=37\), with
  \(\Phi=T^4-6T^3+10T^2-222T+1369\).

Both have automorphism group of order two.

## Exact masks and height-100000 exclusion

The script code/contact6_m612_tb_recovery_family.m computes independent
finite extra-3 masks for both branches. Every good finite fiber is checked
to contain the base [2,12]. It is retained only if it contains [6,12];
every bad residue and every denominator residue passes conservatively.

Masks were computed through \(p=251\). The C sieve
code/contact6_m612_tb_recovery_height_sieve.c uses CRT primes 17, 79, 83,
then every mask through 251. Since the family is even in \(u\), it scans
primitive

\[
0\le {\rm numerator}(u)\le100000,\qquad
1\le {\rm denominator}(u)\le100000.
\]

For each branch it processed about 178 million CRT tuples and 96 million
primitive candidates. The only survivors were \(u=0,1,3\), all parameter
boundaries. Hence neither branch has a smooth [6,12] candidate through
projective height 100000.

Principal logs are:

* data/contact6_m612_tb_recovery_family_summary.txt
* data/contact6_m612_tb_recovery_masks_p251.txt
* data/contact6_m612_tb_recovery_sieve_h100000_p251.txt
* data/contact6_m612_tb_recovery_sample_verify.txt
