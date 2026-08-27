# The \(A,B,C\) surface for \((2,2,2,12)\): geometry and exhaustive search

Date: 2026-07-20.

## Headline

The search requested from the two sextic square conditions produced a **third
geometrically simple genus-two Jacobian over \(\mathbf Q\) with torsion exactly
\([2,2,2,12]\)**. It is not isomorphic to either example found on 2026-07-18.

The smallest-height chart in the direct search is

\[
(A,B,C)=(4165,19661,5364),
\]

with

\[
y=6705977678400,\qquad z=2348407118400.
\]

It first appears between direct \(ABC\)-height \(10000\) and \(20000\). A
complete search through height \(20000\) found exactly three surface orbits:
the two old ones and this one.

Exact curve data and certificates are in data/22212_abc_curve3.txt.
Reproduction and verification use code/search_22212_abc.cpp and
code/verify_22212_abc_hit.m.

## 1. Factorization of the two sextics

Put

\[
D=A^2-B^2+C^2,\qquad
L=(A^2-B^2)(B^2-C^2).
\]

The two apparently unrelated sextics factor as

\[
F=DL,\qquad G=D(L+D^2).
\]

In particular,

\[
G-F=D^3.
\]

For an exact final test, put

\[
g=\gcd(|D|,|L|),\qquad d=D/g,\quad \ell=L/g,\quad q=(L+D^2)/g.
\]

Since \(\gcd(d,\ell)=\gcd(d,q)=1\), both \(F\) and \(G\) are squares if and
only if \(d,\ell,q\) have the same sign and their absolute values are integer
squares. This avoids square-testing general degree-six integers.

The two old primitive points are recovered as

\[
(408,437,143),\qquad (120,218,143),
\]

with square roots

\[
(4116840,4108728),\qquad (3371550,3054675),
\]

respectively.

## 2. Birational projective model and Kodaira dimension

On \(D\ne0\), set

\[
t=y/D,\qquad w=z/D.
\]

The main open of the fiber product is birational to the complete intersection
\(X\subset\mathbf P^4_{[a:b:c:t:w]}\)

\[
b^2+w^2=a^2+c^2+t^2,
\]

\[
b^2w^2=a^2c^2+a^2t^2+c^2t^2,
\]

where \(a=A,b=B,c=C\). The inverse on this open is \(y=Dt,z=Dw\).

The same complete intersection has the more symmetric presentation

\[
x^2+y^2+z^2=u^2+v^2,\qquad
x^4+y^4+z^4=u^4+v^4,
\]

under \([x:y:z:u:v]=[a:c:t:b:w]\).  Most importantly, its associated
genus-two curve simplifies to

\[
W^2=T(T-x^2)(T-y^2)(T-z^2)(T-u^2)(T-v^2).
\]

A direct torsion proof and the explicit map back to \((A,B,C)\) and
\((s,m,n)\) are recorded in
\(\texttt{notes/22212\_symmetric\_model.md}\).

An exact Jacobian-ideal calculation gives 36 isolated rational singular
points on \(X\), all ordinary double points: 12 lie in the evident symmetry
orbit of \([1:1:0:0:0]\), and 24 in the orbit of
\([1:1:0:1:1]\). Thus \(X\) is a normal nodal \((2,4)\) complete
intersection. Adjunction and crepant resolution of the nodes give

\[
\kappa=2,\qquad K^2=8,\qquad p_g=5,\qquad q=0,\qquad
\chi(\mathcal O)=6,\qquad c_2=64.
\]

So the birational main surface is of **general type**. This does not prove
that the nondegenerate rational points are finite, but it gives no reason to
expect a rational parametrization or a Zariski-dense rational family.

The evident effective projective symmetry subgroup is

\[
\bigl((\mathbf Z/2)^5/\langle-1\rangle\bigr)\rtimes(S_3\times S_2),
\]

of order 192. Here \(S_3\) permutes \(\{a,c,t\}\), \(S_2\) permutes
\(\{b,w\}\), and signs are independent modulo overall projective sign.
This is an evident subgroup, not a claim about the full automorphism group.
The supplied discriminant open is not invariant in a fixed \(A,B,C\) chart,
so all 12 chart assignments are retested.

## 3. Quadratic-height divisor search

Let

\[
U=a^2+c^2,\qquad R=w^2-U,\qquad S=b^2-U.
\]

The complete-intersection equations are exactly equivalent to

\[
RS=(ac)^2
\]

and

\[
w^2=U+R,\qquad b^2=U+S,\qquad t^2=U+R+S.
\]

After using signs, \(S_3\), and \(S_2\), take
\(0<a\le c\le t\) and \(0<b\le w\). Enumerating the signed divisor pairs
\(RS=(ac)^2\) is complete for bounded canonical projective height, with
streaming memory and average \(T^2\operatorname{polylog}T\) work.

The complete canonical-height search through \(T=10000\) examined
14,828,863,891 unordered divisor pairs. It found only

\[
[120,143,266]\mid[218,241]
\]

and

\[
[143,408,1015]\mid[437,1013].
\]

There was no third nondegenerate orbit through canonical height 10000. This
run used about 4 MB of memory.

Canonical height and \(ABC\)-height are different. From integral
\((A,B,C,y,z)\), an integral representative of the new model is

\[
[AD:BD:CD:y:z],
\]

followed by division by the common gcd. Consequently, a canonical-height
search through \(T\) is not an exhaustive \(ABC\)-height-\(T\) search.

## 4. Direct exhaustive \(ABC\)-height search

Because both sextics depend only on squares, take \(A,B,C>0\), primitive.
The equations are invariant under \(A\leftrightarrow C\), so the program
scans \(A<C\) and checks the full discriminant separately in both
orientations.

With \(A<C\), positivity is possible only when \(B>A\) and \(C\) lies in one
of the two exact ranges

\[
B<C\le H
\]

or

\[
A<C<\sqrt{B^2-A^2}.
\]

The second range is present only when \(B^2>2A^2\). These conditions leave
asymptotically

\[
\frac{4+\pi}{24}H^3
\]

positive triples.

The implementation processes \(C\) in blocks of 64. Precomputed local
double-square masks for \(p=79,73,71,59,43,31\) reject almost every block;
the few set bits are tested modulo \(109,97,89,83,64\), then by the exact
gcd criterion above. The sieve tables occupy about 16 MB independent of the
number of threads.

Complete results:

| \(ABC\) height | valid chart representatives | surface orbits |
|---:|---:|---:|
| 500 | 14 | 2 old |
| 2,000 | 24 | 2 old |
| 5,000 | 24 | 2 old |
| 10,000 | 24 | 2 old |
| 20,000 | 26 | 2 old + 1 new |

The height-20000 run statistics were:

- 2,380,237,958,933 positivity-eligible triples;
- 37,359,400,237 blocks of 64;
- 7,170,101 candidates reaching the exact test;
- 7,775 primitive double-square points before the discriminant;
- 26 valid oriented charts;
- 84.219 seconds on 16 threads;
- 16 MB peak resident memory.

Thus the result is exhaustive for primitive rational projective
\((A:B:C)\) with \(\max(|A|,|B|,|C|)\le20000\), modulo the harmless sign
normalization, and subject to the full nonvanishing condition in the question.

## 5. The third orbit and curve

The direct-search point has

\[
D=-340435200,
\]

\[
(s:m:n)=(254097487:-10481401502:555111200).
\]

Its canonical surface key is

\[
[16660,21456,78793]\mid[27593,78644].
\]

For this key,

\[
\frac{R}{ac}=\frac{320}{21},\qquad \frac{S}{ac}=\frac{21}{320}.
\]

A symmetry chart giving smaller curve parameters is

\[
(A,B,C)=(16660,78644,78793),
\]

\[
(s,m,n)=(1500518600,253638081,138777800).
\]

All 12 nondegenerate symmetry charts have identical absolute
\(G_2\)-invariants. A reduced minimal model is

\[
\begin{aligned}
y^2+(x^2+x)y={}&3703062294195264x^6
-360079374491052216x^5\\
&+8901721379573296848x^4
-5397945250386334945x^3\\
&-86737535708373850908x^2
+36346694984390901540x\\
&+43035470132681030400.
\end{aligned}
\]

Magma gives

\[
J(\mathbf Q)_{\mathrm{tors}}\simeq
(\mathbf Z/2\mathbf Z)^3\times\mathbf Z/12\mathbf Z.
\]

The torsion order is exactly 96. Frobenius root-power certificates at
\(p=37,127,131,179\) prove geometric simplicity, and the absolute
\(G_2\)-invariants differ from both old curves.

## 6. Reproduction

~~~sh
g++ -O3 -march=native -std=c++17 -fopenmp \
  code/search_22212_abc.cpp -o /tmp/search_22212_abc

/tmp/search_22212_abc 500 --threads 16
/tmp/search_22212_abc 20000 --threads 16

g++ -O3 -march=native -std=c++17 -fopenmp \
  code/search_22212_canonical.cpp -o /tmp/search_22212_canonical

/tmp/search_22212_canonical --T 10000 --threads 16 --first-chart-only

magma -b code/verify_22212_abc_hit.m
magma -b code/verify_22212_symmetric_model.m
~~~

The direct search supports --a-min and --a-max for deterministic sharding.
The new point alone is reproduced by

~~~sh
/tmp/search_22212_abc 20000 --threads 16 --a-min 4165 --a-max 4165
~~~

## 7. A plausible next search

Writing \(x=a/c\) and \(\lambda=R/(ac)\), after scaling by \(c\) the three
square conditions become

\[
w^2=x^2+1+x\lambda,\qquad
b^2=x^2+1+x/\lambda,
\]

\[
t^2=x^2+1+x(\lambda+\lambda^{-1}).
\]

The generic fiber over \(\lambda\) has genus 5, while the most obvious
low-genus values \(\lambda=\pm1\) lie on the excluded discriminant boundary.
A useful next route is therefore a targeted search over low-height
\(\lambda\), using local solubility and the three quadratic covers, rather
than another undirected cubic box.
