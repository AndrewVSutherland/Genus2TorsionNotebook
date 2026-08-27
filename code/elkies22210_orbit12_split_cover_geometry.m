/*
Geometry of the complementary-cubic splitting cover on a fixed fiber s
of the first orbit-12 radicand quotient.

The script verifies the generic equations and computes exact genera and
elliptic ranks at the known split fiber s=59/49.

Usage:

  magma -b code/elkies22210_orbit12_split_cover_geometry.m
*/

Q := Rationals();
A<s,q,x> := PolynomialRing(Q,3);

f := x^3+(1+q)*x^2+s*x-(1+q)*(q-s);
G := q*(s-q)*(2-q^2+(q+2)*s);
D := (1+q)^2*s^2-4*s^3+
     4*(1+q)^4*(q-s)-27*(1+q)^2*(q-s)^2-
     18*(1+q)^2*s*(q-s);
assert D eq -Resultant(f,Derivative(f,3),3);

// If x is a root, q satisfies this quadratic.  Its discriminant is H.
q_equation := q^2-(x^2+s-1)*q-(x+1)*(x^2+s);
assert q_equation eq -f;
H := x^4+4*x^3+2*(s+1)*x^2+4*s*x+(s+1)^2;
assert H eq (x^2+s-1)^2+4*(x+1)*(x^2+s);

// After selecting x, the other two roots split iff delta is a square.
delta := (q+x+1)^2-4*(s+x*(q+x+1));
assert delta eq q^2-2*q*x+2*q-3*x^2-2*x+1-4*s;

Ks<ss> := FunctionField(Q);
KsX<X> := PolynomialRing(Ks);
Hs := X^4+4*X^3+2*(ss+1)*X^2+4*ss*X+(ss+1)^2;
assert Discriminant(Hs) eq
       256*(16*ss^4-16*ss^3-8*ss^2-24*ss-11);

printf "ELKIES22210_ORBIT12_SPLIT_COVER_GEOMETRY\n";
printf "complementary_cubic %o\n", f;
printf "first_radicand %o\n", G;
printf "cubic_discriminant %o\n", D;
printf "one_root_quartic %o\n", H;
printf "remaining_quadratic_discriminant %o\n", delta;
printf "branch_overlap_resultant_factorization %o\n",
       Factorization(Resultant(D,G,2));

// The known smooth split-cubic point on the first-radicand cover.
s0 := Q!59/49;
qseed := Q!8/7;
vseed := Q!192/343;
roots := [Q!-9/7,Q!-5/7,Q!-1/7];
assert &and[Evaluate(f,[s0,qseed,a]) eq 0 : a in roots];
assert vseed^2 eq Evaluate(G,[s0,qseed,0]);
assert IsSquare(Evaluate(D,[s0,qseed,0]));

Pq<qq> := PolynomialRing(Q);
D0 := Evaluate(D,[s0,qq,0]);
G0 := qq*(s0-qq)*(2-qq^2+(qq+2)*s0);
H0 := qq^4+4*qq^3+2*(s0+1)*qq^2+4*s0*qq+(s0+1)^2;
assert IsSquarefree(D0) and IsSquarefree(G0) and IsSquarefree(H0);
assert GCD(D0,G0) eq 1;
assert #Factorization(D0) eq 1 and Degree(Factorization(D0)[1][1]) eq 5;

// Function-field towers.  R is the one-root curve.  Adding v imposes
// the first radicand; adding w splits the remaining quadratic.
Kx<xx> := FunctionField(Q);
KxT<T> := PolynomialRing(Kx);
H0fun := xx^4+4*xx^3+2*(s0+1)*xx^2+4*s0*xx+(s0+1)^2;
R<dd> := ext<Kx | T^2-H0fun>;
qR := (R!xx^2+s0-1+dd)/2;
GR := qR*(s0-qR)*(2-qR^2+(qR+2)*s0);
deltaR := (qR+R!xx+1)^2-4*(s0+(R!xx)*(qR+R!xx+1));

RT<V> := PolynomialRing(R);
Froot<v> := ext<R | V^2-GR>;
Ssplit<w> := ext<R | V^2-deltaR>;
FrootT<W> := PolynomialRing(Froot);
Hfull<ww> := ext<Froot | W^2-(Froot!deltaR)>;

assert Genus(R) eq 1;
assert Genus(Ssplit) eq 4;
assert Genus(Froot) eq 7;
assert Genus(Hfull) eq 19;
printf "fiber_s %o\n", s0;
printf "genus_one_root %o\n", Genus(R);
printf "genus_complete_split_without_G %o\n", Genus(Ssplit);
printf "genus_G_plus_one_root %o\n", Genus(Froot);
printf "genus_G_plus_complete_split %o\n", Genus(Hfull);

// The direct one-root base is unexpectedly high-rank at this fiber.
Croot := HyperellipticCurve(H0);
proot := Croot![Q!-9/7,Q!3/7,1];
Eroot, rootmap := EllipticCurve(Croot,proot);
rlo,rhi := RankBounds(Eroot);
assert rlo eq 4 and rhi eq 4;
Troot := TorsionSubgroup(Eroot);
assert #Troot eq 1;
printf "one_root_elliptic_rank_bounds %o %o\n", rlo,rhi;
printf "one_root_elliptic_torsion %o\n", Invariants(Troot);

// By contrast the first-radicand quotient itself has rank two.
Eq := EllipticCurve([Q | 0,s0^2-2*s0-2,0,
                     -4*s0^2*(s0+1),4*s0^2*(s0+1)^2]);
elo,ehi := RankBounds(Eq);
assert elo eq 2 and ehi eq 2;
TEq := TorsionSubgroup(Eq);
assert Invariants(TEq) eq [4];
printf "first_radicand_elliptic_rank_bounds %o %o\n", elo,ehi;
printf "first_radicand_elliptic_torsion %o\n", Invariants(TEq);

// A lower-genus necessary filter for complete splitting is the V4 cover
// v^2=G0(q), z^2=D0(q).  Its three hyperelliptic quotients have genera
// 1, 2, and 4, so the fiber product has genus 7.
CD := HyperellipticCurve(D0);
CGD := HyperellipticCurve(D0*G0);
assert Genus(CD) eq 2 and Genus(CGD) eq 4;
printf "discriminant_filter_quotient_genera 1 %o %o\n",
       Genus(CD),Genus(CGD);
printf "discriminant_polynomial_factorization_at_seed %o\n",
       Factorization(D0);
printf "seed_s %o seed_q %o seed_v %o seed_roots %o\n",
       s0,qseed,vseed,roots;
printf "DONE\n";

quit;
