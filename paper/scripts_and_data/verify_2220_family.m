// verify_2220_family.m -- the [2,2,20] subsection ("A second order-80
// group") of the paper.
//
// Verifies, in the paper's coordinates:
//   (1) over Q(t), the family
//         F_t : y^2 = x(x-2)(x^4+(4t+2)x^3+8t(t+1)x^2+8t(t+1)^2x+4(t+1)^4)
//       has the class [inf_+ - (0,0)] of exact order 5 and the class
//       [(-t-1,(t+1)^2(t+3)) - (0,0)] of exact order 20;
//   (2) with t = -(z^4+4z+4)/(z^4+4z^3+8z^2+8z+4), the quartic factor
//       acquires a rational root: it factors as (linear)(cubic) over Q(z);
//   (3) the condition that the remaining cubic has a root defines a plane
//       curve Ytilde (degree-3 cover of P^1_z) of geometric genus 2 whose
//       normalization is isomorphic to
//         Y : y^2 + (x^2+1)y = x^5 - x
//       of conductor 2528 (LMFDB 2528.a.20224.1);
//   (4) Jac(Y) has Mordell-Weil rank exactly 1 (RankBounds = [1,1]) and
//       torsion Z/4; Magma's Chabauty method, from an infinite-order
//       generator, provably computes all rational points of Y (9 points);
//   (5) transporting Y(Q) back to Ytilde (plus the rational singular points
//       of the plane model and its points at infinity) gives z in
//       {-1, 0, -1/7, -7/9, infinity}; the values -1, 0, infinity make the
//       sextic degenerate, so z = -1/7 and z = -7/9 are the ONLY
//       nondegenerate rational fibers;
//   (6) both surviving fibers are isomorphic to the displayed curve
//       Eq. (eq:2220), whose exact torsion is [2,2,20] -- so that curve is
//       the unique [2,2,20] member of this one-parameter family.
// (End(Jac_Qbar) = Z for the eq:2220 curve is certified in
// verify_simple_certificates.m; it is the same curve as the Table 1
// [2,2,20] row, which is checked here as well.)
// Run from this directory:  magma -b verify_2220_family.m
SetColumns(0);
SetMemoryLimit(8*10^9);
Q := Rationals();
t0c := Cputime();

// ---- (1) the family over Q(t) ----
Kt<t> := FunctionField(Q);
Pt<xt> := PolynomialRing(Kt);
fam := xt*(xt-2)*(xt^4 + (4*t+2)*xt^3 + 8*t*(t+1)*xt^2
                  + 8*t*(t+1)^2*xt + 4*(t+1)^4);
C := HyperellipticCurve(fam);
J := Jacobian(C);
inf := PointsAtInfinity(C);
P0 := C![0, 0];
D5 := J![inf[1], P0];
assert D5 ne J!0 and 5*D5 eq J!0;                       // exact order 5
D20 := J![C![-t-1, (t+1)^2*(t+3)], P0];
assert 20*D20 eq J!0 and 10*D20 ne J!0 and 4*D20 ne J!0; // exact order 20
printf "(1) generic orders over Q(t): [inf+ - (0,0)] = 5, [(-t-1,.) - (0,0)] = 20\n";

// ---- (2) the z-parameterization of the quartic-root condition ----
Kz<z> := FunctionField(Q);
Pz<x> := PolynomialRing(Kz);
tz := -(z^4+4*z+4)/(z^4+4*z^3+8*z^2+8*z+4);
quart := x^4 + (4*tz+2)*x^3 + 8*tz*(tz+1)*x^2
         + 8*tz*(tz+1)^2*x + 4*(tz+1)^4;
fac := Factorization(quart);
assert Sort([ Degree(ff[1]) : ff in fac ]) eq [1, 3];
cub := [ ff[1] : ff in fac | Degree(ff[1]) eq 3 ][1];
printf "(2) at t = t(z) the quartic factors as (linear)(cubic) over Q(z)\n";

// ---- (3) the cubic-root cover and its normalization ----
R2 := PolynomialRing(Q, 2);
den := LCM([ Denominator(Coefficient(cub, i)) : i in [0..3] ]);
cnum := &+[ Evaluate(Numerator(Coefficient(cub, i)*(Kz!den)), R2.1) * R2.2^i
            : i in [0..3] ];
A2 := AffineSpace(Q, 2);
Ytil := Curve(A2, Evaluate(cnum, [A2.1, A2.2]));
assert Genus(Ytil) eq 2;
isodd, H, mp := IsHyperelliptic(Ytil);
assert isodd;
Hs, sm := SimplifiedModel(H);
Pu<u> := PolynomialRing(Q);
Y := HyperellipticCurve(u^5 - u, u^2 + 1);       // y^2 + (x^2+1)y = x^5 - x
assert Conductor(Y) eq 2528;
Yp, smY := SimplifiedModel(Y);
iso, phi := IsIsomorphic(Hs, Yp);
assert iso;
printf "(3) the cover is a singular plane genus-2 curve with normalization Y (conductor 2528)\n";

// ---- (4) rank 1 and Chabauty ----
JY := Jacobian(Yp);
rlo, rhi := RankBounds(JY);
assert rlo eq 1 and rhi eq 1;
assert Invariants(AbelianGroup(TorsionSubgroup(JY))) eq [4];
gen := JY!0;
for D in Points(JY : Bound := 500) do
    if D ne JY!0 and Order(D) eq 0 then gen := D; break; end if;
end for;
assert gen ne JY!0;
YQ := Chabauty(gen);
assert #YQ eq 9;
printf "(4) rank(Jac Y) = 1, torsion Z/4; Chabauty: Y(Q) has exactly %o points\n", #YQ;

// ---- (5) transport back to the z-line ----
F := mp * sm * phi;                              // Ytil -> Yp, birational
// Every smooth rational point of the plane model lifts uniquely (hence
// rationally) to the normalization, so the rational points of Ytil are
// contained in: the rational preimages of Y(Q) under F (scheme pullback --
// a superset, as it includes the base locus), the rational singular points,
// and the points at infinity of the projective closure.
cand := {};
for pt in YQ do
    for r in RationalPoints(pt @@ F) do
        Include(~cand, Coordinates(r));
    end for;
end for;
for s in SingularPoints(Ytil) do
    Include(~cand, Coordinates(s));
end for;
zvals := { c[1] : c in cand };
assert zvals eq { Q!-1, Q!0, Q!-1/7, Q!-7/9 };
// t(z) is defined at every rational z (the denominator has no rational
// root), and a rational root of the cubic is an affine point, so the
// candidates above are complete; z = infinity has t = -1 (degenerate)
Pu2<w> := PolynomialRing(Q);
assert #Roots(w^4+4*w^3+8*w^2+8*w+4) eq 0;
// the degenerate fibers: t = -1 (from z = 0, -1, infinity)
famt := func< tv | Pu!(u*(u-2)*(u^4 + (4*tv+2)*u^3 + 8*tv*(tv+1)*u^2
                   + 8*tv*(tv+1)^2*u + 4*(tv+1)^4)) >;
for zv in [ Q!-1, Q!0 ] do
    tv := Evaluate(tz, zv);
    assert tv eq -1;
    assert Discriminant(famt(tv)) eq 0;          // degenerate sextic
end for;
printf "(5) all rational points of the cover have z in {-1, 0, -1/7, -7/9, oo}; z = -1, 0, oo give t = -1, a degenerate fiber\n";

// ---- (6) the two surviving fibers are Eq. (eq:2220) ----
Pq<xq> := PolynomialRing(Q);
eq2220 := SimplifiedModel(HyperellipticCurve(
    -391671*xq^6+1894851*xq^5+6846924*xq^4-15133525*xq^3
    +3904068*xq^2+2625336*xq+254016, xq^2+xq));
table_row := SimplifiedModel(HyperellipticCurve(
    -(xq-1)*(6*xq+1)*(2*xq+7)*(6217*xq+1008)*(21*xq^2-161*xq+144)));
assert IsIsomorphic(eq2220, table_row);          // same curve as Table 1 row
for zv in [ Q!-1/7, Q!-7/9 ] do
    tv := Evaluate(tz, zv);
    fz := famt(tv);
    assert Discriminant(fz) ne 0;                // nondegenerate
    dz := LCM([ Denominator(co) : co in Coefficients(fz) ]);
    Cz := SimplifiedModel(HyperellipticCurve(dz^2*fz));   // integral model
    assert IsIsomorphic(Cz, eq2220);
    assert Invariants(AbelianGroup(TorsionSubgroup(Jacobian(Cz))))
        eq [2, 2, 20];
    printf "(6) z = %o: fiber isomorphic to Eq.(eq:2220), exact torsion [2,2,20]\n", zv;
end for;

printf "[2,2,20] FAMILY VERIFIED: Eq.(eq:2220) is the unique [2,2,20] member (%.1o s)\n", Cputime() - t0c;
quit;
