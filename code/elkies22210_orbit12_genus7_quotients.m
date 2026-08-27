/*
Arithmetic of the three quotient factors of the genus-7 V4 cover at
s=59/49:

  C_G  : v^2 = G(q),
  C_D  : z^2 = D(q),
  C_GD : y^2 = G(q)D(q).

The class-group work in the hyperelliptic 2-descents is run under GRH.
All polynomial identities, genera, displayed rational points, orders, and
the elliptic Mordell--Weil computation are unconditional.

Usage:

  magma -b code/elkies22210_orbit12_genus7_quotients.m
*/

SetClassGroupBounds("GRH");
Q := Rationals();
P<q> := PolynomialRing(Q);
s := Q!59/49;

G := q*(s-q)*(2-q^2+(q+2)*s);
D := (1+q)^2*s^2-4*s^3+4*(1+q)^4*(q-s)
     -27*(1+q)^2*(q-s)^2-18*(1+q)^2*s*(q-s);
FD := 7^6*D;
FGD := 7^10*G*D;
assert &and[Denominator(a) eq 1:a in Coefficients(FD)];
assert &and[Denominator(a) eq 1:a in Coefficients(FGD)];
assert IsSquarefree(FD) and IsSquarefree(FGD);
assert FGD eq
    1129900996*q^9-7188706847*q^8+5744565372*q^7
    +33182708370*q^6-41331102120*q^5-48218619491*q^4
    +64606391416*q^3+26876076480*q^2-35080459776*q;

CD := HyperellipticCurve(FD);
CGD := HyperellipticCurve(FGD);
JD := Jacobian(CD);
JGD := Jacobian(CGD);
assert Genus(CD) eq 2 and Genus(CGD) eq 4;

// The first quotient is the already certified rank-two elliptic curve.
E := EllipticCurve([Q|0,s^2-2*s-2,0,
                    -4*s^2*(s+1),4*s^2*(s+1)^2]);
Emin,minmap := MinimalModel(E);
elo,ehi := RankBounds(Emin);
MW,mwmap := MordellWeilGroup(Emin);
assert elo eq 2 and ehi eq 2 and Invariants(MW) eq [4,0,0];

// Two rational points on C_D give independent divisor classes.
A8 := JD![q-Q!8/7,Q!128];
As := JD![q-s,Q!590/7];
assert Order(A8) eq 0 and Order(As) eq 0;
HD := HeightPairingMatrix([A8,As]);
assert Determinant(HD) gt 0;

// The seed also supplies a non-torsion class on C_GD.
AGD := JGD![q-Q!8/7,Q!24576/7];
assert Order(AGD) eq 0;

// A bounded search on the integral genus-4 model.  These are the three
// Weierstrass q-values infinity, 0, s and the two signs of the seed.
CGDPointBound := 1000;
searched_GD := Points(CGD : Bound:=CGDPointBound);
expected_GD := {
    CGD![Q!1,Q!0,Q!0],
    CGD![Q!0,Q!0,Q!1],
    CGD![Q!59,Q!0,Q!49],
    CGD![Q!8,Q!59006976,Q!7],
    CGD![Q!8,Q!-59006976,Q!7]
};
assert searched_GD eq expected_GD;

// Hyperelliptic 2-descent upper bounds.  These use the GRH class-group mode
// declared above; the explicit lower bounds do not.
dlo,dhi := RankBounds(JD);
gdlo,gdhi := RankBounds(JGD);
assert dhi eq 3 and gdhi eq 2;

// The full fake 2-Selmer set of C_D consists of five classes.  Three are
// represented by the known x-coordinates infinity, 8/7, and 59/49.
Sel,toSel := TwoCoverDescent(CD);
A := Domain(toSel);
theta := A.1;
h_inf := toSel(A!1);
h_8 := toSel(A!(Q!8/7)-theta);
h_s := toSel(A!s-theta);
assert #Sel eq 5;
assert h_inf in Sel and h_8 in Sel and h_s in Sel;
assert #{h_inf,h_8,h_s} eq 3;

AutD,autmap,action := AutomorphismGroup(CD);
assert #AutD eq 2; // only identity and the hyperelliptic involution.

printf "ELKIES22210_ORBIT12_GENUS7_QUOTIENTS\n";
printf "s %o\n",s;
printf "G %o\n",G;
printf "D %o\n",D;
printf "D_factorization %o\n",Factorization(D);
printf "GD_factor_degrees %o\n",[Degree(a[1]):a in Factorization(G*D)];
printf "quotient_genera 1 %o %o\n",Genus(CD),Genus(CGD);
printf "E_rank_bounds %o %o torsion %o\n",elo,ehi,[4];
printf "CD_known_non_torsion_orders %o %o\n",Order(A8),Order(As);
printf "CD_height_pairing %o\n",HD;
printf "CD_height_determinant %o\n",Determinant(HD);
printf "CGD_seed_order %o\n",Order(AGD);
printf "CGD_point_search_bound %o\n",CGDPointBound;
printf "CGD_point_search_points %o\n",searched_GD;
printf "class_group_mode GRH\n";
printf "CD_2descent_rank_bounds_from_intrinsic %o %o\n",dlo,dhi;
printf "CD_rank_bounds_with_explicit_points 2 %o\n",dhi;
printf "CGD_2descent_rank_bounds_from_intrinsic %o %o\n",gdlo,gdhi;
printf "CGD_rank_bounds_with_seed 1 %o\n",gdhi;
printf "genus7_total_rank_interval 5 7\n";
printf "CD_fake_2Selmer_classes %o\n",#Sel;
printf "CD_known_descent_classes_distinct %o\n",#{h_inf,h_8,h_s};
printf "CD_automorphism_group_order %o\n",#AutD;
printf "DONE\n";
quit;
