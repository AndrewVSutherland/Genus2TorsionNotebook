/*
Arithmetic of the genus-2 discriminant quotient at the known orbit-12
seed fiber s=59/49.

This script deliberately uses unconditional class-group bounds.  On the
current machine the 2-descent takes about two minutes and the point search
to height 10^6 takes about three minutes.

Usage:

  magma -b code/elkies22210_orbit12_genus2_discriminant.m
*/

SetColumns(0);
SetSeed(1);

Q := Rationals();
P<q> := PolynomialRing(Q);
s := Q!59/49;

// Discriminant in U of
// U^3+(1+q)U^2+sU-(1+q)(q-s).
D := -4*s^3-8*s^2*q^2-16*s^2*q-8*s^2
     -4*s*q^4+20*s*q^3+48*s*q^2+20*s*q-4*s
     +4*q^5-11*q^4-30*q^3-11*q^2+4*q;
assert D eq 4*q^5-775/49*q^4-290/49*q^3+84509/2401*q^2
               +11728/2401*q-2752704/117649;
assert IsSquarefree(D);

// Put X=14q and Y=1372z.  Since D has leading coefficient 4,
// this gives the primitive integral odd-degree model Y^2=F(X).
PX<X> := PolynomialRing(Q);
F := 14*X^5-775*X^4-4060*X^3+338036*X^2
     +656768*X-44043264;
assert F eq 14^6*Evaluate(D/4,X/14);
assert GCD([Integers()!Coefficient(F,i) : i in [0..5]]) eq 1;

C := HyperellipticCurve(F);
assert Genus(C) eq 2;

// Magma's globally minimal Weierstrass model, followed by SL_2(Z)
// reduction.  The second model is still globally minimal.
Cglob, toglob := MinimalWeierstrassModel(C);
Cmin, reduce_map := ReducedModel(Cglob : Smallest:=true, Height:=true);
fmin,hmin := HyperellipticPolynomials(Cmin);
assert fmin eq 7*X^5+86*X^4-2228*X^3-29613*X^2+7680*X+4096;
assert hmin eq X^2+X;

// A direct map from z^2=D(q) to the reduced globally minimal model is
//   x=7q-8,
//   y=(343z-49q^2+105q-56)/2.
A<qa,za> := PolynomialRing(Q,2);
ua := 7*qa-8;
wa := (343*za-49*qa^2+105*qa-56)/2;
Da := Evaluate(D,qa);
fa := Evaluate(fmin,ua);
ha := Evaluate(hmin,ua);
rel := wa^2+ha*wa-fa;
assert rel-Coefficient(rel,2,2)*(za^2-Da) eq 0;

// The seed, the G=0 boundary q=s, and infinity.
seed_plus := Cmin![Q!0,Q!64,Q!1];
seed_minus := Cmin![Q!0,Q!-64,Q!1];
boundary_plus := Cmin![Q!3/7,Q!2050/49,Q!1];
boundary_minus := Cmin![Q!3/7,Q!-2080/49,Q!1];
infinity := Cmin![Q!1,Q!0,Q!0];
known_min := {infinity,seed_plus,seed_minus,boundary_plus,boundary_minus};
assert #known_min eq 5;

// On the primitive y^2=F(X) model these points are especially useful
// for constructing divisor classes.
seed_C := C![Q!16,Q!512,Q!1];
boundary_C := C![Q!118,Q!115640,Q!7];
J := Jacobian(C);
Dseed := J![X-16,512];
Dboundary := J![X-118/7,2360/7];
assert Order(Dseed) eq 0 and Order(Dboundary) eq 0;
reg := Regulator([Dseed,Dboundary]);
assert reg gt 20 and reg lt 22;

printf "ELKIES22210_ORBIT12_GENUS2_DISCRIMINANT\n";
printf "s %o\n",s;
printf "D(q) %o\n",D;
printf "primitive_integral_F(X) %o\n",F;
printf "globally_minimal_reduced_model y^2+(%o)y=%o\n",hmin,fmin;
printf "minimal_discriminant %o\n",Discriminant(Cmin);
printf "minimal_discriminant_factorization %o\n",
       Factorization(Integers()!Discriminant(Cmin));
printf "known_minimal_points %o\n",known_min;
printf "known_divisor_regulator %.20o\n",reg;

// Saturate the explicit rank-two subgroup at all primes <=31.
B := [Dseed,Dboundary];
sat_primes := [2,3,5,7,11,13,17,19,23,29,31];
for p in sat_primes do
    Bp := Saturation(B,p);
    assert #Bp eq 2;
    assert Abs(Regulator(Bp)-Regulator(B)) lt 10^-20;
    B := Bp;
end for;
printf "explicit_subgroup_p_saturated %o\n",sat_primes;

// This is an unconditional 2-descent.  The reduced minimal model lets
// RankBounds find the two displayed independent classes automatically.
Jmin := Jacobian(Cmin);
time rlo,rhi := RankBounds(Jmin);
assert rlo eq 2 and rhi eq 3;
T,incT := TorsionSubgroup(J);
assert #T eq 1;
S2,selmap := TwoSelmerGroup(Jmin);
assert Invariants(S2) eq [2,2,2];
square_sha_parity := HasSquareSha(Jmin);
assert square_sha_parity;
printf "rank_bounds_unconditional %o %o\n",rlo,rhi;
printf "torsion_invariants %o\n",Invariants(T);
printf "two_selmer_invariants %o\n",Invariants(S2);
printf "has_square_sha_assuming_sha_finite %o\n",square_sha_parity;
printf "conditional_on_sha_finite_rank 3\n";

// A single good reduction proves Q-simplicity: an abelian surface that
// split over Q would have reducible Frobenius polynomial at every good
// prime.
k := GF(11);
Pk<t> := PolynomialRing(k);
C11 := HyperellipticCurve(Pk!F);
L11 := LPolynomial(C11);
assert #Factorization(PolynomialRing(Q)!L11) eq 1;
printf "L_polynomial_p11 %o\n",L11;
printf "L_polynomial_p11_factorization_over_Q %o\n",
       Factorization(PolynomialRing(Q)!L11);

// This is a bounded search, not a proof of completeness.
PointBound := 10^6;
time searched := Points(C : Bound:=PointBound);
expected := {
    C![Q!1,Q!0,Q!0],
    C![Q!16,Q!512,Q!1], C![Q!16,Q!-512,Q!1],
    C![Q!118,Q!115640,Q!7], C![Q!118,Q!-115640,Q!7]
};
assert searched eq expected;
printf "point_search_bound %o\n",PointBound;
printf "point_search_points %o\n",searched;
printf "classical_chabauty_available false rank_lower_%o_genus_%o\n",
       rlo,Genus(Cmin);
printf "DONE\n";

quit;
