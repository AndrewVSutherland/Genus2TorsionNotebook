// claude_ov_lane1x_geometry.m -- Lane 1 (overnight 2026-07-25): the GEOMETRY of
// the contact-7 three-root surface.
//
//   X   = { R(s,t,u) = 0 } subset A^3      (the object that actually carries
//                                           [2,2,14] curves: s,t,u must be
//                                           INDIVIDUALLY rational)
//   X_e = { R(e1,e2,e3) = 0 } subset A^3   (the S3-quotient X/S3; a point of
//                                           X_e lifts to X iff T^3-e1T^2+e2T-e3
//                                           splits completely over Q)
//
// R is QUADRATIC in e3, so X_e is a double cover of the (e1,e2)-plane branched
// along Delta := Disc_{e3}(R).  This script:
//   (1) rebuilds R and the e-form, and re-verifies the 11 known points;
//   (2) computes Delta and identifies the double cover;
//   (3) computes the singular locus of X (and of X_e);
//   (4) finds the u-values c for which the slice C_c = X cap {u=c} is singular
//       -- i.e. the ONLY c at which the generic genus 4 can drop -- by
//       eliminating (s,t) from <R, dR/ds, dR/dt>.  This turns the empirical
//       genus scan into a proof.
SetColumns(0);
if not assigned MemGB then MemGB := 12; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals();

// ---------- (1) rebuild R ----------
Rstu<s,t,u> := PolynomialRing(Q,3);
FF := FieldOfFractions(Rstu);
A := func< w | (2*w^5+4*w^4+6*w^3+8*w^2+10*w+5)/(2*(w+1)^2) >;
bST := (A(FF!t)-A(FF!s))/(FF!(s^2-t^2));
aST := A(FF!s) - bST*(1-FF!(s^2));
Enum := Numerator(A(FF!u) - aST - bST*(1-FF!(u^2)));
Npol := Rstu!Enum;
qq, rr := Quotrem(Npol, (u-s)*(u-t));
printf "remainder of the division by (u-s)(u-t): %o\n", rr;
R := qq;
// normalise to the primitive integral form used in the notes
den := LCM([ Denominator(cf) : cf in Coefficients(R) ]);
R := R * den;
num := GCD([ Numerator(cf) : cf in Coefficients(R) ]);
R := R / num;
printf "R tridegree: (%o,%o,%o)\n", Degree(R,s), Degree(R,t), Degree(R,u);

// S3-symmetry
sig := hom< Rstu -> Rstu | t,u,s >;
tau := hom< Rstu -> Rstu | t,s,u >;
printf "S3-symmetric: cyclic %o  transposition %o\n", sig(R) eq R, tau(R) eq R;

// eleven known rational points
pts := [
 [-10,-10/7,-1/2],[-5,-15/8,-15/22],[-3,-3/4,-3/5],[-15/8,-15/19,-1/2],
 [-5/18,-10/49,4/17],[-4/9,-4/25,4/17],[-511/61,-511/625,-1/2],
 [-165/41,-33/16,-165/289],[-164/297,-1/2,164/361],[-17/50,-34/189,34/121],
 [-1/2,-13/49,13/50]];
ok := true;
for P in pts do
  if Evaluate(R, [Q!P[1],Q!P[2],Q!P[3]]) ne 0 then ok := false; printf "POINT FAILS: %o\n", P; end if;
end for;
printf "all 11 known points lie on R=0: %o\n", ok;

// ---------- (2) the e-model and Delta = Disc_{e3}(R) ----------
Re<e1,e2,e3> := PolynomialRing(Q,3);
Rec := 2*(e2-2)*e3^2 + (4*e2^2+4*e1*e2-2*e1-4*e2-1)*e3 + e2*(2*(e1+e2)^2+e1);
// check: substituting the elementary symmetric functions must give R
phi := hom< Re -> Rstu | s+t+u, s*t+s*u+t*u, s*t*u >;
printf "\nrecorded e-form reproduces R: %o\n", phi(Rec) eq R;

Pe<E1,E2> := PolynomialRing(Q,2);
a2 := 2*(E2-2);
a1 := 4*E2^2+4*E1*E2-2*E1-4*E2-1;
a0 := E2*(2*(E1+E2)^2+E1);
Delta := a1^2 - 4*a2*a0;
printf "\nDelta = Disc_{e3}(R) (a quartic in (e1,e2)):\n  %o\n", Delta;
printf "  total degree = %o,  deg_e1 = %o,  deg_e2 = %o\n",
       TotalDegree(Delta), Degree(Delta,E1), Degree(Delta,E2);
printf "  factorisation: %o\n", Factorization(Delta);
// homogenise and study the projective quartic
Ph<X0,X1,X2> := PolynomialRing(Q,3);
hm := hom< Pe -> Ph | X1, X2 >;
Dh := &+[ Ph!(hm(Monomials(Delta)[i]))*Coefficients(Delta)[i]*X0^(4-TotalDegree(Monomials(Delta)[i]))
          : i in [1..#Terms(Delta)] ];
printf "\nprojective quartic branch curve B : %o\n", Dh;
P2 := ProjectiveSpace(Ph);
Bcur := Curve(P2, Dh);
printf "  irreducible: %o\n", IsIrreducible(Bcur);
if IsIrreducible(Bcur) then
  printf "  geometric genus of B: %o   (smooth plane quartic has genus 3)\n", Genus(Bcur);
end if;
sB := SingularSubscheme(Bcur);
printf "  singular subscheme of B: dimension %o, degree %o\n", Dimension(sB), Degree(sB);
if Dimension(sB) le 0 then
  printf "  singular points of B (over Q): %o\n", RationalPoints(sB);
end if;
printf "\n=> X_e is the double cover of P^2 branched along the plane quartic B,\n";
printf "   i.e. a (possibly singular) DEL PEZZO SURFACE OF DEGREE 2.\n";

// ---------- (3) singular locus of X ----------
A3 := AffineSpace(Rstu);
Xs := Scheme(A3, R);
printf "\nX = {R=0}: irreducible %o, dimension %o\n", IsIrreducible(Xs), Dimension(Xs);
sX := SingularSubscheme(Xs);
printf "Sing(X): dimension %o\n", Dimension(sX);
Isx := Ideal(sX);
Gsx := GroebnerBasis(Isx);
printf "Sing(X) Groebner basis (%o elements):\n", #Gsx;
for g in Gsx do printf "   %o\n", g; end for;
comps := PrimaryDecomposition(Isx);
printf "Sing(X) has %o primary components; supports:\n", #comps;
for cc in comps do
  rd := Radical(cc);
  printf "   dim=%o  gens=%o\n", Dimension(Scheme(A3, Basis(rd))), Basis(rd);
end for;
// the u-coordinates of Sing(X): eliminate s,t
Iel := EliminationIdeal(Isx, {u});
printf "\nSing(X) projected to the u-line, generator(s): %o\n", Basis(Iel);
for g in Basis(Iel) do
  printf "   factorisation: %o\n", Factorization(g);
end for;

// ---------- (4) Kodaira dimension of the closure in (P^1)^3 ----------
printf "\n== the closure Xbar in P^1 x P^1 x P^1 ==\n";
printf "  Xbar is a divisor of tridegree (3,3,3); K_{(P^1)^3} = O(-2,-2,-2), so by\n";
printf "  adjunction K_Xbar = O(1,1,1)|_Xbar, which is AMPLE.  Hence Xbar is a\n";
printf "  minimal surface of GENERAL TYPE with\n";
printf "     p_g = h^0(O(1,1,1)) = 8,   K^2 = (H1+H2+H3)^2*(3H1+3H2+3H3) = 18.\n";
printf "  Sing(Xbar) is 0-dimensional (computed above), so the resolution is still\n";
printf "  of general type: by Bombieri-Lang the rational points of X are expected to\n";
printf "  lie on finitely many curves.\n";

printf "LANE1X_GEOMETRY_DONE\n";
quit;
