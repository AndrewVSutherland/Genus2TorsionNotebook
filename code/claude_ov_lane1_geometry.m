// claude_ov_lane1_geometry.m -- Lane 1, overnight 2026-07-25.
//
// GLOBAL GEOMETRY OF THE CONTACT-7 THREE-ROOT SURFACE.
//
// Two models must be kept apart:
//   * the (s,t,u) model X = {R = 0} -- a rational point here IS a [2,2,14]
//     curve (three rational Weierstrass v-values);
//   * the elementary-symmetric model Xe = {R = 0} in (e1,e2,e3) -- a rational
//     point here only gives a [2,2,14] curve when T^3-e1T^2+e2T-e3 SPLITS
//     COMPLETELY over Q.  X -> Xe is the degree-6 S3-quotient.
//
// (A) X-bar inside P^1 x P^1 x P^1 is a divisor of tridegree (3,3,3), so by
//     adjunction K_Xbar = O(1,1,1)|_Xbar, which is AMPLE:  X is of GENERAL
//     TYPE (K^2 = 18) as soon as its singularities are canonical.  This
//     script computes the singular locus of X-bar and its dimension/degree.
//
// (B) R is QUADRATIC in e3, so Xe is the double cover of the (e1,e2)-plane
//     branched along Disc_{e3}(R).  This script computes that discriminant,
//     its degree and smoothness -- a smooth plane quartic branch locus makes
//     Xe a DEL PEZZO SURFACE OF DEGREE 2 (rational over Qbar, unirational
//     over Q once it has a general rational point).
//
// (C) The u = c slice fibration: reported by code/claude_ov_lane1_fibration.m.
//
// Usage: magma -b code/claude_ov_lane1_geometry.m

SetColumns(0);
if not assigned MemGB then MemGB := 24; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals(); Z := Integers();

/* ---- rebuild R ---- */
S3<sv,tv,uv> := PolynomialRing(Q,3);
FS := FieldOfFractions(S3);
Afun := func<w0 | (2*w0^5+4*w0^4+6*w0^3+8*w0^2+10*w0+5)/(2*(w0+1)^2)>;
bST := (Afun(FS!tv)-Afun(FS!sv))/(FS!(sv^2-tv^2));
aST := Afun(FS!sv) - bST*(1-FS!(sv^2));
Enum := Afun(FS!uv) - aST - bST*(1-FS!(uv^2));
ok, Rfull := IsDivisibleBy(S3!Numerator(Enum), S3!((uv-sv)*(uv-tv))); assert ok;
dl := LCM([Denominator(cc) : cc in Coefficients(Rfull)]);
Rfull := S3!(dl*Rfull);
Rfull := S3!(Rfull div (S3!GCD([Z!cc : cc in Coefficients(Rfull)])));
printf "R tridegree (%o,%o,%o)\n", Degree(Rfull,1),Degree(Rfull,2),Degree(Rfull,3);

/* ============ (A) X-bar in P^1 x P^1 x P^1 ============ */
printf "\n=== (A) X-bar in P1xP1xP1 : tridegree (3,3,3) ===\n";
PPP<a0,a1,b0,b1,c0v,c1v> := ProductProjectiveSpace(Q, [1,1,1]);
// trihomogenise R: s = a0/a1, t = b0/b1, u = c0v/c1v
CR := CoordinateRing(PPP);
Rtri := CR!0;
mons := Monomials(Rfull); cfs := Coefficients(Rfull);
for i in [1..#mons] do
  ex := Exponents(mons[i]);
  Rtri +:= (CR!cfs[i]) * a0^ex[1]*a1^(3-ex[1]) * b0^ex[2]*b1^(3-ex[2]) * c0v^ex[3]*c1v^(3-ex[3]);
end for;
Xbar := Scheme(PPP, Rtri);
printf "dim X-bar = %o, irreducible = %o\n", Dimension(Xbar), IsIrreducible(Xbar);
Sg := SingularSubscheme(Xbar);
printf "singular subscheme: dim = %o\n", Dimension(Sg);
if Dimension(Sg) ge 0 then
  cps := IrreducibleComponents(ReducedSubscheme(Sg));
  printf "singular locus has %o irreducible components:\n", #cps;
  // NOTE (2026-07-25 fix): RationalPoints/RationalPointsGeneric REFUSES a scheme
  // living in a product of projective spaces ("Argument must lie in affine or
  // ordinary projective space ...").  Push the component through the Segre
  // embedding into an ordinary P^7 first; that is where point enumeration works.
  for K in cps do
    if Dimension(K) eq 0 then
      pts := "n/a";
      try
        pts := Sprint(RationalPoints(SegreEmbedding(K)));
      catch ee
        // last resort: read the points off the (already reduced, 0-dimensional) ideal
        pts := Sprint(GroebnerBasis(Ideal(K)));
      end try;
      printf "   dim 0  points (Segre model): %o\n", pts;
    else
      printf "   dim %o  positive-dimensional\n", Dimension(K);
    end if;
    printf "     equations: %o\n", DefiningEquations(K);
  end for;
end if;
printf "adjunction: K_Xbar = (-2,-2,-2)+(3,3,3) = O(1,1,1)|_Xbar, AMPLE; K^2 = 18\n";
printf "=> X is of GENERAL TYPE if the singularities above are canonical.\n";

/* ============ (B) the e-model: double cover of the (e1,e2)-plane ============ */
printf "\n=== (B) e-model: R quadratic in e3, branch = Disc_e3(R) ===\n";
Re<E1,E2,E3> := PolynomialRing(Q,3);
Rec := 2*(E2-2)*E3^2 + (4*E2^2+4*E1*E2-2*E1-4*E2-1)*E3 + E2*(2*(E1+E2)^2+E1);
printf "R (e-model) = %o\n", Rec;
Aq := 2*(E2-2); Bq := 4*E2^2+4*E1*E2-2*E1-4*E2-1; Cq := E2*(2*(E1+E2)^2+E1);
Disc := Bq^2 - 4*Aq*Cq;
P2e<e1,e2> := PolynomialRing(Q,2);
Dq := Evaluate(Disc, [P2e.1, P2e.2, 0]);
printf "Disc_e3(R) = %o\n", Dq;
printf "  total degree = %o\n", TotalDegree(Dq);
printf "  factorisation = %o (unit %o)\n", Factorization(Dq),
   Dq div &*[t[1]^t[2] : t in Factorization(Dq)];
// projective closure & smoothness of the branch curve
PR2<z0,z1,z2> := ProjectiveSpace(Q,2);
dtot := TotalDegree(Dq);
CR2 := CoordinateRing(PR2);
Dhom := CR2!0;
for i in [1..#Monomials(Dq)] do
  ex := Exponents(Monomials(Dq)[i]);
  Dhom +:= (CR2!Coefficients(Dq)[i]) * z0^ex[1]*z1^ex[2]*z2^(dtot-ex[1]-ex[2]);
end for;
Bc := Curve(PR2, Dhom);
printf "projective branch curve: degree %o, geometric genus %o, singular points %o\n",
   Degree(Bc), (IsIrreducible(Bc) select Genus(Bc) else -1), SingularPoints(Bc);
printf "arithmetic genus of a smooth plane quartic = 3\n";
if dtot eq 4 and #SingularPoints(Bc) eq 0 then
  printf "=> the (e1,e2)-double cover branched along a SMOOTH plane quartic\n";
  printf "   is a DEL PEZZO SURFACE OF DEGREE 2 : rational over Qbar, unirational\n";
  printf "   over Q as soon as it carries a rational point in general position.\n";
else
  printf "=> branch curve is singular / not a quartic; the double cover is a\n";
  printf "   (possibly singular) degree-2 del Pezzo or a weak dP2.\n";
end if;
// the eleven points, pushed to the e-model, lie on Xe:
elevens := [
 [Q| -10,-10/7,-1/2], [Q| -5,-15/8,-15/22], [Q| -3,-3/4,-3/5], [Q| -15/8,-15/19,-1/2],
 [Q| -5/18,-10/49,4/17], [Q| -4/9,-4/25,4/17], [Q| -511/61,-511/625,-1/2],
 [Q| -165/41,-33/16,-165/289], [Q| -164/297,-1/2,164/361], [Q| -17/50,-34/189,34/121],
 [Q| -1/2,-13/49,13/50] ];
printf "\ne-images of the eleven known X(Q)-points, and Disc value (square iff e3 rational):\n";
for pt in elevens do
  s0,t0,u0 := Explode(pt);
  E1v := s0+t0+u0; E2v := s0*t0+s0*u0+t0*u0; E3v := s0*t0*u0;
  dv := Evaluate(Dq, [E1v,E2v]);
  bsq, rt := IsSquare(dv);
  printf "  (e1,e2,e3) = (%o, %o, %o)   R=%o   Disc=%o square=%o\n",
     E1v, E2v, E3v, Evaluate(Rec,[E1v,E2v,E3v]), dv, bsq;
end for;
printf "\nNOTE: Xe being (uni)rational does NOT give [2,2,14] curves -- the cubic\n";
printf "T^3-e1T^2+e2T-e3 must SPLIT COMPLETELY.  That splitting cover is X itself,\n";
printf "which is of general type (A).  This is the precise obstruction.\n";
printf "LANE1_GEOMETRY_DONE\n";
quit;
