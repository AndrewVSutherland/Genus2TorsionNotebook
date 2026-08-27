// claude_ov_lane1_fibration.m -- Lane 1, overnight 2026-07-25.
//
// GEOMETRY OF THE CONTACT-7 THREE-ROOT SURFACE X = {R(s,t,u)=0} subset A^3.
//
// R is S3-symmetric of tridegree (3,3,3).  Fix u = c and symmetrise the
// remaining pair via p = s+t, q = st.  Then
//
//   Rpq(p,q;c) = 2c(c+1)^2 p^3 + 2(c+1)^2(2c+1) p^2 q + 2(c+1)^2(c+2) p q^2
//              + 2(c+1)^2 q^3 + c(2c+1)^2 p^2 + (4c^3+4c^2+2c+1) p q
//              + c^2(2c+1) p
//
// a PLANE CUBIC through O = (p,q) = (0,0) for every c.  So
//     Y := {Rpq = 0} subset A^2 x A^1_c
// is an ELLIPTIC SURFACE over the c-line with the zero section O, and
//     X --2:1--> Y   branched over p^2 = 4q  (i.e. s = t).
// A rational point of X = a rational point of Y at which p^2-4q is a square.
//
// This script:
//   (1) rebuilds Rpq and checks the factored coefficient shape;
//   (2) builds the generic fibre over Q(c) as an elliptic curve, prints its
//       Weierstrass model, discriminant, j-invariant, and the c-values with
//       SINGULAR fibres (the candidate extra rational slices);
//   (3) for every rational singular c, examines the slice: is the cubic
//       nodal/cuspidal, rational?  and if so what is the (s,t) double cover?
//   (4) reports the same data for a list of rational c of small height, and
//       computes the geometric genus of the (s,t) slice curve.
//
// Usage: magma -b code/claude_ov_lane1_fibration.m

SetColumns(0);
if not assigned MemGB then MemGB := 12; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals(); Z := Integers();

/* ---------------- rebuild R and Rpq ---------------- */
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
printf "R tridegree (%o,%o,%o), S3-symmetric %o\n", Degree(Rfull,1),Degree(Rfull,2),Degree(Rfull,3),
  Evaluate(Rfull,[tv,uv,sv]) eq Rfull;

// Rpq: substitute t = p - s and reduce mod s^2 - p s + q
Rc<pp,qq,cc0> := PolynomialRing(Q,3);
Rst<ss> := PolynomialRing(Rc);
Rsub := Evaluate(Rfull, [Rst!ss, Rst!(pp) - Rst!ss, Rst!cc0]);
Rrem := Rsub mod (ss^2 - pp*ss + qq);
assert Degree(Rrem) le 0;
Rpq := Rc!Coefficient(Rrem,0);
Rpq := Rc!(Rpq div (Rc!GCD([Z!cc : cc in Coefficients(Rpq)])));
printf "Rpq = %o\n", Rpq;
// factored coefficient check
PC<c> := PolynomialRing(Q);
chkR := 2*cc0*(cc0+1)^2*pp^3 + 2*(cc0+1)^2*(2*cc0+1)*pp^2*qq + 2*(cc0+1)^2*(cc0+2)*pp*qq^2
      + 2*(cc0+1)^2*qq^3 + cc0*(2*cc0+1)^2*pp^2 + (4*cc0^3+4*cc0^2+2*cc0+1)*pp*qq + cc0^2*(2*cc0+1)*pp;
printf "factored-coefficient form matches: %o\n", chkR eq Rpq;

/* ---------------- generic fibre over Q(c) ---------------- */
Fc<cg> := FunctionField(Q);
PP2<X0,X1,X2> := ProjectiveSpace(Fc,2);   // X0=p, X1=q, X2=z
cub := 2*cg*(cg+1)^2*X0^3 + 2*(cg+1)^2*(2*cg+1)*X0^2*X1 + 2*(cg+1)^2*(cg+2)*X0*X1^2
     + 2*(cg+1)^2*X1^3 + cg*(2*cg+1)^2*X0^2*X2 + (4*cg^3+4*cg^2+2*cg+1)*X0*X1*X2
     + cg^2*(2*cg+1)*X0*X2^2;
Ccub := Curve(PP2, cub);
Opt := Ccub ! [0,0,1];
printf "generic fibre is a plane cubic through O=(0:0:1); building Weierstrass model...\n";
Egen, mp := EllipticCurve(Ccub, Opt);
printf "Egen = %o\n", Egen;
Dg := Discriminant(Egen);
jg := jInvariant(Egen);
printf "Discriminant(Egen) = %o\n", Dg;
printf "jInvariant(Egen)   = %o\n", jg;

numD := PC!Numerator(Fc!Dg); denD := PC!Denominator(Fc!Dg);
printf "disc numerator factored: %o\n", Factorization(numD);
printf "disc denominator factored: %o\n", Factorization(denD);
printf "SINGULAR FIBRES at c = roots of the numerator:\n";
for tt in Factorization(numD) do
  rr := Roots(tt[1]);
  printf "  factor %o (mult %o) rational roots %o\n", tt[1], tt[2], [r[1] : r in rr];
end for;

/* ---------------- slice analysis ---------------- */
// For a given rational c, the (s,t) curve R(s,t,c)=0 is a (3,3) curve.
// Compute (a) whether the plane cubic Rpq(.,.,c) is singular and where,
//         (b) the geometric genus of the (s,t) slice curve.
SliceReport := procedure(c0)
  A2<pa,qa> := AffineSpace(Q,2);
  fc := Evaluate(Rpq, [A2.1, A2.2, c0]);
  if fc eq 0 then printf "c=%o : Rpq identically zero\n", c0; return; end if;
  Ca := Curve(A2, fc);
  sing := SingularPoints(Ca);
  // (s,t) slice
  A2b<sa,ta> := AffineSpace(Q,2);
  fs := Evaluate(Rfull, [A2b.1, A2b.2, c0]);
  if fs eq 0 then printf "c=%o : R(s,t,c) identically zero\n", c0; return; end if;
  Cs := Curve(A2b, fs);
  ircomps := IrreducibleComponents(Cs);
  gs := [];
  for K in ircomps do
    try
      Append(~gs, <Degree(K), Genus(Curve(K))>);
    catch e
      Append(~gs, <Degree(K), -1>);
    end try;
  end for;
  printf "c=%o : cubic sing pts %o ; (s,t) slice components (deg,genus) = %o\n",
     c0, sing, gs;
end procedure;

printf "\n=== slices at the singular-fibre c values and at small-height c ===\n";
cand := [Q| -1/2, 0, -1, -2, 1, 2, -1/3, 1/3, -2/3, 2/3, -3/2, 3/2, -1/4, 1/4, -3, 3, 4/17, -3/5, -10/7 ];
for c0 in cand do SliceReport(c0); end for;

printf "LANE1_FIBRATION_DONE\n";
quit;
