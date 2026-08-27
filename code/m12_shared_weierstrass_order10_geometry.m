//////////////////////////////////////////////////////////////////////
// Exact shared-Weierstrass boundary of compact M(12)+5.
//
// The compact odd model is
//
//   y^2 = F(x) = L*(L*H0^2+4*b*(1+x)^2*(w*L-x^2)),
//   L=b+(2*b-1)*x.
//
// Its visible Weierstrass point is W=(omega,0),
// omega=-b/(2*b-1).  Put X=L and scale y by (2*b-1)^2.  Then
//
//   Y^2 = Fhat(X)=(2*b-1)^4 F((X-b)/(2*b-1))=X*G(X).
//
// For P=(R,H(R)), the exact relation
//
//   5*(P-infinity) = W-infinity
//
// is equivalent, on the stated smooth open, to a cubic H with H(0)=0
// and
//
//   Fhat-H^2 = kappa*X*(X-R)^5.
//
// Write H=X*(a*X^2+c*X+d), s=a^2, C=a*c, D=a*d.  Cancelling X and
// comparing coefficients gives the five small sign-quotient equations
// below.  A rational H exists precisely when s is a nonzero square.
//
// Modes:
//   summary  reconstruct and print the exact formulas (default);
//   quotient build the cheap saturated quotient ideal (bounded memory).
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned mode then mode := "summary"; end if;
if mode notin {"summary","quotient"} then
    error "mode must be summary or quotient";
end if;
if not assigned MemGB then MemGB := 3;
elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
R<b,w,Rp,s,C,D> := PolynomialRing(Q,6,"grevlex");
K := FieldOfFractions(R);

PX<X> := PolynomialRing(K);
q := 2*b-1;
x := (X-b)/q;
L := b+q*x;
H0 := x+w*(1+b*x);
F := L*(L*H0^2+4*b*(1+x)^2*(w*L-x^2));
Fhat := PX!(q^4*F);
assert Fhat mod X eq 0;
G := ExactQuotient(Fhat,X);
assert Degree(G) eq 4;
g := [R!Coefficient(G,i) : i in [0..4]];
g0:=g[1]; g1:=g[2]; g2:=g[3]; g3:=g[4]; g4:=g[5];

// Coefficients X^4,...,X^0 of
// G-X*(aX^2+cX+d)^2+s*(X-Rp)^5, after multiplying those
// involving c^2,cd,d^2 by s and putting C=ac,D=ad.
E4 := g4-2*C-5*s*Rp;
E3 := s*g3-C^2-2*s*D+10*s^2*Rp^2;
E2 := s*g2-2*C*D-10*s^2*Rp^3;
E1 := s*g1-D^2+5*s^2*Rp^4;
E0 := g0-s*Rp^5;
eqs := [E4,E3,E2,E1,E0];

// Full reconstruction before quotienting by the simultaneous sign.
RF<bf,wf,Rf,a,cf,df> := PolynomialRing(Q,6,"grevlex");
PF<Z> := PolynomialRing(RF);
toRF := hom<R -> RF | bf,wf,Rf,a^2,a*cf,a*df>;
h2 := a*Z^2+cf*Z+df;
GRF := &+[toRF(g[i+1])*Z^i : i in [0..4]];
reconstruction := GRF-Z*h2^2+a^2*(Z-Rf)^5;
coeffs := [RF!Coefficient(reconstruction,i) : i in [0..5]];
assert coeffs[6] eq 0;
assert coeffs[5] eq toRF(E4);
assert a^2*coeffs[4] eq toRF(E3);
assert a^2*coeffs[3] eq toRF(E2);
assert a^2*coeffs[2] eq toRF(E1);
assert coeffs[1] eq toRF(E0);
print "M12_SHARED_W_GEOMETRY_SELF_TEST_PASS";
print "FHAT",Fhat;
print "G_COEFFICIENTS";
for i in [0..4] do
    print " g",i,"degree",TotalDegree(g[i+1]),"terms",#Terms(g[i+1]),g[i+1];
end for;
print "QUOTIENT_EQUATION_SHAPES",
      [<TotalDegree(f),#Terms(f)> : f in eqs];
print "EQUATIONS",eqs;
print "RECOVERY","a^2=s, c=C/a, d=D/a, H=X*(a*X^2+c*X+d)";
print "POINT","P=(R,H(R)); 5*(P-infinity)=W-infinity";

if mode eq "summary" then quit; end if;

I := ideal<R|eqs>;
cheap := b*w*(b-1)*q*Rp*s*g0;
print "RAW_BASIS_BEGIN";
time BI := Basis(I);
print "RAW_BASIS_LEN",#BI;
print "CHEAP_SATURATION_BEGIN";
time J := Saturation(I,ideal<R|cheap>);
BJ := Basis(J);
print "CHEAP_SATURATION_BASIS_LEN",#BJ;
try
    dim,degs := Dimension(J);
    print "CHEAP_SATURATION_DIMENSION",dim,"COMPONENT_DEGREES",degs;
catch err
    print "DIMENSION_FAILED",err`Object;
end try;
quit;
