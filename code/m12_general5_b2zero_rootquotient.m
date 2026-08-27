//////////////////////////////////////////////////////////////////////
// Birational sign quotient of the b2=0 general order-5 norm cover.
//
// This is a substantially smaller replacement for a direct six-variable
// Groebner elimination in m12_general5_b2zero_geometry.m.  It covers the
// open locus
//
//   b1 * q(-b0/b1) * Disc(q) != 0.
//
// Put c=b0/b1 and z=x+c.  The norm identity at z=0 implies, uniquely,
//
//   d = A(-c)/q(-c)^2,        q(-c)=d^2.
//
// After z=d*Z, write
//
//   q=d^2*(Z^2+e*Z+1),       b1=lambda*d^4.
//
// The involution B -> -B is lambda -> -lambda.  Its quotient coordinate
// is tau=lambda^2.  If G(Z)=F(d*Z-c)=sum(g_i Z^i), comparison of the
// norm identity gives four small equations in (b,w,c,d,e); tau is then
// recovered rationally.  On the further chart g0 != 0, the first equation
// is linear in e, leaving only three equations in (b,w,c,d).
//
// Modes are intentionally static: this file builds and checks the exact
// formulas, and optionally prepares the small quotient/signed ideals.  It
// does not attempt the expensive primary decomposition by default.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned mode then mode := "summary"; end if;
if mode notin {"summary","quotient","signed"} then
    error "mode must be summary, quotient, or signed";
end if;
if not assigned MemGB then MemGB := 8;
elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
R<b,w,c,d,e> := PolynomialRing(Q,5,"grevlex");
P<Z> := PolynomialRing(R);

// Compact M(12) quintic, followed by x=d*Z-c.
x := d*Z-c;
L := b+(2*b-1)*x;
H := x+w*(1+b*x);
Fshift := L*(L*H^2+4*b*(1+x)^2*(w*L-x^2));
G := P!Fshift;
g := [R!Coefficient(G,i) : i in [0..5]];
g0:=g[1]; g1:=g[2]; g2:=g[3]; g3:=g[4]; g4:=g[5]; g5:=g[6];

pcoef := 5*e/2;
qcoef := 5/2+15*e^2/8;
delta := (2-e)^3;
r3 := 5*delta/8;
r4 := 5*delta*(19*e-6)/64;
r5 := delta*(32*e^2-33*e+58)/32;

ell3 := pcoef*g0-g1;
ell4 := qcoef*g0-g2;
ell5 := qcoef*g0-g3;
ell6 := pcoef*g0-g4;
ell7 := g0-g5;

// The five remaining coefficients of
//
// Atilde^2-tau*Z^2*G-(Z^2+eZ+1)^5
//
// after coefficients 0,1,2,8,9,10 have been forced.
E3 := r3; // overwritten below after tau is introduced in a fraction field

// Eliminating tau with E3=r3+tau*ell3 gives these four equations.
C1 := ell3-ell7;
C2 := (32*e^2-33*e+58)*ell3-20*ell5;
C3 := (19*e-6)*ell3-8*ell6;
C4 := 4*(19*e-6)*ell3^2-32*ell4*ell3+5*delta*g0^2;

// Direct symbolic regression of the coefficient formulas.
Rt<t> := PolynomialRing(R);
PRt<Y> := PolynomialRing(Rt);
At := Y^5+pcoef*Y^4+qcoef*Y^3+
      (qcoef+t*g0/2)*Y^2+pcoef*Y+1;
Qt := Y^2+e*Y+1;
Gt := &+[Rt!g[i+1]*Y^i : i in [0..5]];
Nt := At^2-t*Y^2*Gt-Qt^5;
assert &and[Coefficient(Nt,i) eq 0 : i in [0,1,2,8,9,10]];
assert R!Coefficient(Coefficient(Nt,3),0) eq r3;
assert R!Coefficient(Coefficient(Nt,3),1) eq ell3;
assert R!Coefficient(Coefficient(Nt,4),0) eq r4;
assert R!Coefficient(Coefficient(Nt,4),1) eq ell4;
assert R!Coefficient(Coefficient(Nt,4),2) eq g0^2/4;
assert R!Coefficient(Coefficient(Nt,5),0) eq r5;
assert R!Coefficient(Coefficient(Nt,5),1) eq ell5;
assert R!Coefficient(Coefficient(Nt,6),0) eq r4;
assert R!Coefficient(Coefficient(Nt,6),1) eq ell6;
assert R!Coefficient(Coefficient(Nt,7),0) eq r3;
assert R!Coefficient(Coefficient(Nt,7),1) eq ell7;

print "B2ZERO_ROOTQUOTIENT_FORMULAS_PASS";
print "QUOTIENT_5VAR_SHAPES",
      [<TotalDegree(f),#Terms(f)> : f in [C1,C2,C3,C4]];
print "TAU_RECOVERY","tau=-5*(2-e)^3/(8*ell3)";

// Generic g0!=0 chart.  C1 is linear in e.  Put
//
//   h=g0+g1-g5, L0=g0-g5, e=2h/(5g0).
//
// Clearing only powers of g0 gives three equations in four variables.
S<bb,ww,cc,dd> := PolynomialRing(Q,4,"grevlex");
PS<ZZ> := PolynomialRing(S);
xx := dd*ZZ-cc;
LL := bb+(2*bb-1)*xx;
HH := xx+ww*(1+bb*xx);
GG := PS!(LL*(LL*HH^2+4*bb*(1+xx)^2*(ww*LL-xx^2)));
gg := [S!Coefficient(GG,i) : i in [0..5]];
G0:=gg[1]; G1:=gg[2]; G2:=gg[3]; G3:=gg[4]; G4:=gg[5]; G5:=gg[6];
h := G0+G1-G5;
L0 := G0-G5;

K3 := (38*h-30*G0)*L0-40*G0*(h-G4);
K2 := (128*h^2-330*h*G0+1450*G0^2)*L0
      -1250*G0^3-150*G0*h^2+500*G0^2*G3;
K4 := (95*h-75*G0)*L0^2
      -10*(25*G0^2+3*h^2)*L0+100*G0*G2*L0
      +(5*G0-h)^3;

// Check that these are exactly the numerators obtained by substituting
// e=2h/(5G0) in C3,C2,C4, respectively.
KS := FieldOfFractions(S);
erat := KS!(2*h)/KS!(5*G0);
phi := hom<R -> KS | bb,ww,cc,dd,erat>;
assert phi(C3) eq KS!K3/KS!(5*G0);
assert phi(C2) eq KS!K2/KS!(25*G0^2);
assert phi(C4) eq KS!(8*K4)/KS!(25*G0);

print "GENERIC_G0_NONZERO_SHAPES",
      [<TotalDegree(f),#Terms(f)> : f in [K2,K3,K4]];
print "GENERIC_RECOVERY",
      "e=2*h/(5*G0)",
      "tau=-(5*G0-h)^3/(25*G0^3*(G0-G5))";
print "SIGNED_DOUBLE_COVER",
      "25*G0^3*(G0-G5)*lambda^2+(5*G0-h)^3=0";

if mode eq "summary" then
    print "EXCEPTIONAL_G0_ZERO",
          "G0=0, G1=G5!=0, G2=G4, e=(6*G1+8*G4)/(19*G1), plus C2";
    print "REPEATED_Q_WARNING",
          "e=+/-2 is excluded here and must be audited separately";
    quit;
end if;

// Cheap structural saturation only.  The large Disc(F) and Res(q,F)
// should be tested component-by-component after decomposition, instead of
// multiplying them into the first Groebner computation.
cheap := bb*ww*(bb-1)*(2*bb-1)*dd*G0*L0*(h-5*G0)*(h+5*G0);
Iq := ideal<S|K2,K3,K4>;
print "GENERIC_RAW_BASIS_BEGIN";
time Bq := Basis(Iq);
print "GENERIC_RAW_BASIS_LEN",#Bq;
print "GENERIC_CHEAP_SATURATION_BEGIN";
time Jq := Saturation(Iq,ideal<S|cheap>);
print "GENERIC_CHEAP_SATURATION_BASIS_LEN",#Basis(Jq);
try print "GENERIC_CHEAP_SATURATION_DIMENSION",Dimension(Jq);
catch err print "GENERIC_DIMENSION_FAILED",err`Object; end try;

if mode eq "quotient" then quit; end if;

T<bbb,www,ccc,ddd,lambda> := PolynomialRing(Q,5,"grevlex");
psi := hom<S -> T | bbb,www,ccc,ddd>;
signed := 25*psi(G0)^3*psi(L0)*lambda^2+psi(5*G0-h)^3;
Is := ideal<T|[psi(f):f in [K2,K3,K4]] cat [signed]>;
cheapSigned := psi(cheap)*lambda;
print "SIGNED_COVER_CHEAP_SATURATION_BEGIN";
time Js := Saturation(Is,ideal<T|cheapSigned>);
print "SIGNED_COVER_BASIS_LEN",#Basis(Js);
try print "SIGNED_COVER_DIMENSION",Dimension(Js);
catch err print "SIGNED_DIMENSION_FAILED",err`Object; end try;
quit;
