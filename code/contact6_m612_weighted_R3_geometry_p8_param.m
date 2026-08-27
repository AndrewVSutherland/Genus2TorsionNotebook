//////////////////////////////////////////////////////////////////////
// Clean one-parameter parametrization of the rational P8 component of
// the exact endpoint R3 halving curve.
//
// The normalization conic is y^2=100-6t^2.  Parametrizing through
// (t,y)=(4,-2) gives t=t(u), y=y(u).  The formulas e(t), mu(t) come from
// the normalization of the P8 projection; nu is recovered from its
// monic quadratic gcd using sqrt(discriminant)=h(t)*y(t).
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q:=Rationals();
A<e,mu,nu>:=PolynomialRing(Q,3,"grevlex");
PX<X>:=PolynomialRing(A);
D:=3+5*e;
A1:=(-2*e-3*e^2)*X^2+(6*e+10*e^2)*X+(1-3*e^2);
A2:=2*e*X^2-(1+3*e);
R3:=2-3*X^2;
S:=R3*(mu*X+nu)^2-D*A1*A2;
s4:=Coefficient(S,4); s3:=Coefficient(S,3);
s2:=Coefficient(S,2); s1:=Coefficient(S,1); s0:=Coefficient(S,0);
H1:=8*s4^2*s1-s3*(4*s4*s2-s3^2);
H0:=64*s4^3*s0-(4*s4*s2-s3^2)^2;

K<u>:=FunctionField(Q);
t:=4*(u^2+u-6)/(u^2+6);
y:=2*(u^2-24*u-6)/(u^2+6);
assert y^2 eq 100-6*t^2;

et:=-(Q!25/3)*t^2/(t^4-25*t^2+Q!1250/3);
mt:=(-5*t^7+Q!1750/9*t^5-Q!6250/3*t^3)/
    (t^8-50*t^6+Q!4375/3*t^4-Q!62500/3*t^2+Q!1562500/9);
Bnu:=(-Q!500/3*t^5+Q!175000/27*t^3-Q!625000/9*t)/
    (t^8-50*t^6+Q!4375/3*t^4-Q!62500/3*t^2+Q!1562500/9);
hdisc:=(t^7-Q!500/9*t^5+Q!28750/27*t^3-Q!62500/9*t)/
    (t^8-50*t^6+Q!4375/3*t^4-Q!62500/3*t^2+Q!1562500/9);
nt:=(-Bnu+hdisc*y)/2;

toK:=hom<A -> K | et,mt,nt>;
assert toK(H1) eq 0;
assert toK(H0) eq 0;
assert toK(s4) ne 0;

print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_P8_PARAM";
print "NORMALIZATION_CONIC","y^2=100-6*t^2";
print "T_OF_U",t;
print "Y_OF_U",y;
print "E_OF_T","-(25/3)t^2/(t^4-25t^2+1250/3)";
print "MU_OF_T","(-5t^7+(1750/9)t^5-(6250/3)t^3)/(t^8-50t^6+(4375/3)t^4-(62500/3)t^2+1562500/9)";
print "BNU_OF_T","(-(500/3)t^5+(175000/27)t^3-(625000/9)t)/den_mu";
print "HDISC_OF_T","(t^7-(500/9)t^5+(28750/27)t^3-(62500/9)t)/den_mu";
print "NU_OF_T_Y","(-Bnu+hdisc*y)/2";
print "COMPOSED_DEGREES_E_MU_NU",
      [<Degree(Numerator(q)),Degree(Denominator(q))>:q in [et,mt,nt]];
print "S4_NONZERO_RATIONAL_FUNCTION",toK(s4) ne 0;

u0:=Q!12;
point:=[Q!Evaluate(et,u0),Q!Evaluate(mt,u0),Q!Evaluate(nt,u0)];
assert point eq [Q!-200/409,Q!-36320/167281,Q!38136/167281];
assert Evaluate(H1,point) eq 0 and Evaluate(H0,point) eq 0;
assert Evaluate(s4,point) ne 0;
print "U12_EXACT_OPEN_POINT_e_mu_nu",point;
print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_P8_PARAM_DONE";
quit;
