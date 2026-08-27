//////////////////////////////////////////////////////////////////////
// Modular connectivity certificate for the degree-40 support cover
// (J[3]-{0})/{+1,-1} on the cancelled s=1 contact-7 branch.
//
// The normalized contact equations are
//
//   w*(1+r1*x+r2*x^2+r3*x^3)^2 = f + k*(x^2+U*x+V)^3.
//
// Saturating w removes only a nongeneric affine chart.  A degree-40
// irreducible resolvent over F_p(t) certifies connectivity over Q(t).
//////////////////////////////////////////////////////////////////////

SetColumns(0); SetSeed(1);
if not assigned p then p:=101;
elif Type(p) eq MonStgElt then p:=StringToInteger(p); end if;
Fp:=GF(p); K<t>:=FunctionField(Fp);
a:=K!(Fp!35/Fp!8);
b:=-(8*t^4+24*t^3+48*t^2+45*t+15)/(8*(t+1)^3);
f0:=K!0;
f1:=2*b-7*a+35;
f2:=a^2-7*b-35;
f3:=2*a*b+21;
f4:=b^2-7;
f5:=K!1; f6:=K!0;

P<iv,r1,r2,r3,w,k,V,U>:=PolynomialRing(K,8,"grevlex");
C0:=w-f0-k*V^3;
C1:=2*w*r1-f1-3*k*U*V^2;
C2:=w*(r1^2+2*r2)-f2-3*k*(U^2*V+V^2);
C3:=2*w*(r3+r1*r2)-f3-k*(U^3+6*U*V);
C4:=w*(r2^2+2*r1*r3)-f4-3*k*(U^2+V);
C5:=2*w*r2*r3-f5-3*k*U;
C6:=w*r3^2-f6-k;
I:=ideal<P|C0,C1,C2,C3,C4,C5,C6,iv*w-1>;
print "CONTACT7_TWO_ROOT_S1_PLUS3_MODULAR_RESOLVENT p",p;
print "equation_shapes",[<TotalDegree(e),#Terms(e)>:e in [C0,C1,C2,C3,C4,C5,C6,iv*w-1]];
time G0:=GroebnerBasis(I);
dim,pars:=Dimension(I); deg:=Dimension(quo<P|I>);
print "grevlex_basis_size",#G0,"dimension",dim,"parameters",pars,"degree",deg;
time IL:=ChangeOrder(I,"lex");
time G:=GroebnerBasis(IL);
print "lex_basis_size",#G;
univ:=[g:g in G|&and[Degree(g,j) eq 0:j in [1..7]] and Degree(g,8) gt 0];
assert #univ gt 0;
R:=UnivariatePolynomial(univ[#univ]); R/:=LeadingCoefficient(R);
fac:=Factorization(R);
print "resolvent_degree",Degree(R);
print "factor_degrees",[<Degree(q[1]),q[2]>:q in fac];
print "connected_certificate",Degree(R) eq 40 and #fac eq 1 and fac[1][2] eq 1;
quit;
