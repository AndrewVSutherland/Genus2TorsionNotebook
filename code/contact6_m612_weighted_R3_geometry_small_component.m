//////////////////////////////////////////////////////////////////////
// Normalize the smaller endpoint R3 projection component.
//
// Its (e,mu)-projection is a genus-zero degree-8 plane curve and the
// generic gcd in nu is quadratic.  Parametrize the plane curve, pull the
// quadratic gcd to Q(t), and reduce its discriminant to a hyperelliptic
// squarefree model for the normalization of the full space component.
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

B<ee,mm>:=PolynomialRing(Q,2,"grevlex");
toB:=hom<A -> B | ee,mm,B!0>;
fac:=Factorization(B!toB(Resultant(H1,H0,nu)));
P:=[f[1]:f in fac|Degree(f[1],mm) eq 4 and GCD(f[1],B!toB(s4)) eq 1][1];
Cp:=ProjectiveClosure(Curve(AffineSpace(B),P));
assert Genus(Cp) eq 0;

print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_SMALL_COMPONENT";
print "PROJECTION_DEGREE",Degree(Cp),"GENUS",Genus(Cp);
Ccon,mc:=Conic(Cp);
print "NORMALIZATION_CONIC",Ccon;
print "CONIC_MAP_DOMAIN",Domain(mc),"CONIC_MAP_CODOMAIN",Codomain(mc);
print "CONIC_LOCALLY_SOLVABLE",IsLocallySolvable(Ccon);
haspoint,pcon:=HasRationalPoint(Ccon);
print "CONIC_HAS_RATIONAL_POINT",haspoint;
if haspoint then print "CONIC_POINT",pcon; end if;
assert haspoint;

psi:=Inverse(mc);
eqs:=DefiningEquations(psi);
print "CONIC_TO_PLANE_COORDINATE_DEGREES",[Degree(g):g in eqs];

Kt<t>:=FunctionField(Q);
// The returned conic is Y^2=X*Z, parametrized by (X:Y:Z)=(t^2:t:1).
et:=Kt!Evaluate(eqs[1],[t^2,t,Kt!1])/
    (Kt!Evaluate(eqs[3],[t^2,t,Kt!1]));
mt:=Kt!Evaluate(eqs[2],[t^2,t,Kt!1])/
    (Kt!Evaluate(eqs[3],[t^2,t,Kt!1]));
print "PARAM_E",et;
print "PARAM_MU",mt;

KtN<n>:=PolynomialRing(Kt);
function PullToKtN(h)
    ans:=KtN!0;
    for j in [0..Degree(h,nu)] do
        c:=Coefficient(h,nu,j);
        cb:=B!toB(c);
        ans+:=Evaluate(cb,[et,mt])*n^j;
    end for;
    return ans;
end function;
g:=GCD(PullToKtN(H1),PullToKtN(H0));
g:=g/LeadingCoefficient(g);
assert Degree(g) eq 2;
disc:=Discriminant(g);
issquare,sqrtdisc:=IsSquare(disc);
print "NU_GCD",g;
print "NU_DISCRIMINANT_IS_SQUARE",issquare;
print "NU_DISCRIMINANT_NUM_FACTORIZATION",Factorization(Numerator(disc));
print "NU_DISCRIMINANT_DEN_FACTORIZATION",Factorization(Denominator(disc));
disc_square_factor:=t*(t^2-Q!50/3)*
    (t^4-Q!350/9*t^2+Q!1250/3)/
    (t^4-25*t^2+Q!1250/3)^2;
disc_constant:=disc/((t^2-Q!50/3)*disc_square_factor^2);
assert disc_constant in Q;
print "NU_DISCRIMINANT_CONSTANT_SQUARECLASS",disc_constant;
print "NU_DISCRIMINANT_SQUARE_FACTOR",disc_square_factor;
P2c<Xc,Yc,Zc>:=ProjectiveSpace(Q,2);
dc:=Q!disc_constant;
Cdisc:=Conic(P2c,Yc^2-dc*(Xc^2-Q!50/3*Zc^2));
print "FULL_COMPONENT_NORMALIZATION_CONIC",Cdisc;
print "FULL_COMPONENT_CONIC_LOCALLY_SOLVABLE",IsLocallySolvable(Cdisc);
full_haspoint,full_point:=HasRationalPoint(Cdisc);
print "FULL_COMPONENT_CONIC_HAS_RATIONAL_POINT",full_haspoint;
if full_haspoint then
    print "FULL_COMPONENT_CONIC_POINT",full_point;
    conic_param:=Parametrization(Cdisc,full_point);
    print "FULL_COMPONENT_CONIC_PARAMETRIZATION",
          DefiningEquations(conic_param);
end if;
// The displayed conic point has t=4 and normalized square root -2.
// Recover one exact open point on the original H1=H0 curve.
t0:=Q!4; y0:=Q!-2;
e0:=Q!Evaluate(et,Kt!t0);
mu0:=Q!Evaluate(mt,Kt!t0);
b0:=Q!Evaluate(Coefficient(g,1),Kt!t0);
sqrt_disc0:=Q!Evaluate(disc_square_factor,Kt!t0)*y0;
nu0:=(-b0+sqrt_disc0)/2;
assert Evaluate(H1,[e0,mu0,nu0]) eq 0;
assert Evaluate(H0,[e0,mu0,nu0]) eq 0;
assert Evaluate(s4,[e0,mu0,nu0]) ne 0;
print "EXACT_OPEN_RATIONAL_POINT_e_mu_nu",[e0,mu0,nu0];

// y^2=N/D is birational to Y^2=N*D.  Remove even factor powers.
Qtt:=PolynomialRing(Q); tt:=Qtt.1;
ND:=Qtt!(Numerator(disc)*Denominator(disc));
sf:=Qtt!1;
for fe in Factorization(ND) do
    if IsOdd(fe[2]) then sf*:=fe[1]; end if;
end for;
if LeadingCoefficient(sf) lt 0 then sf:=-sf; end if;
print "DISCRIMINANT_SQUAREFREE_MODEL",sf;
print "DISCRIMINANT_SQUAREFREE_DEGREE",Degree(sf);
if Degree(sf) ge 3 then
    Ch:=HyperellipticCurve(sf);
    print "FULL_COMPONENT_NORMALIZATION_GENUS",Genus(Ch);
else
    print "FULL_COMPONENT_NORMALIZATION_GENUS",0;
end if;
print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_SMALL_COMPONENT_DONE";
quit;
