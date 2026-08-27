//////////////////////////////////////////////////////////////////////
//  Natural V=1 discovery slice of the exact cubic-contact cover on
//
//    g_s=(1+(1+s)x+2s*x^2)^2-4(1-s)x^5.
//
//  The full coefficient equations form a zero-dimensional quotient ideal
//  in (s,M,U).  After exact open saturation and rational-point recovery we
//  retain only M=L^2 and exact-check every surviving curve.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(6*10^9);
Q := Rationals();
P<s,M,U> := PolynomialRing(Q,3,"grevlex");
PX<x> := PolynomialRing(P);

V := P!1;
h := 1+(1+s)*x+2*s*x^2;
g := h^2-4*(1-s)*x^5;
c := [P!Coefficient(g,i) : i in [0..5]];
assert c eq [P!1,2*(s+1),s^2+6*s+1,4*s*(s+1),4*s^2,4*(s-1)];
A := (M*c[6]+3*U)/2;
B := (M*c[5]+3*(U^2+V)-A^2)/2;
E := (M*c[4]+U^3+6*U*V-2*A*B)/2;
q := x^2+U*x+V;
H := x^3+A*x^2+B*x+E;
identity := H^2-q^3-M*g;
G2raw := B^2+2*A*E-3*(U^2*V+V^2)-M*c[3];
G1raw := 2*B*E-3*U*V^2-M*c[2];
G0raw := E^2-V^3-M;
assert &and[Coefficient(identity,i) eq 0 : i in [3..6]];
assert [Coefficient(identity,i) : i in [0..2]] eq [G0raw,G1raw,G2raw];
print "CONTACT_V1_SELF_TEST_PASS";

function PrimitiveP(f)
    if f eq 0 then return f; end if;
    den := LCM([Denominator(a) : a in Coefficients(f)]);
    ints := [Integers()!(den*a) : a in Coefficients(f)];
    cont := GCD(ints); if cont eq 0 then cont := 1; end if;
    return P!((Q!den/Q!Abs(cont))*f);
end function;

G2 := PrimitiveP(G2raw); G1 := PrimitiveP(G1raw); G0 := PrimitiveP(G0raw);
discq := U^2-4;
resqg := PrimitiveP(P!Resultant(q,g));
Ds := 8*s^3-59*s^2-18*s+197;
boundary := M*(s-1)*(s-2)*Ds*discq*resqg;
Iraw := ideal<P | G2,G1,G0>;
print "V1 equation_degrees",[TotalDegree(f) : f in [G2,G1,G0]],
      "terms",[#Terms(f) : f in [G2,G1,G0]],
      "raw_dimension",Dimension(Iraw);
time Isat := Saturation(Iraw,ideal<P | boundary>);
print "V1_SAT dimension",Dimension(Isat),"basis_size",#Basis(Isat),
      "basis_degrees",[TotalDegree(f) : f in Basis(Isat)],
      "basis_terms",[#Terms(f) : f in Basis(Isat)];

RL<UL,ML,SL> := PolynomialRing(Q,3,"lex");
mp := hom<P -> RL | SL,ML,UL>;
Jraw := ideal<RL | mp(G2),mp(G1),mp(G0)>;
time J := Saturation(Jraw,ideal<RL | mp(boundary)>);
time gb := GroebnerBasis(J);
elim := [f : f in gb | Degree(f,1) eq 0 and Degree(f,2) eq 0];
print "V1_LEX dimension",Dimension(J),"basis_size",#gb,
      "variable_degrees",[[Degree(f,i) : i in [1..3]] : f in gb],
      "elim_generators",#elim;

QS<zS> := PolynomialRing(Q);
toQS := hom<RL -> QS | QS!0,QS!0,zS>;
spols := [QS!toQS(f) : f in elim | f ne 0];
if #spols eq 0 then print "V1_NO_S_ELIMINATION"; quit; end if;
s_core := spols[1];
for i in [2..#spols] do s_core := GCD(s_core,spols[i]); end for;
s_core := s_core/LeadingCoefficient(s_core);
fac := Factorization(s_core);
print "V1_S_ELIM degree",Degree(s_core),"terms",#Coefficients(s_core),
      "factors",#fac;
for fe in fac do
    print " S_FACTOR multiplicity",fe[2],"degree",Degree(fe[1]),
          "terms",#Coefficients(fe[1]),fe[1];
end for;
print "V1_RATIONAL_S_ROOTS",Roots(s_core);

pts := Variety(J);
print "V1_RATIONAL_POINTS",#pts,pts;
QX<X> := PolynomialRing(Q);
square_candidates := 0; hits := 0;
for pp in pts do
    u0 := Q!pp[1]; m0 := Q!pp[2]; s0 := Q!pp[3];
    if Evaluate(mp(boundary),[u0,m0,s0]) eq 0 then continue; end if;
    issq,l0 := IsSquare(m0);
    print "V1_OPEN_POINT s",s0,"M",m0,"U",u0,"M_square",issq;
    if not issq then continue; end if;
    square_candidates +:= 1;
    g0 := (1+(1+s0)*X+2*s0*X^2)^2-4*(1-s0)*X^5;
    cc := [Coefficient(g0,i) : i in [0..5]];
    a0 := (m0*cc[6]+3*u0)/2;
    b0 := (m0*cc[5]+3*(u0^2+1)-a0^2)/2;
    e0 := (m0*cc[4]+u0^3+6*u0-2*a0*b0)/2;
    q0 := X^2+u0*X+1;
    H0 := X^3+a0*X^2+b0*X+e0;
    assert H0^2-q0^3 eq m0*g0;
    C0 := HyperellipticCurve(g0); J0 := Jacobian(C0);
    TG,phi := TorsionSubgroup(J0); inv := Invariants(TG);
    t0 := (1+s0)/(1-s0);
    print "V1_EXACT_CANDIDATE s",s0,"t",t0,"M",m0,"L",l0,
          "U",u0,"g",g0,"q",q0,"H",H0,"torsion",inv,
          "cyclic60",inv eq [60];
    if inv eq [60] then hits +:= 1; end if;
end for;
print "V1_DONE square_candidates",square_candidates,"cyclic60_hits",hits;
quit;
