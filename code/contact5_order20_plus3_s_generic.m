//////////////////////////////////////////////////////////////////////
//  Low-degree square cover for exact cyclic Z/60 on the full
//  contact-5/order-20 family.
//
//    s=(t-1)/(t+1), x_old=(1-s)*x,
//    g_s=(1+(1+s)x+2*s*x^2)^2-4*(1-s)*x^5.
//
//  For H^2-q^3=L^2*g_s, q=x^2+U*x+V, the constant equation is
//  E^2-L^2=V^3.  On V*L != 0 it is parameterized by
//
//    L=V*(w^2-V)/(2*w), E=V*(w^2+V)/(2*w).
//
//  The exceptional E=0 chart is globally closed in
//  contact5_order20_plus3_s_ezero.m.  Here E != 0, so the x^1 equation
//  solves linearly for B.  The x^5 equation solves for A.  The remaining
//  x^4,x^3,x^2 equations are the three generators below.  This is smaller
//  than recursively forming B and E and then clearing their denominators.
//
//  Modes: summary (default), saturate, decompose.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned mode then mode := "summary"; end if;
if mode notin {"summary","saturate","decompose"} then
    error "mode must be summary, saturate, or decompose";
end if;
if not assigned MemGB then MemGB := 8;
elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
S<s,U,V,w> := PolynomialRing(Q,4,"grevlex");
K := FieldOfFractions(S);
KX<x> := PolynomialRing(K);

h := 1+(1+s)*x+2*s*x^2;
g := h^2-4*(1-s)*x^5;
c := [K!Coefficient(g,i) : i in [0..5]];
assert c eq [K!1,2*(s+1),s^2+6*s+1,4*s*(s+1),4*s^2,4*(s-1)];

Delta := w^2-V;
Sigma := w^2+V;
L := K!(V*Delta/(2*w));
E := K!(V*Sigma/(2*w));
M := L^2;
assert E^2-L^2 eq V^3;

// Triangular reconstruction from x^5 and x^1.
A := (M*c[6]+3*U)/2;
B := (3*U*V^2+M*c[2])/(2*E);
q := x^2+U*x+V;
H := x^3+A*x^2+B*x+E;
identity := H^2-q^3-M*g;

assert Coefficient(identity,6) eq 0;
assert Coefficient(identity,5) eq 0;
assert Coefficient(identity,1) eq 0;
assert Coefficient(identity,0) eq 0;

function PrimitiveS(f)
    if f eq 0 then return f; end if;
    den := LCM([Denominator(a) : a in Coefficients(f)]);
    ints := [Integers()!(den*a) : a in Coefficients(f)];
    cont := GCD(ints);
    if cont eq 0 then cont := 1; end if;
    return S!((Q!den/Q!Abs(cont))*f);
end function;

F4 := PrimitiveS(S!Numerator(Coefficient(identity,4)));
F3 := PrimitiveS(S!Numerator(Coefficient(identity,3)));
F2pre := PrimitiveS(S!Numerator(Coefficient(identity,2)));
// V is nonzero on this chart and is a visible common factor in F2pre.
vpow := 0;
F2 := F2pre;
while IsDivisibleBy(F2,V) do
    F2 := ExactQuotient(F2,V);
    vpow +:= 1;
end while;

// Independent direct formulas are a regression control for the triangular
// construction.  They use A=AN/(2w^2), B=BN/(2wSigma), E=EN/(2w).
AN := (s-1)*V^2*Delta^2+3*U*w^2;
BN := V*(6*U*w^2+(s+1)*Delta^2);
EN := V*Sigma;
F4direct := PrimitiveS(AN^2*Sigma + 4*w^3*BN
             -12*w^4*Sigma*(U^2+V)
             -4*s^2*V^2*Delta^2*w^2*Sigma);
F3direct := PrimitiveS(2*w^2*Sigma*EN + AN*BN
             -2*w^3*Sigma*(U^3+6*U*V)
             -2*s*(s+1)*V^2*Delta^2*w*Sigma);
F2direct := PrimitiveS(w*BN^2 + 2*Sigma^2*AN*EN
             -12*w^3*Sigma^2*(U^2*V+V^2)
             -w*Sigma^2*V^2*Delta^2*(s^2+6*s+1));
while IsDivisibleBy(F2direct,V) do
    F2direct := ExactQuotient(F2direct,V);
end while;
assert F4 eq F4direct or F4 eq -F4direct;
assert F3 eq F3direct or F3 eq -F3direct;
assert F2 eq F2direct or F2 eq -F2direct;
print "CONTACT_S_GENERIC_SELF_TEST_PASS";

discq := U^2-4*V;
resqg := PrimitiveS(S!Resultant(KX!q,KX!g));
Ds := 8*s^3-59*s^2-18*s+197;

// Sigma=0 is E=0, already closed exactly in the companion file.  The
// repeated-q branch is genuine 3-torsion geometry, but on a smooth odd
// quintic it factors g as 2+3; together with the fixed root x=1/(1-s),
// this forces rational 2-rank >=2.  It is therefore excluded only because
// this file targets exact cyclic [60], not because it is invalid torsion.
boundary := V*w*Delta*Sigma*(s-1)*(s-2)*Ds*discq*resqg;
Iraw := ideal<S | F4,F3,F2>;

print "CONTACT_S_GENERIC mode",mode;
print "equation_degrees",[TotalDegree(f) : f in [F4,F3,F2]],
      "terms",[#Terms(f) : f in [F4,F3,F2]],"removed_V_power_F2",vpow;
print "boundary_degree",TotalDegree(boundary),"terms",#Terms(boundary);
print "raw_dimension",Dimension(Iraw),
      "basis_size",#Basis(Iraw),
      "basis_degrees",[TotalDegree(f) : f in Basis(Iraw)];

if mode eq "summary" then
    // Audit s=infinity, corresponding to the singular original fiber t=-1.
    // With S0=1/s and division by s^2 (a square), the limiting polynomial is
    // [S0+(S0+1)x+2x^2]^2+4S0(1-S0)x^5; at S0=0 it is a square.
    QX<X> := PolynomialRing(Q);
    ginf := (X+2*X^2)^2;
    assert ginf eq X^2*(1+2*X)^2;
    print "S_INFINITY_TMINUS1_BOUNDARY",ginf,"singular",Discriminant(ginf) eq 0;
    quit;
end if;

print "SATURATING_GENERIC_COVER";
time Isat := Saturation(Iraw,ideal<S | boundary>);
print "sat_dimension",Dimension(Isat),
      "basis_size",#Basis(Isat),
      "basis_degrees",[TotalDegree(f) : f in Basis(Isat)],
      "basis_terms",[#Terms(f) : f in Basis(Isat)];

if mode eq "decompose" then
    print "DECOMPOSING_GENERIC_COVER";
    time comps := PrimaryDecomposition(Isat);
    print "components",#comps;
    for i in [1..#comps] do
        C := comps[i];
        print "COMPONENT",i,"dimension",Dimension(C),"prime",IsPrime(C),
              "basis_size",#Basis(C),
              "basis_degrees",[TotalDegree(f) : f in Basis(C)],
              "basis_terms",[#Terms(f) : f in Basis(C)];
    end for;
end if;

print "CONTACT_S_GENERIC_DONE";
quit;
