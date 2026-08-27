//////////////////////////////////////////////////////////////////////
//  E=0 slice of the exact cubic-contact cover for the cyclic-order-20
//  contact-5 family, in the low-degree s-chart
//
//    s=(t-1)/(t+1),  x_old=(1-s)*x,
//    g_s=(1+(1+s)x+2*s*x^2)^2-4*(1-s)*x^5.
//
//  A rational 3-class has
//
//    H^2-q^3=M*g_s,   q=x^2+U*x+V,   M=L^2.
//
//  On E=H(0)=0 the constant equation is -V^3=M.  Hence write
//
//    V=-R, M=R^3, U=-(2/3)*R*(s+1).
//
//  The formula for U is the x-coefficient equation.  The square condition
//  is exactly R=r^2.  After the recursive x^5,x^4,x^3 equations, only the
//  x^2 equation and Erec=0 remain, so this is a zero-dimensional slice in
//  (s,R).  We saturate every genuine open boundary before extracting the
//  rational points.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(4*10^9);

Q := Rationals();
P<s,R> := PolynomialRing(Q,2,"grevlex");
PX<x> := PolynomialRing(P);

h := 1+(1+s)*x+2*s*x^2;
g := h^2-4*(1-s)*x^5;
c := [P!Coefficient(g,i) : i in [0..5]];

assert c[1] eq 1;
assert c[2] eq 2*(s+1);
assert c[3] eq s^2+6*s+1;
assert c[4] eq 4*s*(s+1);
assert c[5] eq 4*s^2;
assert c[6] eq 4*(s-1);

V := -R;
M := R^3;
U := -(Q!2/3)*R*(s+1);
A := (M*c[6]+3*U)/2;
B := (M*c[5]+3*(U^2+V)-A^2)/2;
Erec := (M*c[4]+U^3+6*U*V-2*A*B)/2;

G2raw := B^2-3*(U^2*V+V^2)-M*c[3];
G1raw := -3*U*V^2-M*c[2];

// Positive algebraic controls: all matched coefficients are reconstructed
// directly, and G1 vanishes identically after the E=0 substitution.
q := x^2+U*x+V;
H := x^3+A*x^2+B*x;
identity := H^2-q^3-M*g;
assert Coefficient(identity,5) eq 0;
assert Coefficient(identity,4) eq 0;
assert Coefficient(identity,3) eq -2*Erec;
assert Coefficient(identity,2) eq G2raw;
assert Coefficient(identity,1) eq G1raw;
assert Coefficient(identity,0) eq 0;
assert G1raw eq 0;
print "CONTACT_EZERO_SELF_TEST_PASS";

function PrimitiveP(f)
    if f eq 0 then return f; end if;
    den := LCM([Denominator(a) : a in Coefficients(f)]);
    ints := [Integers()!(den*a) : a in Coefficients(f)];
    cont := GCD(ints);
    if cont eq 0 then cont := 1; end if;
    return P!((Q!den/Q!Abs(cont))*f);
end function;

FE := PrimitiveP(Erec);
F2 := PrimitiveP(G2raw);
discq := PrimitiveP(U^2-4*V);
resqg := PrimitiveP(P!Resultant(q,g));
Ds := 8*s^3-59*s^2-18*s+197;
boundary := R*(s-1)*(s-2)*Ds*discq*resqg;

print "EZERO equation_degrees",[TotalDegree(FE),TotalDegree(F2)],
      "terms",[#Terms(FE),#Terms(F2)];
print "EZERO boundary_degree",TotalDegree(boundary),
      "discq",Factorization(discq);

common2 := GCD(FE,F2);
print "EZERO_COMMON_FACTOR degree",TotalDegree(common2),
      "terms",#Terms(common2),"factorization",Factorization(common2);
FEred := ExactQuotient(FE,common2);
F2red := ExactQuotient(F2,common2);

Iraw := ideal<P | FE,F2>;
print "EZERO_RAW dimension",Dimension(Iraw),
      "basis_degrees",[TotalDegree(f) : f in Basis(Iraw)];
time Isat := Saturation(Iraw,ideal<P | boundary>);
print "EZERO_SAT dimension",Dimension(Isat),
      "basis_degrees",[TotalDegree(f) : f in Basis(Isat)],
      "basis_terms",[#Terms(f) : f in Basis(Isat)];

// Exact elimination in R, with s retained as the univariate base variable.
QS<ss> := PolynomialRing(Q);
PR<rr> := PolynomialRing(QS);
mp := hom<P -> PR | PR!ss,rr>;
fe := mp(FEred);
f2 := mp(F2red);
time elim_s := Resultant(fe,f2);
print "EZERO_RESULTANT_S degree",Degree(elim_s),
      "terms",#Coefficients(elim_s);
fac_s := Factorization(elim_s);
print "EZERO_RESULTANT_FACTORS",#fac_s;
for z in fac_s do
    print " S_FACTOR multiplicity",z[2],"degree",Degree(z[1]),
          "terms",#Coefficients(z[1]),z[1];
end for;

QR<zR> := PolynomialRing(Q);
rat_s := Roots(elim_s);
print "EZERO_RATIONAL_S_ROOTS",rat_s;

candidates := [];
for sr in rat_s do
    s0 := sr[1];
    fer := QR![Evaluate(Coefficient(fe,i),s0) : i in [0..Degree(fe)]];
    f2r := QR![Evaluate(Coefficient(f2,i),s0) : i in [0..Degree(f2)]];
    common := GCD(fer,f2r);
    for zr in Roots(common) do
        R0 := zr[1];
        if Evaluate(FE,[s0,R0]) ne 0 or Evaluate(F2,[s0,R0]) ne 0 then
            continue;
        end if;
        open := Evaluate(boundary,[s0,R0]) ne 0;
        is_square,r0 := IsSquare(R0);
        print "EZERO_RATIONAL_POINT s",s0,"R",R0,
              "open",open,"R_square",is_square;
        if open and is_square then
            Append(~candidates,<s0,R0,r0>);
        end if;
    end for;
end for;

print "EZERO_OPEN_SQUARE_CANDIDATES",#candidates,candidates;

// Any square candidate is rebuilt from scratch and exact-checked.  This
// block is normally empty but makes a positive result immediately usable.
QX<X> := PolynomialRing(Q);
for tup in candidates do
    s0,R0,r0 := Explode(tup);
    U0 := -(Q!2/3)*R0*(s0+1);
    V0 := -R0;
    M0 := R0^3;
    g0 := (1+(1+s0)*X+2*s0*X^2)^2-4*(1-s0)*X^5;
    cc := [Coefficient(g0,i) : i in [0..5]];
    A0 := (M0*cc[6]+3*U0)/2;
    B0 := (M0*cc[5]+3*(U0^2+V0)-A0^2)/2;
    H0 := X^3+A0*X^2+B0*X;
    q0 := X^2+U0*X+V0;
    assert H0^2-q0^3 eq M0*g0;
    C := HyperellipticCurve(g0);
    J := Jacobian(C);
    TG,phi := TorsionSubgroup(J);
    inv := Invariants(TG);
    t0 := (1+s0)/(1-s0);
    print "EZERO_EXACT_CANDIDATE s",s0,"t",t0,"R",R0,
          "r",r0,"g",g0,"q",q0,"H",H0,
          "torsion",inv,"cyclic60",inv eq [60];
end for;

print "EZERO_DONE";
quit;
