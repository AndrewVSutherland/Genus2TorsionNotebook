//////////////////////////////////////////////////////////////////////
//  Algebraic extra-2 conditions on the contact-6 core [6,6] surface.
//
//  The contact-6 polynomial factors identically as
//
//      f = x * Qp * Qm,
//      Qp = (b+3)*x^2 + (a-3)*x + 2,
//      Qm = 2*x^2 + (b-3)*x + (a+3).
//
//  The core [6,6] construction supplies two independent 3-parts and the
//  two 2-parts coming from the three rational factors x,Qp,Qm.  To enlarge
//  from [6,6] to exact [2,6,6], force exactly one of Qp,Qm to split over Q.
//
//  On the generic cubic-contact branch F1=0 gives a=N/D in variables
//  (b,M,U,v), M=L^2.  The two square conditions are:
//
//      yp^2 = (N - 3D)^2 - 8*(b+3)*D^2,
//      ym^2 = D*((b-3)^2*D - 8*N - 24*D).
//
//  Equivalently, use a rational root r of one quadratic:
//
//      Qp(r)=0: a =  3 - (b+3)*r - 2/r,
//      Qm(r)=0: a = -3 - 2*r^2 - (b-3)*r.
//
//  These are the two branches of h6(r)=eps*(r-1)^3.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
Z := Integers();

function PrimitivePolynomial(f)
    if f eq 0 then
        return f;
    end if;
    coeffs := Coefficients(f);
    den := LCM([Denominator(c) : c in coeffs]);
    g := Parent(f)!(den*f);
    nums := [Z!c : c in Coefficients(g)];
    cont := GCD([Abs(n) : n in nums | n ne 0]);
    if cont gt 1 then
        g := Parent(f)!(g/cont);
    end if;
    return g;
end function;

function FactorSummary(f)
    if f eq 0 then
        return "0";
    end if;
    return Sprint([<TotalDegree(fe[1]), #Terms(fe[1]), fe[2]> : fe in Factorization(f)]);
end function;

S<b,M,U,v> := PolynomialRing(Q, 4);
K := FieldOfFractions(S);
PA<a> := PolynomialRing(K);

c1 := 2*a + 6;
c2 := a^2 + 2*b - 15;
c3 := 2*a*b + 22;
c4 := 2*a + b^2 - 15;
c5 := 2*b + 6;

B := c5*M + 3*U;
Delta := 4*c4*M + 12*(U^2 + v^2) - B^2;
F3 := B*Delta + 16*v^3 - 8*c3*M - 8*U^3 - 48*U*v^2;
F2 := Delta^2 + 64*B*v^3 - 64*c2*M - 192*(U^2*v^2 + v^4);
F1 := Delta*v^3 - 4*c1*M - 12*U*v^4;

D := Coefficient(F1, 1);
N := -Coefficient(F1, 0);
assert D eq 8*M*(v^3 - 1);

function SubstituteANumerator(P)
    d := Degree(P);
    value := K!0;
    for i in [0..d] do
        value +:= Coefficient(P, i)*N^i*D^(d-i);
    end for;
    return PrimitivePolynomial(S!Numerator(value));
end function;

G2 := SubstituteANumerator(F2);
G3 := SubstituteANumerator(F3);
common := GCD(G2, G3);
G2red := ExactQuotient(G2, common);
G3red := ExactQuotient(G3, common);

SplitP := PrimitivePolynomial(S!((N - 3*D)^2 - 8*(b+3)*D^2));
SplitM := PrimitivePolynomial(S!(D*((b-3)^2*D - 8*N - 24*D)));

print "Contact-6 core extra-2 algebraic conditions";
print "Variables: b,M,U,v with M=L^2";
print "D =", D;
print "a = N/D";
print "N degree", TotalDegree(PrimitivePolynomial(S!N)),
      "degrees", [Degree(PrimitivePolynomial(S!N), i) : i in [1..4]],
      "terms", #Terms(PrimitivePolynomial(S!N));
print "G2red degree", TotalDegree(G2red), "terms", #Terms(G2red),
      "factor", FactorSummary(G2red);
print "G3red degree", TotalDegree(G3red), "terms", #Terms(G3red),
      "factor", FactorSummary(G3red);
print "";
print "Qp = (b+3)*x^2 + (a-3)*x + 2";
print "Qm = 2*x^2 + (b-3)*x + (a+3)";
print "Extra 2 for exact [2,6,6]: exactly one of Qp,Qm splits.";
print "";
print "SplitP square-cover equation: yp^2 = SplitP";
print "SplitP degree", TotalDegree(SplitP),
      "degrees", [Degree(SplitP, i) : i in [1..4]],
      "terms", #Terms(SplitP), "factor", FactorSummary(SplitP);
print "SplitP =", SplitP;
print "";
print "SplitM square-cover equation: ym^2 = SplitM";
print "SplitM degree", TotalDegree(SplitM),
      "degrees", [Degree(SplitM, i) : i in [1..4]],
      "terms", #Terms(SplitM), "factor", FactorSummary(SplitM);
print "SplitM =", SplitM;
print "";
print "Root parameter branches:";
print "  Qp(r)=0: a = 3 - (b+3)*r - 2/r";
print "  Qm(r)=0: a = -3 - 2*r^2 - (b-3)*r";
print "  unified: a = (eps*(r-1)^3 - 1 - r^3 - b*r^2)/r";

quit;
