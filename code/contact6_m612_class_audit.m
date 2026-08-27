//////////////////////////////////////////////////////////////////////
// Rational 2-torsion classes on the simple contact-6 [6,6] hit,
// and a sanity check on the earlier M(2,12) chart.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
P<x> := PolynomialRing(Q);

a := 133/39;
b := -7/13;
h := 1+a*x+b*x^2+x^3;
B := h-(x-1)^3;
C3 := h+(x-1)^3;
C := ExactQuotient(C3,x);
f := x*B*C;
J := Jacobian(HyperellipticCurve(f));
D := J![x-1,Evaluate(h,1)];
T0 := J![x,0];
TB := J![B/LeadingCoefficient(B),0];
TC := J![C/LeadingCoefficient(C),0];

print "SIMPLE_HIT_FACTORS", x, B, C;
print "T0", T0, "order", Order(T0);
print "TB", TB, "order", Order(TB);
print "TC", TC, "order", Order(TC);
print "T0+TB=TC", T0+TB eq TC;
print "3D=T0", 3*D eq T0;
print "3D=TB", 3*D eq TB;
print "3D=TC", 3*D eq TC;
df := LCM([Denominator(Coefficient(f,i)) : i in [0..Degree(f)]]);
JI := Jacobian(HyperellipticCurve(df^2*f));
DI := JI![x-1,df*Evaluate(h,1)];
T0I := JI![x,0];
TBI := JI![B/LeadingCoefficient(B),0];
TCI := JI![C/LeadingCoefficient(C),0];
print "D_HALVES", IsDivisibleBy(DI,2);
print "T0_HALVES", IsDivisibleBy(T0I,2);
print "TB_HALVES", IsDivisibleBy(TBI,2);
print "TC_HALVES", IsDivisibleBy(TCI,2);

z := -8/3;
r := 5;
m := (1-z^2)/(4*(r+1));
T := m*x^2-x+r;
hh := (x-r)*(T+1);
W := hh^2+4*m*x^2*T*(T+1);
print "M212_SAMPLE", "z",z,"r",r,"m",m;
print "M212_FACTOR", Factorization(W);
dW := LCM([Denominator(Coefficient(W,i)) : i in [0..Degree(W)]]);
WI := dW^2*W;
GG, mp := TorsionSubgroup(Jacobian(HyperellipticCurve(WI)));
print "M212_TORSION", Invariants(GG);

quit;
