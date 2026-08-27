
SetColumns(0);
Q := Rationals();
K<R,w> := RationalFunctionField(Q, 2);
P<x> := PolynomialRing(K);
t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
A := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
B := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);
f := x*A*B;
c4 := Coefficient(A*B,4);
print "c4 num", Factorization(Numerator(c4));
print "c4 den", Factorization(Denominator(c4));
D := Discriminant(f);
print "disc num", Factorization(Numerator(D));
print "disc den", Factorization(Denominator(D));
print "disc total degree num", TotalDegree(Numerator(D));
quit;
