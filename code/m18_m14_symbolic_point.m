
//////////////////////////////////////////////////////////////////////
//  Symbolic point on the M_1(8,4) family controlling the order-8 chain.
//////////////////////////////////////////////////////////////////////
SetColumns(0);
Q := Rationals();
K<R,w> := RationalFunctionField(Q, 2);
P<x> := PolynomialRing(K);

t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
A := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
B := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);
f := x*A*B;
print "t =", t;
procedure PrintRatFactor(label, q)
    print label, "num", Factorization(Numerator(q));
    print label, "den", Factorization(Denominator(q));
end procedure;
PrintRatFactor("A(-R)", Evaluate(A, -R));
PrintRatFactor("B(-R)", Evaluate(B, -R));
PrintRatFactor("f(-R)", Evaluate(f, -R));
print "sqrt_candidate_squared_ratio_to_f";
Y := -2*R*(R-1)^2*(R^2 - (1/2)*R*w^2 + (1/2)*R - w^2)/(w^2-1);
print "Y =", Y;
PrintRatFactor("f(-R)-Y^2", Evaluate(f,-R) - Y^2);

// Derive a simple square root automatically if possible by factoring f(-R).
quit;
