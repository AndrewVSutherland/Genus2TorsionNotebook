//////////////////////////////////////////////////////////////////////
// Exact equivalence of the two order-11/order-22 source presentations.
//////////////////////////////////////////////////////////////////////

Q := Rationals();
Ku<u> := FunctionField(Q);
P<x> := PolynomialRing(Ku);
t := 4*u;

F := x^6 + 2*x^5 + (2*t+3)*x^4 + 2*x^3
     + (t^2+1)*x^2 + 2*t*(1-t)*x + t^2;
G := x^6 - 4*x^5 + 8*(1+u)*x^4 - (10+32*u)*x^3
     + 8*(1+6*u+2*u^2)*x^2
     - 4*(1+6*u+16*u^2)*x + 64*u^2+1;

assert Evaluate(F,x-1) eq G;
print "VERIFIED G_u(x) = F_(4u)(x-1)";
print "coordinate map x_Flynn = x_DS - 1, y_Flynn = y_DS";
print "marked root s^2 maps to 1+s^2";

quit;
