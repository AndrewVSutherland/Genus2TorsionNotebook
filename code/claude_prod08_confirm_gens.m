// Confirm the Mumford generators on the family (integral scaled model):
// y^2 = F(x;c^2,1): D3 = (x^2-4c^2, 2cx) order 3, D5 = (x^2-4c^2, 4c^2) order 5, D3+D5 order 15.
SetColumns(0);
Qx<x> := PolynomialRing(Rationals());
q2 := -1/4;
for c in [1, 2, 3, 1/2, 7, 27/8] do
  Q := q2*x^2 + x + c^2;
  Qp := Q - x;
  h := Qp - x*Q;
  F := h^2 + 4*Q^2*Qp;
  den := LCM([Denominator(co) : co in Coefficients(F)]);
  Fi := F*den^2;   // y scales by den
  J := Jacobian(HyperellipticCurve(Fi));
  D3 := elt<J | x^2 - 4*c^2, Qx!(2*c*den*x)>;
  D5 := elt<J | x^2 - 4*c^2, Qx!(4*c^2*den)>;
  printf "c=%o: Order(D3)=%o Order(D5)=%o Order(D3+D5)=%o\n", c, Order(D3), Order(D5), Order(D3+D5);
end for;
quit;
