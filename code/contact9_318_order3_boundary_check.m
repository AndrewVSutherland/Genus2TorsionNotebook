// Independent Magma check of the normalized order-3 contact fiber at t=4 mod 7.

F := GF(7);
P<u,v,h0,h1,h2,h3> := PolynomialRing(F, 6);
R<x> := PolynomialRing(P);
f := x^5 + 6*x^4 + x^3 + 2*x + 5;
q := x^2 + u*x + v;
H := h0 + h1*x + h2*x^2 + h3*x^3;
contact := H^2 - f - h3^2*q^3;
equations := [Coefficient(contact, i) : i in [0..5]];
I := ideal<P | equations>;
points := Variety(I);
J := Matrix(P, 6, 6, [Derivative(equation, variable) :
    equation in equations, variable in [u,v,h0,h1,h2,h3]]);

print "contact9 [3,18] order-3 cover check";
print "field", F;
print "f", f, "factorization", Factorization(f);
print "ideal_dimension", Dimension(I);
print "solution_count", #points;
for point in Sort(Setseq(Seqset(points))) do
    values := [point[index] : index in [1..6]];
    evaluated_jacobian := Matrix(F, 6, 6,
        [Evaluate(entry, values) : entry in Eltseq(J)]);
    q_at_node := Evaluate(Evaluate(q, F!6), values);
    print "solution", values,
          "q_at_node", q_at_node,
          "fiber_rank", Rank(evaluated_jacobian),
          "fiber_det", Determinant(evaluated_jacobian);
end for;

Q := Rationals();
Qx<X> := PolynomialRing(Q);
function IntegralPolynomialModel(poly)
    denominator_lcm := 1;
    for index in [0..Degree(poly)] do
        denominator_lcm := LCM(
            denominator_lcm, Denominator(Coefficient(poly, index))
        );
    end for;
    return Parent(poly)!(denominator_lcm^2*poly), denominator_lcm;
end function;
f4 := X^5 + 1313043849/100000000*X^4 - 411903/16000*X^3
      + 5305461/160000*X^2 - 787923/40000*X + 179469/40000;
f4_integral, scale4 := IntegralPolynomialModel(f4);
C4 := HyperellipticCurve(f4_integral);
J4 := Jacobian(C4);
T4, torsion_map := TorsionSubgroup(J4);
print "exact_t4";
print "f4", f4;
print "integral_model", f4_integral, "y_scale", scale4;
print "discriminant", Discriminant(f4);
print "torsion_invariants", Invariants(T4);
