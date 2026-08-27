// Unsigned order-3 q-cover for the rational-root contact-9 family.
//
// Put G=H/h3 and z=h3^(-2).  Then
//     G^2-z*f=(x^2+u*x+v)^3,
// with G=x^3+a*x^2+b*x+c.  The coefficients of x^5,x^4,x^3
// solve a,b,c triangularly, leaving three equations in u,v,z.

SetColumns(0);

Q := Rationals();
Kt<t> := FunctionField(Q);
Kx<x> := PolynomialRing(Kt);

r := 1-t^2;
h_at_r := 1-9/2*r+63/8*r^2-105/16*r^3;
contact_a := (t^9-h_at_r)/r^4;
marked_h := 1-9/2*x+63/8*x^2-105/16*x^3+contact_a*x^4;
quotient, remainder := Quotrem(marked_h^2+(x-1)^9, x^4);
assert remainder eq 0;
f := quotient;
assert Degree(f) eq 5 and Coefficient(f, 5) eq 1;

P<u,v,z> := PolynomialRing(Kt, 3, "grevlex");
f0 := P!Coefficient(f, 0);
f1 := P!Coefficient(f, 1);
f2 := P!Coefficient(f, 2);
f3 := P!Coefficient(f, 3);
f4 := P!Coefficient(f, 4);
f5 := P!Coefficient(f, 5);

a := (3*u+z*f5)/2;
b := (z*f4+3*u^2+3*v-a^2)/2;
c := (z*f3+u^3+6*u*v-2*a*b)/2;

R2 := b^2+2*a*c-z*f2-3*(u^2*v+v^2);
R1 := 2*b*c-z*f1-3*u*v^2;
R0 := c^2-z*f0-v^3;
equations := [R2, R1, R0];
Iraw := ideal<P | equations>;

print "contact9 [3,18] unsigned q-cover elimination";
print "base_field Q(t)";
print "f", f;
print "triangular_a", a;
print "triangular_b", b;
print "triangular_c", c;
for index in [1..3] do
    print "equation", index,
          "total_degree", TotalDegree(equations[index]),
          "term_count", #Terms(equations[index]);
    print equations[index];
end for;
print "raw_ideal_dimension", Dimension(Iraw);
PX<X> := PolynomialRing(P);
f_over_P := &+[P!Coefficient(f, index)*X^index : index in [0..5]];
q_over_P := X^2+u*X+v;
q_f_resultant := Resultant(q_over_P, f_over_P);
boundary := z*(u^2-4*v)*q_f_resultant;
I := Saturation(Iraw, ideal<P | boundary>);
print "saturation_boundary", boundary;
print "ideal_dimension", Dimension(I);
print "quotient_degree", Degree(I);

GrobnerBasis(~I);
L<U,V,Z> := PolynomialRing(Kt, 3, "lex");
Ilex, order_map := ChangeOrder(I, L);
basis := Basis(Ilex);
print "groebner_basis_count", #basis;
for index in [1..#basis] do
    print "basis", index,
          "total_degree", TotalDegree(basis[index]),
          "term_count", #Terms(basis[index]);
    print basis[index];
end for;

univariate := [polynomial : polynomial in basis |
    Degree(polynomial, 1) eq 0 and Degree(polynomial, 2) eq 0];
print "z_univariate_count", #univariate;
for polynomial in univariate do
    print "z_eliminant_degree", Degree(polynomial, 3);
    factors := Factorization(polynomial);
    print "z_factor_degrees", [<Degree(pair[1], 3), pair[2]> : pair in factors];
    for pair in factors do
        print "z_factor", Degree(pair[1], 3), pair[2], pair[1];
    end for;
end for;
