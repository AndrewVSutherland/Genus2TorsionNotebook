//////////////////////////////////////////////////////////////////////
// Independent extraction audit for the r=4 [2,6,6] plane images.
//
// Run under an external timeout:
//   magma -b code/audit_contact6_m36_266_plane_extract.m
//
// This rebuilds only the structural saturation (known to equal the full
// open saturation in the recorded exact run), changes independently to a
// lexicographic ring, and prints the two primitive plane factors in (U,v).
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
Z := Integers();
R<b,M,U,v> := PolynomialRing(Q, 4, "grevlex");
Px<x> := PolynomialRing(R);

r := Q!4;
a := 3 - (b+3)*r - 2/r;
h := 1 + a*x + b*x^2 + x^3;
f := h^2 - (x-1)^6;
c := [Coefficient(f,i) : i in [0..5]];
B := c[6]*M + 3*U;
Delta := 4*c[5]*M + 12*(U^2+v^2) - B^2;
F1 := Delta*v^3 - 4*c[2]*M - 12*U*v^4;
F2 := Delta^2 + 64*B*v^3 - 64*c[3]*M
      - 192*(U^2*v^2+v^4);
F3 := B*Delta + 16*v^3 - 8*c[4]*M - 8*U^3 - 48*U*v^2;
structural := (b+3)*M*v*(U^2-4*v^2)*(a+b+2);

Iraw := ideal<R | F1,F2,F3>;
I := Saturation(Iraw, ideal<R | structural>);

L<bb,MM,UU,vv> := PolynomialRing(Q,4,"lex");
phi := hom<R -> L | bb,MM,UU,vv>;
IL := ideal<L | [phi(g) : g in Basis(I)]>;
GL := GroebnerBasis(IL);
plane := [g : g in GL | Degree(g,1) eq 0 and Degree(g,2) eq 0];
require #plane eq 1: "expected a unique (U,v) elimination polynomial";
fac := Factorization(plane[1]);
require #fac eq 2 and {TotalDegree(t[1]) : t in fac} eq {11,21}:
    "unexpected plane factorization";

print "AUDIT_EXTRACT_START";
print "saturated_dimension_degree", Dimension(I), Degree(Homogenization(I));
print "lex_basis_length", #GL;
print "plane_total_degree_terms", TotalDegree(plane[1]), #Terms(plane[1]);
print "plane_factor_summary",
      [<TotalDegree(t[1]),Degree(t[1],3),Degree(t[1],4),#Terms(t[1]),t[2]>
       : t in fac];

for t in fac do
    g := t[1];
    d := TotalDegree(g);
    print "FACTOR_BEGIN", d;
    print "irreducible_over_Q", IsIrreducible(g), "multiplicity", t[2];
    print "arithmetic_plane_genus", (d-1)*(d-2) div 2;
    print "value_mod13_at_U4_v6",
          Evaluate(ChangeRing(g,GF(13)), <GF(13)!0,GF(13)!0,GF(13)!4,GF(13)!6>);
    print "value_mod19_at_U9_v3",
          Evaluate(ChangeRing(g,GF(19)), <GF(19)!0,GF(19)!0,GF(19)!9,GF(19)!3>);
    print "factor_polynomial", g;
    print "FACTOR_END", d;
end for;

print "AUDIT_EXTRACT_DONE";
