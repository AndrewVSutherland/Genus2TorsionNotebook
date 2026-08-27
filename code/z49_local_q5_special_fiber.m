//////////////////////////////////////////////////////////////////////
// Exact special-fiber certificate for the contact-7 Q5-only branch.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
F := GF(3);
P<x> := PolynomialRing(F);
a := F!1;
b := F!1;
h := 1-(F!7/F!2)*x+a*x^2+b*x^3;
f := ExactQuotient(h^2+(x-1)^7,x^2);
g := x^3+2*x+1;
assert f eq x^2*g;
assert Evaluate(g,F!0) eq 1;
E := EllipticCurve([F!0,F!0,F!0,F!2,F!1]);
P7 := E![F!1,F!1,F!1];
assert #E eq 7;
assert Order(P7) eq 7;
print "Q5_ONLY_RESIDUE",<1,1>;
print "h",h,"h(1)",Evaluate(h,F!1);
print "f_factor",Factorization(f);
print "split_node",true,"node_preimages",[<0,1>,<0,2>];
print "normalization",E;
print "E_order",#E,"E_invariants",Invariants(AbelianGroup(E));
print "marked_projection",P7,"marked_projection_order",Order(P7);
print "identity_group_order",(3-1)*#E;
print "conclusion","if 7 does not divide node thickness, D7 is not 7-divisible";
a := F!1;
b := F!0;
hI := 1-(F!7/F!2)*x+a*x^2+b*x^3;
fI := ExactQuotient(hI^2+(x-1)^7,x^2);
assert fI eq x*(x-1)^4;
print "INTERSECTION_RESIDUE",<1,0>;
print "intersection_h",hI,"h(1)",Evaluate(hI,F!1);
print "intersection_f_factor",Factorization(fI);
print "intersection_conclusion","nonordinary four-root cluster; separate blowup required";
quit;
