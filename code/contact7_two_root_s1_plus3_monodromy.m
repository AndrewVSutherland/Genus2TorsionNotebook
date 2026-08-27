//////////////////////////////////////////////////////////////////////
// Exact singular-fiber and transvection audit for the cancelled s=1
// branch of the contact-7 two-root family.
//
// The nonzero J[3] support modulo sign has degree 40.  At every simple
// nodal fiber, Picard--Lefschetz gives a symplectic transvection with
// cycle shape 1^13 3^9 and ramification contribution 18.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q:=Rationals(); K<t>:=RationalFunctionField(Q); P<x>:=PolynomialRing(K);
a:=K!(Q!35/8);
b:=-(8*t^4+24*t^3+48*t^2+45*t+15)/(8*(t+1)^3);
f:=x^5+(b^2-7)*x^4+(2*a*b+21)*x^3
   +(a^2-7*b-35)*x^2+(2*b-7*a+35)*x+(2*a-Q!35/4);
assert Coefficient(f,0) eq 0;
w:=1-t^2;
assert Evaluate(f,w) eq 0;

disc:=Discriminant(f);
num:=Numerator(disc); den:=Denominator(disc);
print "CONTACT7_TWO_ROOT_S1_PLUS3_MONODROMY";
print "b",b;
print "f_factor",Factorization(f);
print "disc_num_degree",Degree(num),"disc_den_degree",Degree(den);
print "disc_num_factors",[<Degree(q[1]),q[2],q[1]>:q in Factorization(num)];
print "disc_den_factors",[<Degree(q[1]),q[2],q[1]>:q in Factorization(den)];

function Canonical(v)
    nv:=[(-z) mod 3:z in v]; return v lt nv select v else nv;
end function;
S:={};
for z in CartesianPower({0,1,2},4) do
    v:=[Integers()!u:u in z];
    if &or[u ne 0:u in v] then Include(~S,Canonical(v)); end if;
end for;
seen:={}; lengths:=[];
for v in S do
    if v in seen then continue; end if;
    u:=v; ell:=0;
    repeat
        Include(~seen,u); ell+:=1;
        // transvection about e1 in (e1,e2,f1,f2) coordinates
        u:=Canonical([(u[1]-u[3]) mod 3,u[2],u[3],u[4]]);
    until u in seen;
    Append(~lengths,ell);
end for;
print "support_degree",#S;
print "transvection_cycles",[<e,#[z:z in lengths|z eq e]>:e in Seqset(lengths)];
print "ramification_per_simple_node",#S-#lengths;

// Every squarefree discriminant factor of multiplicity one contributes
// its degree in simple nodal fibers, provided the residual cubic collision
// is ordinary.  Record the rigorous RH lower bound from those factors.
simple_nodes:=&+[Degree(q[1]):q in Factorization(num)|q[2] eq 1];
ram:=#S-#lengths;
print "simple_discriminant_roots",simple_nodes;
print "connected_cover_genus_lower_bound",1+Ceiling((-2*#S+simple_nodes*ram)/2);
quit;
