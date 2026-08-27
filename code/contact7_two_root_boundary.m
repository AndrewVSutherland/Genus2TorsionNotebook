//////////////////////////////////////////////////////////////////////
// Symbolic boundary factorization for the contact-7 two-root surface.
// Signs are absorbed into signed parameters S,T, so h(1-S^2)=S^7
// and h(1-T^2)=T^7.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
R<S,T> := PolynomialRing(Q,2);
K := FieldOfFractions(R);
PK<X> := PolynomialRing(K);

r := 1-S^2;
w := 1-T^2;
Us := S^7-1+(Q!7/2)*r;
Ut := T^7-1+(Q!7/2)*w;
den := r^2*w^2*(w-r);
a := (Us*w^3-Ut*r^3)/den;
b := (r^2*Ut-w^2*Us)/den;
h := 1-(K!7/2)*X+a*X^2+b*X^3;
f := ExactQuotient(h^2+(X-1)^7,X^2);

assert h^2-X^2*f eq -(X-1)^7;
assert Evaluate(h,K!r) eq K!S^7;
assert Evaluate(h,K!w) eq K!T^7;
assert Evaluate(f,K!r) eq 0 and Evaluate(f,K!w) eq 0;
g3 := ExactQuotient(f,(X-r)*(X-w));
assert Degree(g3) eq 3;

H := S^2*T^2+2*S^2*T+S^2+2*S*T^2+2*S*T+(Q!1/2)*S+T^2+(Q!1/2)*T;
Rst := S^4*T^2+2*S^4*T+S^4+2*S^3*T^3+7*S^3*T^2+8*S^3*T+3*S^3
       +6*S^2*T^3+13*S^2*T^2+8*S^2*T+S^2+6*S*T^3+7*S*T^2
       +2*S*T+2*T^3+T^2;
assert Evaluate(h,1) eq -S^2*T^2*H/((S+1)^2*(T+1)^2*(S+T));

print "parameter_denominator_factorization",Factorization(R!den);
print "reduced_a_denominator_factorization",Factorization(Denominator(a));
print "reduced_b_denominator_factorization",Factorization(Denominator(b));
print "h1_numerator_factorization",Factorization(Numerator(Evaluate(h,1)));
print "h1_denominator_factorization",Factorization(Denominator(Evaluate(h,1)));
print "fprime_r_numerator_factorization",Factorization(Numerator(Evaluate(Derivative(f),r)));
print "fprime_w_numerator_factorization",Factorization(Numerator(Evaluate(Derivative(f),w)));
print "cubic_discriminant_numerator_factorization",Factorization(Numerator(Discriminant(g3)));
print "cubic_discriminant_denominator_factorization",Factorization(Denominator(Discriminant(g3)));
fullfac := Factorization(Numerator(Discriminant(f)));
print "full_discriminant_numerator_factorization",fullfac;
fulldenfac := Factorization(Denominator(Discriminant(f)));
print "full_discriminant_denominator_factorization",fulldenfac;
numprod := R!1;
for q in fullfac do numprod *:= q[1]^q[2]; end for;
denprod := R!1;
for q in fulldenfac do denprod *:= q[1]^q[2]; end for;
numunit := ExactQuotient(Numerator(Discriminant(f)),numprod);
denunit := ExactQuotient(Denominator(Discriminant(f)),denprod);
print "full_discriminant_scalar",Q!numunit/(Q!denunit);
Kfac := [ q[1] : q in fullfac | Degree(q[1],1) eq 8 and Degree(q[1],2) eq 8 ][1];
print "residual_K_bidegree",Degree(Kfac,1),Degree(Kfac,2);
print "residual_K_symmetric",Evaluate(Kfac,[T,S]) eq Kfac;
expected_disc := 8*S^14*T^14*(S-T)^2*H^7*Rst^2*Evaluate(Rst,[T,S])^2*Kfac
                 /((S+1)^26*(T+1)^26*(S+T)^13);
assert Discriminant(f) eq expected_disc;
print "normalized_discriminant_scalar_with_Rst",8;

quit;
