//////////////////////////////////////////////////////////////////////
// Exact audit of the removable u=1 boundary on the signed two-root
// contact-7 surface.
//
// This is the compatible root-at-zero limit.  The incompatible u=-1
// choice would require h(0)=-1 although h(0)=1.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q:=Rationals();
R<t>:=PolynomialRing(Q);
K:=FieldOfFractions(R);
P<x>:=PolynomialRing(K);

V:=8*t^4+24*t^3+48*t^2+45*t+15;
a:=K!(Q!35/8);
b:=-K!V/(8*(t+1)^3);
r:=K!0;
w:=1-t^2;
h:=1-(K!7/2)*x+a*x^2+b*x^3;
f:=ExactQuotient(h^2+(x-1)^7,x^2);

assert Coefficient(f,0) eq 0;
assert Evaluate(f,w) eq 0;
assert Evaluate(h,w) eq t^7;
assert h^2-x^2*f eq -(x-1)^7;
assert Evaluate(h,1) eq -t^2*(8*t^2+9*t+3)/(8*(t+1)^3);
g3:=ExactQuotient(f,x*(x-w));
assert Degree(g3) eq 3;

// Exact projective limit t=infinity.  Divide y by t, hence f by t^2.
function LimitAtInfinityOfScaled(c,k)
    n:=Numerator(c); d:=Denominator(c);
    e:=Degree(n)-Degree(d)-k;
    if e lt 0 then return Q!0; end if;
    assert e eq 0;
    return LeadingCoefficient(n)/LeadingCoefficient(d);
end function;
finf:=&+[LimitAtInfinityOfScaled(Coefficient(f,i),2)*x^i:i in [0..5]];
hinf:=&+[LimitAtInfinityOfScaled(Coefficient(h,i),1)*x^i:i in [0..3]];
assert hinf eq -x^3;
assert finf eq x^4;

print "CONTACT7_TWO_ROOT_PLUS3_U1_BOUNDARY";
print "a",a;
print "b",b;
print "r",r,"w",w;
print "h1 factorization",Factorization(Numerator(Evaluate(h,1))),
      Factorization(Denominator(Evaluate(h,1)));
print "discriminant numerator",Factorization(Numerator(Discriminant(f)));
print "discriminant denominator",Factorization(Denominator(Discriminant(f)));
print "residual cubic discriminant numerator",
      Factorization(Numerator(Discriminant(g3)));
print "residual cubic discriminant denominator",
      Factorization(Denominator(Discriminant(g3)));
print "t=infinity scaled_h_limit",hinf;
print "t=infinity scaled_f_limit",finf,"singular",Discriminant(finf) eq 0;
print "classification";
print " u=1 removable and generically smooth";
print " u=-1 incompatible: h(0)=1 cannot equal (-1)^7";
print " u=0 singular/contact collision";
print " u=v coincident-root singularity";
print " u=-v incompatible opposite values at the same root";
print " t=infinity on u=1 degenerates to y^2=x^4 after y/t scaling";

quit;
