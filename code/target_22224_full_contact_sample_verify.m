// Independent Magma check of one row emitted by
// target_22224_full_contact_audit.py.
if not assigned p then p:=17; end if;
if not assigned vals then vals:=[1,2,3,6]; end if;
if not assigned wit then wit:=[7,11,6]; end if;

F:=GF(p); P<x>:=PolynomialRing(F);
a,b,c,d:=Explode([F!z:z in vals]);
L,U,v:=Explode([F!z:z in wit]);
M:=L^2;
A0:=a^2+b^2+c^2+d^2;
A1:=a^2*b^2+a^2*c^2+a^2*d^2+b^2*c^2+b^2*d^2+c^2*d^2;
A2:=a^2*b^2*c^2+a^2*b^2*d^2+a^2*c^2*d^2+b^2*c^2*d^2;
A3:=a^2*b^2*c^2*d^2;
PP:=4*M*A0+12*(U^2+v^2)-(M+3*U)^2;
m:=1/L;
h:=m*x^3+(M+3*U)/(2*L)*x^2+PP/(8*L)*x+v^3/L;
q:=x^2+U*x+v^2;
f:=x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2);
assert f eq x^5+A0*x^4+A1*x^3+A2*x^2+A3*x;
assert h^2-f eq m^2*q^3;
C:=HyperellipticCurve(f);
n:=Integers()!Evaluate(LPolynomial(C),1);
print "CONTACT_SAMPLE",p,vals,wit,"q",q,"h",h,"Jorder",n,"divisible3",n mod 3 eq 0;
assert n mod 3 eq 0;
