//////////////////////////////////////////////////////////////////////
// Symbolic resultant/norm identities behind the cheap discriminant
// covers for halving the distinguished Richelot dual-kernel classes.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q:=Rationals(); A<a,b>:=PolynomialRing(Q,2);
P<x>:=PolynomialRing(A);
A0:=a+3; S:=b+3; Delta:=A0*S-4;
R1:=(b^2-2*a-3)*x^2+(2*a*b+6*a+6*b+10)*x+(a^2-2*b-3);
R2:=2*x^2-A0;
R3:=2-S*x^2;
DB:=(a-3)^2-8*(b+3);
DC:=(b-3)^2-8*(a+3);

assert R1 eq (S^2-6*S+12-2*A0)*x^2+2*Delta*x
                +(A0^2-6*A0+12-2*S);
assert Resultant(R1,R2) eq Delta^2*DC;
assert Resultant(R1,R3) eq Delta^2*DB;
assert Resultant(R2,R3) eq Delta^2;

// If 2H=[Ri,0], the exact quartic-square identity is
//
//   Ri*L^2-Delta*Rj*Rk = k*u^2.
//
// Reducing modulo Ri and taking norms makes the product of the two
// corresponding resultants a square.  The three necessary square covers
// are therefore
//
//   R1: DB*DC=square,  R2: DC=square,  R3: DB=square.

print "CONTACT6_M612_DUAL_DISCRIMINANT_COVERS";
print "Res(R1,R2)",Factorization(Resultant(R1,R2));
print "Res(R1,R3)",Factorization(Resultant(R1,R3));
print "Res(R2,R3)",Factorization(Resultant(R2,R3));
print "NECESSARY","R1:DB*DC square","R2:DC square","R3:DB square";
quit;
