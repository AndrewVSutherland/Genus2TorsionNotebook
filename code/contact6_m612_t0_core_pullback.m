//////////////////////////////////////////////////////////////////////
// Exact pullback of the contact-6 [3,3] cubic-contact core to the
// reduced T0-halving cover.
//
// Variables:
//   s=b+3, omega^2=(a+3)/(b+3), m^2=W,
//   L (with M=L^2), U, nu (constant term nu^2 of q3).
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
R<s,omega,m,L,U,nu> := PolynomialRing(Q,6);

a := s*omega^2-3;
b := s-3;
W := m^2;
M := L^2;

// f=x*(c4g*x^4+c3g*x^3+c2g*x^2+c1g*x+c0g).
c5 := 2*s;
c4 := b^2+2*a-15;
c3 := 2*a*b+22;
c2 := a^2+2*b-15;
c1 := 2*(a+3);

N := c2-omega*c4+omega*W;
R8 := 8*s*(c3-4*s*omega)-(c4-W)^2;
Fhalf := R8^2-256*s^2*W*N;

B3 := c5*M+3*U;
Delta3 := 4*c4*M+12*(U^2+nu^2)-B3^2;
F3 := B3*Delta3+16*nu^3-8*c3*M-8*U^3-48*U*nu^2;
F2 := Delta3^2+64*B3*nu^3-64*c2*M
      -192*(U^2*nu^2+nu^4);
F1 := Delta3*nu^3-4*c1*M-12*U*nu^4;

print "CONTACT6_M612_T0_CORE_PULLBACK";
print "variables (s,omega,m,L,U,nu)";
for rec in [<"Fhalf",Fhalf>,<"F3",F3>,<"F2",F2>,<"F1",F1>] do
    print rec[1],"total_degree",TotalDegree(rec[2]),
          "terms",#Terms(rec[2]);
end for;
print "Fhalf_degree_m",Degree(Fhalf,3);
print "Fhalf_irreducible",#Factorization(Fhalf) eq 1;
print "boundary_open_product",
      s*m*L*nu*(omega^2-1)*(U^2-4*nu^2)*(a+b+2);

quit;
