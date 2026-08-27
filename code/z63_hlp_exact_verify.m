//////////////////////////////////////////////////////////////////////
// Exact HLP cyclic [63] example (HLP, equation (4), Section 3.6).
// The marked divisor is reconstructed from the conjugate point pair
// printed in the paper, and the full rational torsion is certified.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemMB then MemMB:=220; end if;
if Type(MemMB) eq MonStgElt then MemMB:=StringToInteger(MemMB); end if;
SetMemoryLimit(MemMB*10^6);

Q:=Rationals(); P<x>:=PolynomialRing(Q);
f:=897*x^6-197570*x^4+79136353*x^2-146398496;
assert Degree(f) eq 6 and Discriminant(f) ne 0;
C:=HyperellipticCurve(f); J:=Jacobian(C);

// The two involutions above the even-sextic symmetry give quotient curves
// with coordinates (X,Y)=(897*x^2,897*y) and
// (X,Y)=(-146398496/x^2,-146398496*y/x^3), respectively.
a:=Coefficient(f,6); b4:=Coefficient(f,4);
c2:=Coefficient(f,2); d:=Coefficient(f,0);
E9:=EllipticCurve([Q!0,b4,Q!0,a*c2,a^2*d]);
E7:=EllipticCurve([Q!0,c2,Q!0,d*b4,d^2*a]);
G9:=TorsionSubgroup(E9); G7:=TorsionSubgroup(E7);
assert Invariants(G9) eq [9] and Invariants(G7) eq [7];

// If s^2=4369, HLP's conjugate affine points are
//   ((-69+s)/2, 4515015-68241*s) and its conjugate.
// Thus their Mumford support and interpolation polynomial are:
u:=x^2+69*x+98;
v:=-136482*x-193614;
assert (v^2-f) mod u eq 0;
D:=J![u,v];
assert 63*D eq J!0 and 21*D ne J!0 and 9*D ne J!0;

G,mp:=TorsionSubgroup(J);
assert Invariants(G) eq [63];
T:=mp(G.1);
assert Order(T) eq 63 and Order(D) eq 63;

// HLP's sharp good-reduction certificate: #J(F_5)=63.
F5:=GF(5); P5<X>:=PolynomialRing(F5);
f5:=P5![F5!Coefficient(f,i): i in [0..Degree(f)]];
C5:=HyperellipticCurve(f5); J5:=Jacobian(C5);
assert #J5 eq 63;

print "HLP_Z63_EXACT";
print "f",f;
print "factor",Factorization(f);
print "marked_D",D,"order",Order(D);
print "TORSION",Invariants(G),"order",#G;
print "generator",T,"order",Order(T);
print "J_F5_order",#J5;
print "elliptic_quotient_torsion",Invariants(G9),Invariants(G7);
print "reduced_model",ReducedMinimalWeierstrassModel(C);
print "GEOMETRICALLY_SPLIT",true;
print "HLP_Z63_EXACT_DONE";
quit;
