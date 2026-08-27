//////////////////////////////////////////////////////////////////////
// Deterministic marked-torsion identities at the explicit HLP Z/60 seed.
//
// The torsion-subgroup generator returned by Magma is seed-dependent, so
// this verifier starts from the published Mumford representative P60 and
// checks its order-5, order-3, order-4, and order-2 multiples exactly.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB:=2; end if;
if Type(MemGB) eq MonStgElt then MemGB:=StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q:=Rationals(); P<x>:=PolynomialRing(Q);
f:=-185*(125*x^2-1728)*(2000*x^4-48525*x^2+294912);
assert Degree(f) eq 6 and Discriminant(f) ne 0;
C:=HyperellipticCurve(f); J:=Jacobian(C); O:=J!0;

P60:=J![x^2+74/15*x+128/25,
          -1081510/9*x-1226624/3];
assert Order(P60) eq 60;

// Deterministic primary multiples of this fixed P60.
q5:=x^2-4608/395;       v5:=-164280/79*x;
q3:=x^2-316/25;         v3:=-5476;
u4:=x^2-1506/125;       v4:=-24642/5;
q0:=x^2-1728/125;
D5:=J![q5,v5]; D3:=J![q3,v3]; D4:=J![u4,v4]; T2:=J![q0,0];

assert D5 eq 12*P60 and Order(D5) eq 5;
assert D3 eq 20*P60 and Order(D3) eq 3;
assert D4 eq 15*P60 and Order(D4) eq 4;
assert T2 eq 30*P60 and Order(T2) eq 2;
assert 2*D4 eq T2;

// Order 5: the norm of A5-B5*y is a fifth power of q5.
// A5 == B5*v5 (mod q5), so A5-B5*y selects D5 rather than -D5.
B5:=q0;
A5:=37/320*x*(52705*x^4-1294080*x^2+7962624);
k5:=341553260289/4096;
assert A5^2-f*B5^2 eq k5*q5^5;
assert A5 mod q5 eq (B5*v5) mod q5;
assert GCD(B5,q5) eq 1;

// Order 3: the norm of H3-y is a cube of q3.
H3:=29600-2775*x^2; k3:=46250000;
assert H3^2-f eq k3*q3^3;
assert H3 mod q3 eq v3;

// Order 4: q0 is the rational two-torsion factor and the norm of
// ell4+y is q0 times a square of u4, certifying 2*D4=T2.
ell4:=-2775*q0; k4:=46250000;
assert f mod q0 eq 0;
assert ell4^2-f eq k4*q0*u4^2;
assert ell4 mod u4 eq (-v4) mod u4;

G,mp:=TorsionSubgroup(J);
assert Invariants(G) eq [60];
assert #AutomorphismGroup(C) eq 4;

print "HLP_Z60_MARKED_IDENTITIES_VERIFY";
print "f",f;
print "P60",P60,"order",Order(P60);
print "D5",D5,"order",Order(D5);
print "D3",D3,"order",Order(D3);
print "D4",D4,"order",Order(D4);
print "T2",T2,"order",Order(T2);
print "ORDER5_IDENTITY",A5^2-f*B5^2 eq k5*q5^5;
print "ORDER3_IDENTITY",H3^2-f eq k3*q3^3;
print "ORDER4_IDENTITY",ell4^2-f eq k4*q0*u4^2;
print "TORSION",Invariants(G),"AUTOMORPHISM_GROUP_ORDER",#AutomorphismGroup(C);
print "HLP_Z60_MARKED_IDENTITIES_VERIFY_DONE";
quit;
