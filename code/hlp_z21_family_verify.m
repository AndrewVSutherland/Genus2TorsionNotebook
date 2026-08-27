//////////////////////////////////////////////////////////////////////
// One-parameter cyclic [21] family obtained from HLP Lemma 11 and
// Proposition 4.  The generic construction is checked over Q(t), and
// t=2 is verified exactly over Q.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemMB then MemMB:=200; end if;
if Type(MemMB) eq MonStgElt then MemMB:=StringToInteger(MemMB); end if;
SetMemoryLimit(MemMB*10^6);

Q:=Rationals(); K<t>:=FunctionField(Q);

// Tate normal form on X_1(7):
// y^2 + (1-c)xy - b y = x^3 - b x^2,  P7=(0,0),
// b=t^3-t^2 and c=t^2-t.
b:=t^3-t^2; c:=t^2-t;
a1:=1-c; a2:=-b; a3:=-b;
E7:=EllipticCurve([a1,a2,a3,K!0,K!0]);
P7:=E7![0,0,1];
assert 7*P7 eq E7!0 and P7 ne E7!0;

// Short model y^2=x^3+A*x+B of the same elliptic curve.
b2:=a1^2+4*a2;
b4:=a1*a3;
b6:=a3^2;
c4:=b2^2-24*b4;
c6:=-b2^3+36*b2*b4-216*b6;
A:=-c4/48; B:=-c6/864;
assert A ne 0 and B ne 0 and 4*A^3+27*B^2 ne 0;
E7short:=EllipticCurve([K!0,K!0,K!0,A,B]);
P7short:=E7short![b2/12,a3/2,1];
assert 7*P7short eq E7short!0 and P7short ne E7short!0;

// HLP Lemma 11: this curve has a rational point P3=(0,q) of
// order 3 and its 2-torsion is Galois-isomorphic to E7[2].
q:=B^2/A^3;
E3:=EllipticCurve([K!0,K!1,K!0,-2*q,q^2]);
P3:=E3![0,q,1];
assert 3*P3 eq E3!0 and P3 ne E3!0;

// Applying HLP Proposition 4 to the root map
// beta=-alpha^2/A+(B/A^2)*alpha-1 gives, after a rational
// square rescaling, this compact even sextic.
R<x>:=PolynomialRing(K);
h:=-A*B^3*x^6 + (A^6+3*A^3*B^2)*x^4
   -3*A^5*B*x^2 + A^7;
assert Degree(h) eq 6;
assert Discriminant(h) eq 64*A^52*B^3*(4*A^3+27*B^2)^2;

// The two natural degree-2 quotients are exactly E3 and E7short.
// These identities also certify the nonsquare twist in the sextic.
F:=FieldOfFractions(R); xx:=F!x; hh:=F!h;
X3:=q-(B^3/A^5)*xx^2;
X7:=A/xx^2-B/A;
assert (B^6/A^16)*hh eq X3^3+(X3-q)^2;
assert hh/(A^4*xx^6) eq X7^3+A*X7+B;

print "HLP_Z21_FAMILY_GENERIC";
print "A",Factorization(Numerator(A)),"/",Factorization(Denominator(A));
print "B_degrees",Degree(Numerator(B)),Degree(Denominator(B));
print "P7_order",7,"P3_order",3;
print "family","y^2=A^7-3*A^5*B*x^2+(A^6+3*A^3*B^2)*x^4-A*B^3*x^6";
print "degree2_quotients_verified",true;

// Exact small specialization.  Denominator 2^28 is a square, so the
// displayed integral sextic is Q-isomorphic to h at t=2.
Pz<X>:=PolynomialRing(Q);
h2:=43*(292754944*X^6-93515024*X^4
        +6810251592*X^2-6321363049);
assert Discriminant(h2) ne 0;
C:=HyperellipticCurve(h2); J:=Jacobian(C);
G,mp:=TorsionSubgroup(J);
print "SPECIAL_T",2;
print "SPECIAL_H",h2;
Cmin,minmap:=ReducedMinimalWeierstrassModel(C);
print "SPECIAL_REDUCED_MODEL",Cmin;
print "SPECIAL_TORSION",Invariants(G),"order",#G;
assert Invariants(G) eq [21];
P21:=mp(G.1);
assert Order(P21) eq 21;
print "SPECIAL_P21",P21,"order",Order(P21);
print "GEOMETRICALLY_SPLIT",true;
print "HLP_Z21_FAMILY_VERIFY_DONE";
quit;
