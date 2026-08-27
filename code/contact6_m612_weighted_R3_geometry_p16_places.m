//////////////////////////////////////////////////////////////////////
// Resolve rational singular plane points of P16 into places.  A degree-1
// place would give a Q-point on the normalization and an elliptic model.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q:=Rationals();
A<e,mu,nu>:=PolynomialRing(Q,3,"grevlex");
PX<X>:=PolynomialRing(A);
D:=3+5*e;
A1:=(-2*e-3*e^2)*X^2+(6*e+10*e^2)*X+(1-3*e^2);
A2:=2*e*X^2-(1+3*e);
R3:=2-3*X^2;
S:=R3*(mu*X+nu)^2-D*A1*A2;
s4:=Coefficient(S,4); s3:=Coefficient(S,3);
s2:=Coefficient(S,2); s1:=Coefficient(S,1); s0:=Coefficient(S,0);
H1:=8*s4^2*s1-s3*(4*s4*s2-s3^2);
H0:=64*s4^3*s0-(4*s4*s2-s3^2)^2;
B<ee,mm>:=PolynomialRing(Q,2,"grevlex");
toB:=hom<A -> B | ee,mm,B!0>;
fac:=Factorization(B!toB(Resultant(H1,H0,nu)));
P16:=[f[1]:f in fac|Degree(f[1],ee) eq 16 and Degree(f[1],mm) eq 8][1];
Cp:=ProjectiveClosure(Curve(AffineSpace(B),P16));

print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_P16_PLACES";
print "NORMALIZED_GENUS",Genus(Cp);
sing:=SingularPoints(Cp);
print "RATIONAL_SINGULARS",sing;
found:=false;
for q in sing do
    print "PLACE_RESOLUTION_START",q;
    pls:=Places(q);
    print " PLACE_COUNT",#pls,"DEGREES",[Degree(pl):pl in pls];
    for pl in pls do
        if Degree(pl) eq 1 and not found then
            found:=true;
            print " DEGREE_ONE_PLACE",pl;
            E,mp:=EllipticCurve(Cp,pl);
            Emin:=MinimalModel(E);
            print " ELLIPTIC_MODEL",Emin;
            print " ELLIPTIC_AINVARIANTS",aInvariants(Emin);
            print " ELLIPTIC_RANK_BOUNDS",RankBounds(Emin);
            print " ELLIPTIC_TORSION",Invariants(TorsionSubgroup(Emin));
        end if;
    end for;
end for;
print "DEGREE_ONE_PLACE_FOUND",found;
print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_P16_PLACES_DONE";
quit;
