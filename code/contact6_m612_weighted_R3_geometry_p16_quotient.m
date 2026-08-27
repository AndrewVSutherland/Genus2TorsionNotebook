//////////////////////////////////////////////////////////////////////
// Quotient of the P16 genus-one component by mu -> -mu.
//
// Put z=mu^2.  The quotient again has normalized genus one, but unlike
// P16 it has visible rational plane points.  Test which are smooth and,
// when possible, use one to construct an elliptic model by Riemann-Roch.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned PointBound then PointBound:=500; end if;
if Type(PointBound) eq MonStgElt then PointBound:=StringToInteger(PointBound); end if;
if not assigned DoRank then DoRank:=true; end if;
if Type(DoRank) eq MonStgElt then
    DoRank:=DoRank in {"true","True","1","yes","Yes"};
end if;

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

T<et,z>:=PolynomialRing(Q,2,"grevlex");
coeffToT:=hom<B -> T | et,T!0>;
Pz:=T!0;
for j in [0..4] do Pz+:=coeffToT(B!Coefficient(P16,mm,2*j))*z^j; end for;
Cz:=ProjectiveClosure(Curve(AffineSpace(T),Pz));

print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_P16_QUOTIENT";
print "QUOTIENT_SHAPE",<TotalDegree(Pz),Degree(Pz,et),Degree(Pz,z),#Terms(Pz)>;
print "QUOTIENT_IRREDUCIBLE",IsIrreducible(Pz);
print "QUOTIENT_NORMALIZED_GENUS",Genus(Cz);
assert Genus(Cz) eq 1;

pts:=Points(Cz:Bound:=PointBound);
print "QUOTIENT_PLANE_POINTS_BOUND",PointBound,"COUNT",#pts,"POINTS",pts;
smooth:=[];
for q in pts do
    singular:=IsSingular(q);
    print " QUOTIENT_POINT",q,"SINGULAR",singular;
    if not singular then Append(~smooth,q); end if;
end for;
print "QUOTIENT_SMOOTH_RATIONAL_POINTS",smooth;

if #smooth gt 0 then
    q0:=smooth[1];
    print "ELLIPTIC_MODEL_START_USING",q0;
    Eraw,mp:=EllipticCurve(Cz,q0);
    E,mmin:=MinimalModel(Eraw);
    print "ELLIPTIC_RAW",Eraw;
    print "ELLIPTIC_MINIMAL",E;
    print "ELLIPTIC_AINVARIANTS",aInvariants(E);
    print "ELLIPTIC_CONDUCTOR",Conductor(E);
    print "ELLIPTIC_TORSION",Invariants(TorsionSubgroup(E));
    if DoRank then print "ELLIPTIC_RANK_BOUNDS",RankBounds(E); end if;
    print "ELLIPTIC_MAP",mp;
else
    print "ELLIPTIC_MODEL_SKIPPED_NO_SMOOTH_Q_POINT";
end if;
print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_P16_QUOTIENT_DONE";
quit;
