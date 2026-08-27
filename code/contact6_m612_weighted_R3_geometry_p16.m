//////////////////////////////////////////////////////////////////////
// Genus-one analysis of the P16 projection component carrying the
// open p=5 endpoint branches of the exact R3 halving curve.
//
// This script produces a GenusOneModel and its Jacobian when Magma's
// plane-curve normalization routines support the singular degree-16
// model.  Expensive rank/descent work is opt-in.
//
// Usage:
//   magma -b code/contact6_m612_weighted_R3_geometry_p16.m
//   magma -b PointBound:=1000 DoRank:=true code/..._p16.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned PointBound then PointBound:=200; end if;
if Type(PointBound) eq MonStgElt then PointBound:=StringToInteger(PointBound); end if;
if not assigned DoRank then DoRank:=false; end if;
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
P16:=[f[1]:f in fac|Degree(f[1],ee) eq 16 and
                         Degree(f[1],mm) eq 8][1];
assert IsIrreducible(P16);
Caff:=Curve(AffineSpace(B),P16);
Cp:=ProjectiveClosure(Caff);

print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_P16";
print "P16_SHAPE",<TotalDegree(P16),Degree(P16,ee),Degree(P16,mm),#Terms(P16)>;
print "P16_IRREDUCIBLE",IsIrreducible(P16);
print "P16_PROJECTIVE_DEGREE",Degree(Cp);
print "P16_NORMALIZED_GENUS",Genus(Cp);
print "P16_NONSINGULAR_PLANE_MODEL",IsNonsingular(Cp);
assert Genus(Cp) eq 1;

try
    sing:=SingularPoints(Cp);
    print "P16_RATIONAL_SINGULAR_POINT_COUNT",#sing;
    if #sing le 30 then
        for q in sing do
            print " P16_SINGULAR",q,"MULTIPLICITY",Multiplicity(Cp,q);
        end for;
    end if;
catch err
    print "P16_SINGULAR_POINTS_FAILED",err`Object;
end try;

try
    pts:=Points(Cp:Bound:=PointBound);
    print "P16_PLANE_RATIONAL_POINTS_BOUND",PointBound,"COUNT",#pts;
    if #pts le 50 then print "P16_PLANE_RATIONAL_POINTS",pts; end if;
catch err
    print "P16_PLANE_POINT_SEARCH_FAILED",err`Object;
end try;

// At e=0, P16 has the open factor mu^2+9.  This has two Q_5 roots
// because -9=1 mod 5, but no rational root over Q.
U<u>:=PolynomialRing(Q);
P16e0:=U!Evaluate(P16,[Q!0,u]);
print "P16_E0_FACTORIZATION",Factorization(P16e0);
assert IsDivisibleBy(P16e0,u^2+9);

try
    G1:=GenusOneModel(Cp);
    print "GENUS_ONE_MODEL_DEGREE",Degree(G1);
    print "GENUS_ONE_MODEL_EQUATIONS",Equations(Curve(G1));
    E:=MinimalModel(Jacobian(G1));
    print "JACOBIAN_MINIMAL_MODEL",E;
    print "JACOBIAN_INVARIANTS",aInvariants(E);
    print "JACOBIAN_CONDUCTOR",Conductor(E);
    print "JACOBIAN_TORSION",Invariants(TorsionSubgroup(E));
    if DoRank then
        print "JACOBIAN_RANK_BOUNDS",RankBounds(E);
    end if;
    for pp in [2,3,5,7,11,13,17,19,23,29,31,37,41] do
        try
            ok,locpt:=IsLocallySolvable(G1,pp);
            print "GENUS_ONE_LOCAL",pp,ok;
            if ok then print " LOCAL_POINT",locpt; end if;
        catch err
            print "GENUS_ONE_LOCAL_FAILED",pp,err`Object;
        end try;
    end for;
    try
        okreal,realpt:=IsLocallySolvable(G1,0);
        print "GENUS_ONE_LOCAL_REAL",okreal;
        if okreal then print " REAL_POINT",realpt; end if;
    catch err
        print "GENUS_ONE_LOCAL_REAL_FAILED",err`Object;
    end try;
    try
        gpts:=Points(Curve(G1):Bound:=PointBound);
        print "GENUS_ONE_RATIONAL_POINTS_BOUND",PointBound,"COUNT",#gpts;
        if #gpts le 50 then print "GENUS_ONE_RATIONAL_POINTS",gpts; end if;
    catch err
        print "GENUS_ONE_POINT_SEARCH_FAILED",err`Object;
    end try;
catch err
    print "GENUS_ONE_MODEL_FAILED",err`Object;
end try;

print "P5_BRANCH_STATEMENT",
      "open e=0 roots mu^2=-9 give Q_5 points on P16 but no e=0 Q-point";
print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_P16_DONE";
quit;
