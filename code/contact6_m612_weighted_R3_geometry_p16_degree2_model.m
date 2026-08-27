//////////////////////////////////////////////////////////////////////
// Construct an exact binary-quartic model of the P16 normalization from
// the degree-2 place above (-6/17:0:1).
//
// If D is that degree-2 place, choose x in L(D) and y in
// L(2D)\<1,x,x^2>.  The nine functions
//
//   1,x,x^2,x^3,x^4,y,xy,x^2y,y^2
//
// lie in the 8-dimensional L(4D), giving the generalized binary quartic
// relation.  Completing the square produces Y^2=f_4(x).
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned PointBound then PointBound:=1000; end if;
if Type(PointBound) eq MonStgElt then PointBound:=StringToInteger(PointBound); end if;

Q:=Rationals();
A<e,mu,nu>:=PolynomialRing(Q,3,"grevlex");
PX<X>:=PolynomialRing(A);
D0:=3+5*e;
A1:=(-2*e-3*e^2)*X^2+(6*e+10*e^2)*X+(1-3*e^2);
A2:=2*e*X^2-(1+3*e);
R3:=2-3*X^2;
S:=R3*(mu*X+nu)^2-D0*A1*A2;
s4:=Coefficient(S,4); s3:=Coefficient(S,3);
s2:=Coefficient(S,2); s1:=Coefficient(S,1); s0:=Coefficient(S,0);
H1:=8*s4^2*s1-s3*(4*s4*s2-s3^2);
H0:=64*s4^3*s0-(4*s4*s2-s3^2)^2;
B<ee,mm>:=PolynomialRing(Q,2,"grevlex");
toB:=hom<A -> B | ee,mm,B!0>;
fac:=Factorization(B!toB(Resultant(H1,H0,nu)));
P16:=[f[1]:f in fac|Degree(f[1],ee) eq 16 and Degree(f[1],mm) eq 8][1];
Cp:=ProjectiveClosure(Curve(AffineSpace(B),P16));

q2:=Cp![-6/17,0,1];
pl2:=Places(q2)[1];
assert Degree(pl2) eq 2;
D2:=Divisor([<pl2,1>]);

print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_P16_DEGREE2_MODEL";
print "D2_DEGREE",Degree(D2);
V2,rho2:=RiemannRochSpace(D2);
assert Dimension(V2) eq 2;
f2:=[rho2(V2.i):i in [1..2]];
// Exactly one basis element is constant up to scaling.  Choose a
// nonconstant ratio, which still has poles bounded by D2.
assert f2[2] eq Parent(f2[1])!1;
assert f2[1] ne Parent(f2[1])!1;
xfun:=f2[1];
print "X_FUNCTION_IS_NONCONSTANT",true;

V4,rho4:=RiemannRochSpace(2*D2);
assert Dimension(V4) eq 4;
vone:=(Parent(xfun)!1) @@ rho4;
vx:=xfun @@ rho4;
vx2:=(xfun^2) @@ rho4;
W:=sub<V4|vone,vx,vx2>;
assert Dimension(W) eq 3;
yfun:=Parent(xfun)!0;
for i in [1..4] do
    if V4.i notin W then yfun:=rho4(V4.i); break; end if;
end for;
assert yfun ne 0;
print "Y_FUNCTION_CHOSEN",true;

V8,rho8:=RiemannRochSpace(4*D2);
assert Dimension(V8) eq 8;
funcs:=[Parent(xfun)|1,xfun,xfun^2,xfun^3,xfun^4,
                         yfun,xfun*yfun,xfun^2*yfun,yfun^2];
vecs:=[f @@ rho8:f in funcs];
M:=Matrix(Q,9,8,&cat[Eltseq(v):v in vecs]);
K:=Nullspace(M);
assert Dimension(K) eq 1;
c:=Eltseq(Basis(K)[1]);
assert c[9] ne 0;
c:=[q/c[9]:q in c];
print "GENERALIZED_QUARTIC_RELATION_COEFFICIENTS",c;

Qt<t>:=PolynomialRing(Q);
P2:=c[6]+c[7]*t+c[8]*t^2;
Q4:=c[1]+c[2]*t+c[3]*t^2+c[4]*t^3+c[5]*t^4;
f4:=P2^2-4*Q4;
print "BINARY_QUARTIC",f4;
print "BINARY_QUARTIC_FACTORIZATION",Factorization(f4);
assert Degree(f4) in {3,4};
G1:=GenusOneModel(f4);
Gmin,trmin,failed:=Minimise(G1);
print "MINIMISATION_POSITIVE_LEVEL_PRIMES",failed;
// The minimized binary quartic is already tiny; Reduce() is optional and
// can trigger installation-cache writes in sandboxed Magma installations.
Gred:=Gmin;
gmin_eq:=Equations(Curve(Gmin))[1];
Rg:=Parent(gmin_eq);
assert gmin_eq eq 6*Rg.1^4-2*Rg.2^4+Rg.3^2;
E:=MinimalModel(Jacobian(Gred));
print "TORSOR_MODEL_EQUATION",Equations(Curve(G1));
print "MINIMIZED_TORSOR_MODEL_EQUATION",Equations(Curve(Gmin));
print "WORKING_TORSOR_MODEL_EQUATION",Equations(Curve(Gred));
print "TORSOR_DISCRIMINANT",Discriminant(G1);
print "JACOBIAN",E;
print "JACOBIAN_AINVARIANTS",aInvariants(E);
print "JACOBIAN_RANK_BOUNDS",RankBounds(E);
print "JACOBIAN_TORSION",Invariants(TorsionSubgroup(E));
local_results:=AssociativeArray();
for p in [2,3,5,7,11,13,17,19,23,29,31,37,41] do
    try
        ok,pt:=IsLocallySolvable(Gred,p);
        local_results[p]:=ok;
        print "LOCAL_SOLVABILITY",p,ok;
        if ok then print " LOCAL_POINT",pt; end if;
    catch err
        print "LOCAL_SOLVABILITY_FAILED",p,err`Object;
    end try;
end for;
assert IsDefined(local_results,2) and not local_results[2];
assert IsDefined(local_results,3) and not local_results[3];
assert IsDefined(local_results,5) and local_results[5];
try
    pts:=Points(Curve(Gred):Bound:=PointBound);
    print "RATIONAL_POINTS_BOUND",PointBound,"COUNT",#pts,"POINTS",pts;
catch err
    print "RATIONAL_POINT_SEARCH_FAILED",err`Object;
end try;
print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_P16_DEGREE2_MODEL_DONE";
quit;
