//////////////////////////////////////////////////////////////////////
// Birational projection components of the saturated endpoint R3 curve.
//
// Factor Res_nu(H1,H0), discard endpoint-only and s4 factors, and for
// each of the two remaining irreducible factors P(e,mu) compute the gcd
// of H1,H0 in the function field Q(P)[nu].  A linear gcd proves that the
// space-curve component is birational to P.  Also inspect the quotient
// z=mu^2 and attempt exact plane-curve genus computations.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q:=Rationals(); Z:=Integers();
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
Rmu:=B!toB(Resultant(H1,H0,nu));
fac:=Factorization(Rmu);
Ps:=[f[1]:f in fac |
    (Degree(f[1],mm) in {4,8}) and
    GCD(f[1],B!toB(s4)) eq 1];
assert #Ps eq 2;
Sort(~Ps,func<a,b|Degree(a,mm)-Degree(b,mm)>);

T<et,z>:=PolynomialRing(Q,2,"grevlex");
coeffToT:=hom<B -> T | et,T!0>;

print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_COMPONENTS";
print "PROJECTION_COMPONENT_COUNT",#Ps;
for index in [1..#Ps] do
    P:=Ps[index];
    assert IsIrreducible(P);
    IP:=ideal<B|P>;
    assert IsPrime(IP);
    C:=quo<B|IP>;
    K:=FieldOfFractions(C);
    KN<t>:=PolynomialRing(K);

    function ToKN(h)
        ans:=KN!0;
        for j in [0..Degree(h,nu)] do
            cb:=B!toB(Coefficient(h,nu,j));
            ans+:=K!(C!cb)*t^j;
        end for;
        return ans;
    end function;

    h1:=ToKN(H1); h0:=ToKN(H0);
    gh:=GCD(h1,h0);
    gh:=gh/LeadingCoefficient(gh);
    print "COMPONENT",index,"PROJECTION_SHAPE",
          <TotalDegree(P),Degree(P,ee),Degree(P,mm),#Terms(P)>;
    print " COMPONENT_POLYNOMIAL",P;
    print " FUNCTION_FIELD_GCD_DEGREE_NU",Degree(gh);
    print " FUNCTION_FIELD_GCD",gh;

    // P is even in mu.  Build its quotient equation Pz(e,z).
    assert &and[Coefficient(P,mm,j) eq 0:j in [1..Degree(P,mm)]|IsOdd(j)];
    Pz:=T!0;
    for j in [0..Degree(P,mm) div 2] do
        Pz+:=coeffToT(B!Coefficient(P,mm,2*j))*z^j;
    end for;
    print " MU_SQUARE_QUOTIENT_SHAPE",
          <TotalDegree(Pz),Degree(Pz,et),Degree(Pz,z),#Terms(Pz)>;
    print " MU_SQUARE_QUOTIENT",Pz;
    print " MU_SQUARE_QUOTIENT_FACTORIZATION",Factorization(Pz);
    if Degree(Pz,z) eq 2 then
        disc:=Discriminant(Pz,z);
        print " QUADRATIC_Z_DISCRIMINANT_FACTORIZATION",Factorization(disc);
    end if;

    try
        Caff:=Curve(AffineSpace(B),P);
        Cp:=ProjectiveClosure(Caff);
        print " PROJECTIVE_DEGREE",Degree(Cp),"GENUS",Genus(Cp),
              "NONSINGULAR",IsNonsingular(Cp);
    catch err
        print " PROJECTIVE_GEOMETRY_FAILED",err`Object;
    end try;

    try
        Cz:=Curve(AffineSpace(T),Pz);
        Czp:=ProjectiveClosure(Cz);
        print " QUOTIENT_PROJECTIVE_DEGREE",Degree(Czp),
              "QUOTIENT_GENUS",Genus(Czp),
              "QUOTIENT_NONSINGULAR",IsNonsingular(Czp);
    catch err
        print " QUOTIENT_GEOMETRY_FAILED",err`Object;
    end try;
end for;
print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_COMPONENTS_DONE";
quit;
