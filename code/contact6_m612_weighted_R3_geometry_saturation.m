//////////////////////////////////////////////////////////////////////
// Saturation and (optionally) primary decomposition of the exact
// endpoint R3 halving curve H1=H0=0 on the open square-quartic chart.
//
// Usage:
//   magma -b code/contact6_m612_weighted_R3_geometry_saturation.m
//   magma -b mode:=decompose code/contact6_m612_weighted_R3_geometry_saturation.m
//
// The default mode only computes the saturated ideal and its Groebner
// basis.  The caller should bound the decompose mode externally.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned mode then mode:="summary"; end if;

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

Iraw:=ideal<A|H1,H0>;
print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_SATURATION";
print "MODE",mode;
print "RAW_DIMENSION",Dimension(Iraw);
print "SATURATION_START";
time Isat:=Saturation(Iraw,ideal<A|s4>);
bb:=Basis(Isat);
dim,degs:=Dimension(Isat);
print "SATURATED_DIMENSION",dim,"COMPONENT_DEGREES",degs;
try
    print "SATURATED_DEGREE",Degree(Isat);
catch err
    print "SATURATED_DEGREE","unavailable_for_affine_nonhomogeneous_ideal";
end try;
print "SATURATED_BASIS_SIZE",#bb;
print "SATURATED_BASIS_SHAPES",
      [<TotalDegree(g),Degree(g,e),Degree(g,mu),Degree(g,nu),#Terms(g)>:g in bb];
if #bb le 12 then
    for i in [1..#bb] do
        if #Terms(bb[i]) le 100 then print " BASIS",i,bb[i]; end if;
    end for;
end if;

if mode eq "decompose" then
    print "PRIMARY_DECOMPOSITION_START";
    time comps:=PrimaryDecomposition(Isat);
    print "PRIMARY_COMPONENT_COUNT",#comps;
    for i in [1..#comps] do
        C:=comps[i]; bc:=Basis(C);
        cdim,cdegs:=Dimension(C);
        try
            cdegree:=Sprint(Degree(C));
        catch err
            cdegree:="unavailable";
        end try;
        print "COMPONENT",i,"DIMENSION",cdim,"DEGREES",cdegs,
              "DEGREE",cdegree,"PRIME",IsPrime(C),
              "BASIS_SHAPES",
              [<TotalDegree(g),Degree(g,e),Degree(g,mu),Degree(g,nu),#Terms(g)>:g in bc];
        if #bc le 8 then
            for g in bc do if #Terms(g) le 120 then print " COMPONENT_BASIS",g; end if; end for;
        end if;
    end for;
end if;
print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_SATURATION_DONE";
quit;
