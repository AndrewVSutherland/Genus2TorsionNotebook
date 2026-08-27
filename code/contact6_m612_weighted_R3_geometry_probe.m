//////////////////////////////////////////////////////////////////////
// Bounded exact geometry probe for the endpoint R3 halving cover.
// Computes coefficient shapes and the two low-variable-degree resultants.
// Heavy saturation/normalization is deliberately left to follow-up code
// after the projection factors are known.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q:=Rationals();
Z:=Integers();
A<e,mu,nu>:=PolynomialRing(Q,3,"grevlex");
PX<X>:=PolynomialRing(A);

D:=3+5*e;
A1:=(-2*e-3*e^2)*X^2+(6*e+10*e^2)*X+(1-3*e^2);
A2:=2*e*X^2-(1+3*e);
R3:=2-3*X^2;
S:=R3*(mu*X+nu)^2-D*A1*A2;
s:=[Coefficient(S,i):i in [4,3,2,1,0]];
s4,s3,s2,s1,s0:=Explode(s);
H1:=8*s4^2*s1-s3*(4*s4*s2-s3^2);
H0:=64*s4^3*s0-(4*s4*s2-s3^2)^2;

function Primitive(g)
    if g eq 0 then return g; end if;
    den:=LCM([Denominator(c):c in Coefficients(g)]);
    vals:=[Z!(den*c):c in Coefficients(g)];
    cont:=GCD(vals);
    if cont eq 0 then cont:=1; end if;
    h:=A!((Q!den/Q!Abs(cont))*g);
    if LeadingCoefficient(h) lt 0 then h:=-h; end if;
    return h;
end function;

procedure PrintFac(label,g)
    fac:=Factorization(Primitive(g));
    print label,[<TotalDegree(f[1]),Degree(f[1],e),Degree(f[1],mu),
                   Degree(f[1],nu),#Terms(f[1]),f[2]>:f in fac];
    for f in fac do
        if #Terms(f[1]) le 80 then print label cat "_FACTOR",f; end if;
    end for;
end procedure;

print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_PROBE";
for i in [1..5] do
    print "S_COEFF",5-i,"SHAPE",
          <TotalDegree(s[i]),Degree(s[i],e),Degree(s[i],mu),
           Degree(s[i],nu),#Terms(s[i])>,"POLY",s[i];
end for;
print "H1_SHAPE",<TotalDegree(H1),Degree(H1,e),Degree(H1,mu),
                   Degree(H1,nu),#Terms(H1)>;
print "H0_SHAPE",<TotalDegree(H0),Degree(H0,e),Degree(H0,mu),
                   Degree(H0,nu),#Terms(H0)>;
print "GCD_H1_H0",GCD(Primitive(H1),Primitive(H0));
PrintFac("S4_FACTORIZATION",s4);

print "RESULTANT_NU_START";
Rmu:=Primitive(Resultant(H1,H0,nu));
print "RESULTANT_NU_SHAPE",<TotalDegree(Rmu),Degree(Rmu,e),
                             Degree(Rmu,mu),#Terms(Rmu)>;
PrintFac("RESULTANT_NU_FACTORIZATION",Rmu);

print "RESULTANT_MU_START";
Rnu:=Primitive(Resultant(H1,H0,mu));
print "RESULTANT_MU_SHAPE",<TotalDegree(Rnu),Degree(Rnu,e),
                             Degree(Rnu,nu),#Terms(Rnu)>;
PrintFac("RESULTANT_MU_FACTORIZATION",Rnu);
print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_PROBE_DONE";
quit;
