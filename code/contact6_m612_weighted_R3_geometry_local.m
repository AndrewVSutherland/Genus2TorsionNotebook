//////////////////////////////////////////////////////////////////////
// Finite-field local reconnaissance for the two saturated endpoint R3
// components.  Counts open H1=H0 points labelled by the primitive P8
// and P16 projection factors, with Jacobian ranks.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned primes then primes:=[3,5,7,11,13,19]; end if;
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
fac:=Factorization(B!toB(Resultant(H1,H0,nu)));
P8:=[f[1]:f in fac|Degree(f[1],ee) eq 8 and Degree(f[1],mm) eq 4][1];
P16:=[f[1]:f in fac|Degree(f[1],ee) eq 16 and Degree(f[1],mm) eq 8][1];

function Primitive(g)
    den:=LCM([Denominator(c):c in Coefficients(g)]);
    vals:=[Z!(den*c):c in Coefficients(g)];
    cont:=GCD(vals); if cont eq 0 then cont:=1; end if;
    return Parent(g)!((Q!den/Q!Abs(cont))*g);
end function;
P8:=Primitive(P8); P16:=Primitive(P16);

print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_LOCAL";
for p in primes do
    k:=GF(p);
    RA<a,m,n>:=PolynomialRing(k,3,"grevlex");
    redA:=hom<A -> RA | a,m,n>;
    RB<aa,mm>:=PolynomialRing(k,2,"grevlex");
    redB:=hom<B -> RB | aa,mm>;
    hh1:=redA(H1); hh0:=redA(H0); ss4:=redA(s4);
    pp8:=redB(P8); pp16:=redB(P16);
    counts:=AssociativeArray(); samples:=AssociativeArray();
    for label in ["P8","P16","BOTH","OTHER"] do
        counts[label]:=[0,0,0,0]; samples[label]:=[];
    end for;
    for av in k do for mv in k do for nv in k do
        pt:=[av,mv,nv];
        if Evaluate(ss4,pt) eq 0 then continue; end if;
        if Evaluate(hh1,pt) ne 0 or Evaluate(hh0,pt) ne 0 then continue; end if;
        in8:=Evaluate(pp8,[av,mv]) eq 0;
        in16:=Evaluate(pp16,[av,mv]) eq 0;
        label:=in8 and in16 select "BOTH" else
               (in8 select "P8" else (in16 select "P16" else "OTHER"));
        J:=Matrix(k,2,3,[Evaluate(Derivative(hh1,j),pt):j in [1..3]] cat
                         [Evaluate(Derivative(hh0,j),pt):j in [1..3]]);
        rk:=Rank(J);
        counts[label][1]+:=1; counts[label][rk+2]+:=1;
        if #samples[label] lt 8 then Append(~samples[label],<Z!av,Z!mv,Z!nv,rk>); end if;
    end for; end for; end for;
    print "PRIME",p;
    for label in ["P8","P16","BOTH","OTHER"] do
        print " ",label,"TOTAL_RANK0_RANK1_RANK2",counts[label],"SAMPLES",samples[label];
    end for;
end for;
print "CONTACT6_M612_WEIGHTED_R3_GEOMETRY_LOCAL_DONE";
quit;
