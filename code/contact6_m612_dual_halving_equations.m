//////////////////////////////////////////////////////////////////////
// Universal exact open-chart halving equations for the three rational
// 2-classes on the distinguished contact-6 Richelot dual.
//
// For g = Delta*R1*R2*R3 and u=Ri, put ell=u*(M*x+N).  Then
//
//   (ell^2-g)/u = u*(M*x+N)^2-Delta*Rj*Rk.
//
// Requiring this quartic to be a scalar square gives the exact two
// square-quartic covariants E1=E0=0.  This is the direct algebraic cover
// missing from the earlier resultant-square prefilter for R1.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
load "../../halving_mumford_tools.m";

Q:=Rationals();
A<a,b,M,N>:=PolynomialRing(Q,4,"grevlex");
PX<X>:=PolynomialRing(A);

Delta:=a*b+3*a+3*b+5;
R1:=(b^2-2*a-3)*X^2+(2*a*b+6*a+6*b+10)*X+(a^2-2*b-3);
R2:=2*X^2-(a+3);
R3:=2-(b+3)*X^2;

function PrimitiveA(g)
    if g eq 0 then return g; end if;
    den:=LCM([Denominator(c):c in Coefficients(g)]);
    vals:=[Integers()!(den*c):c in Coefficients(g)];
    cont:=GCD(vals); if cont eq 0 then cont:=1; end if;
    return A!((Q!den/Q!Abs(cont))*g);
end function;

procedure Report(label,u,v,w)
    S:=u*(M*X+N)^2-Delta*v*w;
    E1,E0:=SquareQuarticEquations(S);
    E1:=PrimitiveA(A!E1); E0:=PrimitiveA(A!E0);
    gg:=PrimitiveA(GCD(E1,E0));
    Scoeffs:=[A!Coefficient(S,i):i in [0..Degree(S)]];
    fac1:=Factorization(E1); fac0:=Factorization(E0);
    print "CLASS",label;
    print "S_SHAPE",<Degree(S),Maximum([TotalDegree(c):c in Scoeffs]),
          Maximum([Degree(c,M):c in Scoeffs]),
          Maximum([Degree(c,N):c in Scoeffs]),
          &+[#Terms(c):c in Scoeffs]>;
    print "E1_SHAPE",<TotalDegree(E1),Degree(E1,M),Degree(E1,N),#Terms(E1)>;
    print "E1_FACTOR_SHAPES",
          [<TotalDegree(fe[1]),Degree(fe[1],M),Degree(fe[1],N),
            #Terms(fe[1]),fe[2]>:fe in fac1];
    print "E0_SHAPE",<TotalDegree(E0),Degree(E0,M),Degree(E0,N),#Terms(E0)>;
    print "E0_FACTOR_SHAPES",
          [<TotalDegree(fe[1]),Degree(fe[1],M),Degree(fe[1],N),
            #Terms(fe[1]),fe[2]>:fe in fac0];
    print "COMMON_GCD_DEGREE",TotalDegree(gg),
          "GCD_FACTOR_SHAPES",
          [<TotalDegree(fe[1]),#Terms(fe[1]),fe[2]>:fe in Factorization(gg)];
end procedure;

print "CONTACT6_M612_DUAL_HALVING_EQUATIONS";
Report("R1_mixed",R1,R2,R3);
Report("R2_even",R2,R1,R3);
Report("R3_even",R3,R1,R2);
print "CONTACT6_M612_DUAL_HALVING_EQUATIONS_DONE";
quit;
