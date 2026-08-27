//////////////////////////////////////////////////////////////////////
// Universal direct-halving equations for the three nonzero rational
// 2-classes on the contact-6 family
//
// f = x*Q1*Q2,
// Q1=(b+3)x^2+(a-3)x+2,
// Q2=2x^2+(b-3)x+(a+3).
//
// For u in {x,Q1,Q2}, put ell=u*(M*x+N) and
//
//     S=(ell^2-f)/u.
//
// A rational half exists on the degree-4 open chart precisely when S is a
// scalar square of a quadratic.  The two square-quartic covariants E1,E0
// are printed and factored below.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
load "../../halving_mumford_tools.m";

Q := Rationals();
A<a,b,M,N> := PolynomialRing(Q,4,"grevlex");
PX<X> := PolynomialRing(A);

Q1 := (b+3)*X^2+(a-3)*X+2;
Q2 := 2*X^2+(b-3)*X+(a+3);
f := X*Q1*Q2;

function PrimitiveA(g)
    if g eq 0 then return g; end if;
    den := LCM([Denominator(c) : c in Coefficients(g)]);
    vals := [Integers()!(den*c) : c in Coefficients(g)];
    cont := GCD(vals);
    if cont eq 0 then cont := 1; end if;
    return A!((Q!den/Q!Abs(cont))*g);
end function;

procedure Report(label,u,S)
    assert S eq ExactQuotient((u*(M*X+N))^2-f,u);
    E1,E0 := SquareQuarticEquations(S);
    E1 := PrimitiveA(A!E1);
    E0 := PrimitiveA(A!E0);
    lc := PrimitiveA(A!Coefficient(S,4));
    gg := PrimitiveA(GCD(E1,E0));
    print "CLASS",label;
    print "U",u;
    print "S",S;
    print "S_LC",lc,"factorization",Factorization(lc);
    print "E1","degree",TotalDegree(E1),
          "MN_degree",<Degree(E1,M),Degree(E1,N)>,
          "terms",#Terms(E1);
    print E1;
    print "E1_FACTORIZATION",Factorization(E1);
    print "E0","degree",TotalDegree(E0),
          "MN_degree",<Degree(E0,M),Degree(E0,N)>,
          "terms",#Terms(E0);
    print E0;
    print "E0_FACTORIZATION",Factorization(E0);
    print "COMMON_GCD",gg,"factorization",Factorization(gg);

    // Remove only the literal common factor.  Any factor of the leading
    // coefficient is a separate degree-drop chart, not an open component.
    if TotalDegree(gg) gt 0 then
        E1 := ExactQuotient(E1,gg);
        E0 := ExactQuotient(E0,gg);
    end if;
    print "OPEN_CORE_SHAPES",
          <TotalDegree(E1),#Terms(E1),Degree(E1,M),Degree(E1,N)>,
          <TotalDegree(E0),#Terms(E0),Degree(E0,M),Degree(E0,N)>;
    print "OPEN_CORE_GCD",GCD(E1,E0);
end procedure;

SA := X*(M*X+N)^2-Q1*Q2;
SB := Q1*(M*X+N)^2-X*Q2;
SC := Q2*(M*X+N)^2-X*Q1;

print "CONTACT6_M612_HALVING_EQUATIONS";
Report("T_A_x",X,SA);
Report("T_B_Q1",Q1,SB);
Report("T_C_Q2",Q2,SC);
print "CONTACT6_M612_HALVING_EQUATIONS_DONE";
quit;
