//////////////////////////////////////////////////////////////////////
// Low-degree rational chart for halving the marked 2-class T_B=[Q1,0].
//
// The raw square quartic is
//
//   S = Q1*(M*x+N)^2 - x*Q2.
//
// On a smooth open fiber M,N != 0.  Comparing its leading and constant
// coefficients shows that 2/(b+3) must be a square.  Write
//
//   b = 2*s^2-3,  M=1/m,  N=r/m.
//
// After multiplying by the rational square m^2, the quartic
//
//   Sbar = Q1*(x+r)^2-m^2*x*Q2
//
// has leading coefficient 2*s^2 and constant 2*r^2.  The sign in the
// constant coefficient of its square root is absorbed by s -> -s, so it is
// enough to impose
//
//   Sbar = 2*s^2*(x^2+p*x+r/s)^2.
//
// Coefficients x^3 and x^0 force p and the remaining two coefficients give
// H1=H2=0 below.  H1 is linear in a.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
Z := Integers();
A<a,s,r,m> := PolynomialRing(Q,4,"grevlex");
PX<X> := PolynomialRing(A);

b := 2*s^2-3;
Q1 := (b+3)*X^2+(a-3)*X+2;
Q2 := 2*X^2+(b-3)*X+(a+3);
Sbar := Q1*(X+r)^2-m^2*X*Q2;

A3 := a-3+4*s^2*r-2*m^2;
K := FieldOfFractions(A);
p := K!A3/K!(4*s^2);

H1 := s*((a-3)*r^2+4*r-m^2*(a+3))-r*A3;
H2 := 8*s^2*(2*s^2*r^2+2*(a-3)*r+2-m^2*(2*s^2-6))
      -A3^2-32*s^3*r;

// Exact coefficient reconstruction after clearing powers of s.
assert Coefficient(Sbar,4) eq 2*s^2;
assert Coefficient(Sbar,3) eq A3;
assert Coefficient(Sbar,0) eq 2*r^2;
assert s*Coefficient(Sbar,1)-r*A3 eq H1;
assert 8*s^2*Coefficient(Sbar,2)-A3^2-32*s^3*r eq H2;

aden := Coefficient(H1,a,1);
aconst := Evaluate(H1,[A!0,s,r,m]);
anum := -aconst;
assert H1 eq aden*a-anum;

da := Degree(H2,a);
Hred := &+[Coefficient(H2,a,i)*anum^i*aden^(da-i)
           : i in [0..da]];

print "CONTACT6_M612_TB_LOWCHART";
print "PARAMETERIZATION","b=2*s^2-3","Mraw=1/m","Nraw=r/m";
print "Sbar",Sbar;
print "A3",A3;
print "p",p;
print "H1",H1;
print "H2",H2;
print "H1_SHAPE",<TotalDegree(H1),Degree(H1,a),#Terms(H1)>;
print "H2_SHAPE",<TotalDegree(H2),Degree(H2,a),#Terms(H2)>;
print "A_RECOVERY_NUM",anum;
print "A_RECOVERY_DEN",aden;
print "HRED_SHAPE",<TotalDegree(Hred),Degree(Hred,s),Degree(Hred,r),
                    Degree(Hred,m),#Terms(Hred)>;
fac := Factorization(Hred);
print "HRED_FACTORIZATION",
      [<TotalDegree(fe[1]),Degree(fe[1],s),Degree(fe[1],r),
        Degree(fe[1],m),#Terms(fe[1]),fe[2]> : fe in fac];
for fe in fac do
    if #Terms(fe[1]) le 40 then print "HRED_FACTOR",fe; end if;
end for;
print "OPEN","s*r*m*aden*Disc(f) != 0";
print "MONIC_TOOL_RECOVERY",
      "Mtool=2*s^2/m, Ntool=2*s^2*r/m, G=x^2+p*x+r/s";
print "CONTACT6_M612_TB_LOWCHART_DONE";
quit;
