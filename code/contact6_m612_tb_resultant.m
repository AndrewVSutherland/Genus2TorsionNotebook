//////////////////////////////////////////////////////////////////////
// Fast elimination of N from the marked rational 2-class halving cover.
//
// On the contact-6 family
//
//   f=x*Q1*Q2,
//   Q1=(b+3)x^2+(a-3)x+2,
//   Q2=2x^2+(b-3)x+(a+3),
//
// the marked order-6 class D satisfies 3D=[Q1,0].  A rational half of
// [Q1,0] supplies the 2-primary part needed to halve D, hence is the direct
// 2-primary condition for [6,12] after adding the independent 3-direction.
//
// Both square-quartic equations are quadratic in N.  This script uses their
// exact linear-recovery resultant instead of a four-variable Groebner basis.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
load "../../halving_mumford_tools.m";

Q := Rationals();
Z := Integers();
A<a,b,M,N> := PolynomialRing(Q,4,"grevlex");
PX<X> := PolynomialRing(A);

function PrimitiveA(g)
    if g eq 0 then return g; end if;
    den := LCM([Denominator(c) : c in Coefficients(g)]);
    vals := [Z!(den*c) : c in Coefficients(g)];
    cont := GCD(vals);
    if cont eq 0 then cont := 1; end if;
    return A!((Q!den/Q!Abs(cont))*g);
end function;

Q1 := (b+3)*X^2+(a-3)*X+2;
Q2 := 2*X^2+(b-3)*X+(a+3);
f := X*Q1*Q2;
S := Q1*(M*X+N)^2-X*Q2;
E1,E0 := SquareQuarticEquations(S);
E1 := PrimitiveA(A!E1);
E0 := PrimitiveA(A!E0);

a2:=Coefficient(E1,N,2); a1:=Coefficient(E1,N,1);
a0:=Coefficient(E1,N,0);
b2:=Coefficient(E0,N,2); b1:=Coefficient(E0,N,1);
b0:=Coefficient(E0,N,0);
assert Degree(E1,N) eq 2 and Degree(E0,N) eq 2;

Du := b2*a1-a2*b1;
Nu := a2*b0-b2*a0;
Cu := a1*b0-a0*b1;
Ru := Nu^2+Du*Cu;
assert a2*Nu^2+a1*Nu*Du+a0*Du^2 eq a2*Ru;
assert b2*Nu^2+b1*Nu*Du+b0*Du^2 eq b2*Ru;

D := PrimitiveA(Du);
Num := PrimitiveA(Nu);
Cross := PrimitiveA(Cu);
Rcore := PrimitiveA(Ru);

Res := PrimitiveA(Resultant(E1,E0,N));
assert IsDivisibleBy(Res,Rcore) or IsDivisibleBy(Rcore,Res);

print "CONTACT6_M612_TB_RESULTANT";
print "E1_SHAPE",<TotalDegree(E1),Degree(E1,N),Degree(E1,M),#Terms(E1)>;
print "E0_SHAPE",<TotalDegree(E0),Degree(E0,N),Degree(E0,M),#Terms(E0)>;
print "A2_FACTORIZATION",Factorization(PrimitiveA(a2));
print "B2_FACTORIZATION",Factorization(PrimitiveA(b2));
print "D_SHAPE",<TotalDegree(D),Degree(D,M),#Terms(D)>;
print "NUM_SHAPE",<TotalDegree(Num),Degree(Num,M),#Terms(Num)>;
print "CROSS_SHAPE",<TotalDegree(Cross),Degree(Cross,M),#Terms(Cross)>;
print "D_FACTORIZATION",
      [<TotalDegree(fe[1]),Degree(fe[1],M),#Terms(fe[1]),fe[2]>
       : fe in Factorization(D)];
print "NUM_FACTORIZATION",
      [<TotalDegree(fe[1]),Degree(fe[1],M),#Terms(fe[1]),fe[2]>
       : fe in Factorization(Num)];
print "RCORE_SHAPE",<TotalDegree(Rcore),Degree(Rcore,a),Degree(Rcore,b),
                     Degree(Rcore,M),#Terms(Rcore)>;
print "RESULTANT_SHAPE",<TotalDegree(Res),Degree(Res,a),Degree(Res,b),
                         Degree(Res,M),#Terms(Res)>;
facR := Factorization(Rcore);
facRes := Factorization(Res);
print "RCORE_FACTORIZATION",
      [<TotalDegree(fe[1]),Degree(fe[1],a),Degree(fe[1],b),
        Degree(fe[1],M),#Terms(fe[1]),fe[2]> : fe in facR];
print "RESULTANT_FACTORIZATION",
      [<TotalDegree(fe[1]),Degree(fe[1],a),Degree(fe[1],b),
        Degree(fe[1],M),#Terms(fe[1]),fe[2]> : fe in facRes];

for fe in facR do
    if #Terms(fe[1]) le 30 then
        print "RCORE_FACTOR",fe;
    end if;
end for;

print "GENERIC_RECOVERY","N=Nu/Du (use unnormalized pair)";
print "EXCEPTIONAL_BRANCH","D=Num=0 must be solved in N separately";
print "OPEN_BOUNDARY","M*(b+3)*Disc(f) != 0";

// Plug in every distinct recorded exact [6,6] source.  Linear factors in
// M are the only generic rational points of the univariate specialized
// projection; D=Nu=0 is reported separately.
seeds := [
    <"simple_h14",Q!133/39,Q!-7/13>,
    <"core_h6",Q!-19/9,Q!3/2>,
    <"core_h10_2",Q!-43/25,Q!1/8>,
    <"core_h10_3",Q!-15/8,Q!5/9>
];
UM<t> := PolynomialRing(Q);
for seed in seeds do
    label,av,bv := Explode(seed);
    sp := hom<A -> UM | av,bv,t,UM!0>;
    rr := UM!sp(Rcore);
    dd := UM!sp(Du);
    nn := UM!sp(Nu);
    facr := Factorization(rr);
    linear := [fe : fe in facr | Degree(fe[1]) eq 1];
    exceptional := GCD(dd,nn);
    print "SEED_PROJECTION",label,"degree",Degree(rr),
          "factor_degrees",[<Degree(fe[1]),fe[2]> : fe in facr],
          "rational_linear_factors",#linear,
          "exceptional_gcd_degree",Degree(exceptional);
    for fe in linear do
        mval := -Coefficient(fe[1],0)/Coefficient(fe[1],1);
        dval := Evaluate(dd,mval);
        if dval ne 0 then
            nval := Evaluate(nn,mval)/dval;
            print " SEED_RECOVER",label,"M",mval,"N",nval;
        else
            print " SEED_LINEAR_EXCEPTIONAL",label,"M",mval;
        end if;
    end for;
end for;

print "CONTACT6_M612_TB_RESULTANT_DONE";
quit;
