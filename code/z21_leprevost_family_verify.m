//////////////////////////////////////////////////////////////////////
// Exact verification of Leprévost's 1991 one-parameter family with a
// rational divisor class of order 21, plus the small specialization t=1.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemMB then MemMB:=200; end if;
if Type(MemMB) eq MonStgElt then MemMB:=StringToInteger(MemMB); end if;
SetMemoryLimit(MemMB*10^6);

Q:=Rationals(); K<t>:=FunctionField(Q); R<x>:=PolynomialRing(K);

p2:=t^14+4*t^13+19*t^12+32*t^11+113*t^10+188*t^9
    +379*t^8+448*t^7+379*t^6+188*t^5+113*t^4+32*t^3
    +19*t^2+4*t+1;
p1:=t^10+4*t^9+17*t^8+24*t^7+54*t^6+56*t^5+54*t^4
    +24*t^3+17*t^2+4*t+1;
p0:=t^6+4*t^5+15*t^4+16*t^3+15*t^2+4*t+1;
A21:=p2*x^2-2*(t^2+1)^2*p1*x+(t^2+1)^4*p0;
k21:=64*t^4*(t+1)^2*(t^2+1)^3
     *(t^4+2*t^3+6*t^2+2*t+1)^3;
f21:=A21^2-k21*x^3*(x-1)^2;

R24:=t^24+5*t^23+37*t^22+81*t^21+104*t^20-1211*t^19
    -6479*t^18-23783*t^17-59305*t^16-117982*t^15
    -191446*t^14-249302*t^13-286592*t^12-249302*t^11
    -191446*t^10-117982*t^9-59305*t^8-23783*t^7-6479*t^6
    -1211*t^5+104*t^4+81*t^3+37*t^2+5*t+1;
disc21:=-2^44*t^33*(t+1)^30*(t^2+1)^26*(t^2+t+1)
       *(t^4+2*t^3+6*t^2+2*t+1)^8*p0^3*R24;
assert Degree(f21) eq 5 and Discriminant(f21) eq disc21;

// The unique point at infinity is rational.  The marked point above x=0
// on the y=A21 branch gives the class used by Leprévost.
C21:=HyperellipticCurve(f21); J21:=Jacobian(C21);
D21:=J21![x,Evaluate(A21,0)];
assert 21*D21 eq J21!0;
assert 3*D21 ne J21!0 and 7*D21 ne J21!0;

print "Z21_LEPREVOST_FAMILY";
print "GENERIC_DEGREE",Degree(f21);
print "GENERIC_MARKED_CLASS_ORDER",21;
print "GENERIC_MODEL y^2=A21(x)^2-k21*x^3*(x-1)^2";

// At t=1, divide y and A21 by 128.  This gives a small integral model
// and the marked point (0,7).
P<X>:=PolynomialRing(Q);
f1:=-216*X^5+657*X^4-696*X^3+466*X^2-224*X+49;
f21_at_1:=&+[Q!(Evaluate(Numerator(Coefficient(f21,i)),1)
                  /Evaluate(Denominator(Coefficient(f21,i)),1))*X^i
             : i in [0..Degree(f21)]];
assert f21_at_1 eq 128^2*f1;
C1:=HyperellipticCurve(f1); J1:=Jacobian(C1);
D1:=J1![X,7];
assert Order(D1) eq 21;
G,mp:=TorsionSubgroup(J1);
assert Invariants(G) eq [21];
print "SPECIAL_T",1;
print "SPECIAL_F",f1;
print "SPECIAL_D",D1,"order",Order(D1);
print "SPECIAL_TORSION",Invariants(G),"order",#G;
for p in [5,11,13,17,19] do
    Cp:=ChangeRing(C1,GF(p));
    print "GOOD_PRIME",p,"J_ORDER",#Jacobian(Cp),
          "L_FACTOR",Factorization(LPolynomial(Cp));
end for;
print "Z21_LEPREVOST_FAMILY_VERIFY_DONE";
quit;
