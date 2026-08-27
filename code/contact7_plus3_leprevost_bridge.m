//////////////////////////////////////////////////////////////////////
// A one-parameter curve on the direct contact-7 plus cubic-contact
// 3-torsion cover, obtained by normalizing Leprévost's Z/21 family.
//
// The starting family has
//   f=A21^2-lambda*x^3*(x-1)^2
// and D=[(0,A21(0))-infinity] of order 21.  The function y-A21 gives
//   3D+2P7=0, P7=[(1,A21(1))-infinity],
// hence P7=9D has order 7.  We recover the order-7 principal function
// by a [3/1] Pade expansion at P7 and use its linear factor as the new
// coordinate z.  This puts the family in the direct contact-7 form
//   h7^2-z^2*g=-(z-1)^7.
//
// The class 7D has order 3.  Its transformed Mumford pair is converted
// to the normalized cubic-contact identity
//   H3^2-q3^3=L3^2*g.
//
// Everything generic is checked over Q(t).  The t=1 specialization is
// checked exactly and has full rational torsion [21].
//
// Hard memory cap: 280 MB.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(280*10^6);

Q:=Rationals();
K<t>:=FunctionField(Q);
R<x>:=PolynomialRing(K);

a2:=t^14+4*t^13+19*t^12+32*t^11+113*t^10+188*t^9
    +379*t^8+448*t^7+379*t^6+188*t^5+113*t^4+32*t^3
    +19*t^2+4*t+1;
a1:=-2*(t^2+1)^2*(t^10+4*t^9+17*t^8+24*t^7+54*t^6
    +56*t^5+54*t^4+24*t^3+17*t^2+4*t+1);
a0:=(t^2+1)^4*(t^6+4*t^5+15*t^4+16*t^3+15*t^2+4*t+1);
A21:=a2*x^2+a1*x+a0;
lambda:=64*t^4*(t+1)^2*(t^2+1)^3
       *(t^4+2*t^3+6*t^2+2*t+1)^3;
f:=A21^2-lambda*x^3*(x-1)^2;
assert Degree(f) eq 5 and Discriminant(f) ne 0;

C:=HyperellipticCurve(f);
J:=Jacobian(C);
D:=J![x,a0];
P7:=J![x-1,Evaluate(A21,K!1)];
assert 21*D eq J!0 and 3*D ne J!0 and 7*D ne J!0;
assert P7 eq 9*D;
assert 7*P7 eq J!0 and P7 ne J!0;

// Taylor square root y(u)^2=f(1+u), choosing y(0)=A21(1).
S<u>:=PolynomialRing(K);
fshift:=&+[Coefficient(f,i)*(1+u)^i : i in [0..Degree(f)]];
ycoef:=[K|Evaluate(A21,K!1)];
for n in [1..7] do
    cross:=K!0;
    for i in [1..n-1] do
        cross+:=ycoef[i+1]*ycoef[n-i+1];
    end for;
    Append(~ycoef,(Coefficient(fshift,n)-cross)/(2*ycoef[1]));
end for;
assert &and[
  Coefficient((&+[ycoef[i+1]*u^i:i in [0..7]])^2-fshift,n) eq 0
  : n in [0..7]];

// ell=1+rho*u makes ell*y a cubic through precision u^6.
rho:=-ycoef[5]/ycoef[4];
ell:=1+rho*u;
ytrunc:=&+[ycoef[i+1]*u^i : i in [0..7]];
prod:=ell*ytrunc;
assert &and[Coefficient(prod,i) eq 0 : i in [4..6]];
hraw:=&+[Coefficient(prod,i)*u^i : i in [0..3]];
resid:=hraw^2-ell^2*fshift;
c7:=Coefficient(resid,7);
assert resid eq c7*u^7 and c7 ne 0;

// z=ell sends P7 to z=1 and the zero of ell to z=0.  At that zero,
// hraw is a square root of -c7/rho^7, so no extra square cover appears.
u0:=-1/rho;
s0:=Evaluate(hraw,u0);
assert s0 ne 0 and s0^2 eq -c7/rho^7;

KZ<z>:=PolynomialRing(K);
usub:=(z-1)/rho;
h7:=(&+[Coefficient(hraw,i)*usub^i:i in [0..Degree(hraw)]])/s0;
g:=(&+[Coefficient(fshift,i)*usub^i:i in [0..Degree(fshift)]])/s0^2;
assert Degree(h7) eq 3 and Degree(g) eq 5;
assert h7^2-z^2*g eq -(z-1)^7;
assert Coefficient(h7,0) eq 1 and Coefficient(h7,1) eq -7/2;
aa:=Coefficient(h7,2);
bb:=Coefficient(h7,3);
assert g eq z^5+(bb^2-7)*z^4+(2*aa*bb+21)*z^3
    +(aa^2-7*bb-35)*z^2+(-7*aa+2*bb+35)*z+(2*aa-35/4);

// Transform the order-3 class 7D to the z-model.
D3old:=7*D;
assert 3*D3old eq J!0 and D3old ne J!0;
ux:=D3old[1];
vx:=D3old[2];
xsub:=1+usub;
qraw:=&+[Coefficient(ux,i)*xsub^i:i in [0..Degree(ux)]];
q3:=qraw/LeadingCoefficient(qraw);
v3:=(&+[Coefficient(vx,i)*xsub^i:i in [0..Degree(vx)]])/s0;
assert Degree(q3) eq 2 and IsMonic(q3) and Degree(v3) le 1;
assert (v3^2-g) mod q3 eq 0;

// Recover L3 and H3=q3*(z+w3)+L3*v3.  The z^5 coefficient forces
// w3=(L3^2+U)/2.  Taking the gcd of the remaining coefficient
// equations is a univariate calculation over Q(t).
PL<Lvar>:=PolynomialRing(K);
PLZ<Z>:=PolynomialRing(PL);
qL:=&+[PL!Coefficient(q3,i)*Z^i:i in [0..Degree(q3)]];
vL:=&+[PL!Coefficient(v3,i)*Z^i:i in [0..Degree(v3)]];
gL:=&+[PL!Coefficient(g,i)*Z^i:i in [0..Degree(g)]];
UU:=PL!Coefficient(q3,1);
wL:=(Lvar^2+UU)/2;
HL:=qL*(Z+wL)+Lvar*vL;
contact:=HL^2-qL^3-Lvar^2*gL;
eqL:=[PL!Coefficient(contact,i):i in [0..5] |
     Coefficient(contact,i) ne 0];
core:=eqL[1];
for i in [2..#eqL] do core:=GCD(core,eqL[i]); end for;
core:=core/LeadingCoefficient(core);
validL:=[];
for rr in Roots(core) do
    lv:=rr[1];
    if &and[Evaluate(ff,lv) eq 0:ff in eqL] then
        Append(~validL,lv);
    end if;
end for;
assert #validL ge 1;
L3:=validL[1];
w3:=(L3^2+Coefficient(q3,1))/2;
H3:=q3*(z+w3)+L3*v3;
assert L3 ne 0 and H3^2-q3^3 eq L3^2*g;

Cg:=HyperellipticCurve(g);
Jg:=Jacobian(Cg);
D7g:=Jg![z-1,Evaluate(h7,K!1)];
D3g:=Jg![q3,(H3/L3) mod q3];
assert 7*D7g eq Jg!0 and D7g ne Jg!0;
assert 3*D3g eq Jg!0 and D3g ne Jg!0;

function Spec(c)
    return Q!(Evaluate(Numerator(c),1)/Evaluate(Denominator(c),1));
end function;

PQ<zz>:=PolynomialRing(Q);
g1:=&+[Spec(Coefficient(g,i))*zz^i:i in [0..Degree(g)]];
h1:=&+[Spec(Coefficient(h7,i))*zz^i:i in [0..Degree(h7)]];
q1:=&+[Spec(Coefficient(q3,i))*zz^i:i in [0..Degree(q3)]];
H1:=&+[Spec(Coefficient(H3,i))*zz^i:i in [0..Degree(H3)]];
L1:=Spec(L3);
assert g1 eq zz^5-227/36*zz^4+1271/81*zz^3
             -27733/1458*zz^2+298/27*zz-257/108;
assert h1 eq 1-7/2*zz+86/27*zz^2-5/6*zz^3;
assert q1 eq zz^2-2*zz+73/81;
assert H1 eq zz^3-53/18*zz^2+8/3*zz-997/1458;
// The Pade normalization has s0=-81/2 at t=1, so v3 is the negative
// of the earlier y/(81/2) convention and L3 correspondingly is +1/3.
assert L1 eq 1/3;
assert h1^2-zz^2*g1 eq -(zz-1)^7;
assert H1^2-q1^3 eq L1^2*g1;
assert IsIrreducible(g1) and Discriminant(g1) ne 0;

C1:=HyperellipticCurve(g1);
J1:=Jacobian(C1);
D71:=J1![zz-1,Evaluate(h1,Q!1)];
D31:=J1![q1,(H1/L1) mod q1];
assert Order(D71) eq 7 and Order(D31) eq 3;
yscale:=LCM([Denominator(Coefficient(g1,i)):i in [0..Degree(g1)]]);
g1int:=yscale^2*g1;
C1int:=HyperellipticCurve(g1int);
J1int:=Jacobian(C1int);
D71int:=J1int![zz-1,yscale*Evaluate(h1,Q!1)];
D31int:=J1int![q1,yscale*((H1/L1) mod q1)];
assert Order(D71int) eq 7 and Order(D31int) eq 3;
T1,mp1:=TorsionSubgroup(J1int);
assert Invariants(T1) eq [21];

print "CONTACT7_PLUS3_LEPREVOST_BRIDGE_GENERIC";
print "base_dimension",1,"parameter","t";
print "contact7_a_degrees",Degree(Numerator(aa)),Degree(Denominator(aa));
print "contact7_b_degrees",Degree(Numerator(bb)),Degree(Denominator(bb));
print "q3_coefficient_degrees",
      [<Degree(Numerator(Coefficient(q3,i))),
        Degree(Denominator(Coefficient(q3,i)))>:i in [0..2]];
print "L3_degrees",Degree(Numerator(L3)),Degree(Denominator(L3));
print "generic_D7_order",7,"generic_D3_order",3;
print "SPECIAL_T",1,"a",Spec(aa),"b",Spec(bb),"L",L1;
print "SPECIAL_G",g1;
print "SPECIAL_Q3",q1;
print "SPECIAL_H3",H1;
print "SPECIAL_TORSION",Invariants(T1),"order",#T1,
      "D7_order",Order(D71),"D3_order",Order(D31);
print "CONTACT7_PLUS3_LEPREVOST_BRIDGE_DONE";
quit;
