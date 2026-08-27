//////////////////////////////////////////////////////////////////////
// Exact rank-one certificate for pair-scaling fiber P1.
//////////////////////////////////////////////////////////////////////
SetColumns(0);
Q:=Rationals(); PR<p>:=PolynomialRing(Q);
a:=Q!1; b:=Q!55; c:=Q!99; d:=Q!125; lam1:=Q!1089/25;
Kp:=(d+b)/(a+c); Kq:=(a+b)/(d+c);
Ap:=a*(Kp*c-b*p^2)+b*(d*p^2-Kp*a);
Bp:=-Kq*d*(Kp*c-b*p^2)-Kq*c*(d*p^2-Kp*a);
quartic:=-Ap*Bp;
C:=HyperellipticCurve(quartic);
A0:=Evaluate(Ap,Q!1); P0:=C![Q!1,A0,Q!1];

function ExactSquare(x)
    Z:=Integers(); x:=Q!x;
    okN,n:=IsSquare(Z!Numerator(x)); okD,m:=IsSquare(Z!Denominator(x));
    assert okN and okD; return Q!n/Q!m;
end function;
Ra:=((a*lam1+b)*(a*lam1+c))/((a+b)*(a+c));
Rb:=((a*lam1+b)*(d*lam1+b))/((a+b)*(d+b));
Rc:=((a*lam1+c)*(d*lam1+c))/((a+c)*(d+c));
ya:=ExactSquare(Ra); yb:=ExactSquare(Rb); yc:=ExactSquare(Rc);
p1:=ya/yb; q1:=yc/ya; P1:=C![p1,q1*Evaluate(Ap,p1),Q!1];

Eraw,phi:=EllipticCurve(C,P0); Emin,minmap:=MinimalModel(Eraw);
MW,toE:=MordellWeilGroup(Emin); inv:=Invariants(MW);
Qmin:=minmap(phi(P1)); relation:="not found in range";
for i in [0..1] do for j in [0..1] do for n in [-50..50] do
    if toE(i*MW.1+j*MW.2+n*MW.3) eq Qmin then
        relation:=Sprintf("%o*MW.1 + %o*MW.2 + %o*MW.3",i,j,n);
    end if;
end for; end for; end for;

print "P1_BASE_TUPLE",[a,b,c,d],"SECOND_LAMBDA",lam1;
print "Kp",Kp,"Kq",Kq;
print "A(p)",Ap;
print "B(p)",Bp;
print "QUARTIC",quartic;
print "BASE_POINT",P0,"SECOND_POINT",P1;
print "ERAW",Eraw;
print "EMIN",Emin;
print "RANK_BOUNDS",RankBounds(Emin);
print "MW_INVARIANTS",inv;
print "MW_GENERATORS",[toE(MW.i):i in [1..#inv]];
print "KNOWN_SECOND_ON_EMIN",Qmin;
print "KNOWN_SECOND_IN_MW",relation;
print "ORDER_KNOWN_SECOND",Order(Qmin);
print "CURVE_TO_ERAW_MAP",phi;
print "ERAW_TO_CURVE_MAP",Inverse(phi);
print "ERAW_TO_EMIN_MAP",minmap;
print "EMIN_TO_ERAW_MAP",Inverse(minmap);
quit;
