//////////////////////////////////////////////////////////////////////
// Exact recovery-boundary analysis for the compact T_B-halving chart.
//
// The generic chart recovers
//
//   a = anum/aden,  aden=s*(r^2-K)-r.
//
// This script analyzes aden=anum=0, factors all components, imposes the
// square lift K=m^2, and determines whether H2 has a rational a.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q:=Rationals(); Z:=Integers();
R<a,s,r,K>:=PolynomialRing(Q,4,"grevlex");

A3:=a-3+4*s^2*r-2*K;
H1:=s*((a-3)*r^2+4*r-K*(a+3))-r*A3;
H2:=8*s^2*(2*s^2*r^2+2*(a-3)*r+2-K*(2*s^2-6))
    -A3^2-32*s^3*r;
aden:=Coefficient(H1,a,1);
anum:=-Evaluate(H1,[R!0,s,r,K]);
assert H1 eq aden*a-anum;

// On s!=0, aden=0 gives K=r^2-r/s.  The numerator then splits into
// r=1/s and r=s(2s+3).
RS<S,Rr>:=PolynomialRing(Q,2,"grevlex");
KRS:=FieldOfFractions(RS);
toRS:=hom<R -> KRS |
          KRS!0,KRS!S,KRS!Rr,KRS!(Rr^2)-KRS!Rr/KRS!S>;
nsub:=Numerator(toRS(anum));
print "CONTACT6_M612_TB_RECOVERY_BOUNDARY";
print "ADEN",aden;
print "ANUM",anum;
print "ANUM_ON_ADEN",Factorization(RS!nsub);

// Component C0: r=1/s, K=0.
US<u,z>:=PolynomialRing(Q,2,"grevlex");
KUS:=FieldOfFractions(US);
mapC0:=hom<R -> KUS | KUS!z,KUS!u,KUS!1/KUS!u,KUS!0>;
h2C0:=US!Numerator(mapC0(H2));
print "C0","r=1/s","K=0","H2",Factorization(h2C0);

// Component C1: r=s(2s+3).
rr1:=S*(2*S+3);
kk1:=(S+1)^2*(2*S-1)*(2*S+3);
QA<Avar,Svar>:=PolynomialRing(Q,2,"grevlex");
mapC1:=hom<R -> QA |
          Avar,Svar,Svar*(2*Svar+3),
          (Svar+1)^2*(2*Svar-1)*(2*Svar+3)>;
h2C1:=QA!mapC1(H2);
print "C1","r=s(2s+3)",
      "K=(s+1)^2(2s-1)(2s+3)";
print "C1_H2",Factorization(h2C1);
discA:=Discriminant(h2C1,Avar);
print "C1_H2_DISC",Factorization(discA);

// K is a square iff w^2=(2s-1)(2s+3).  Parametrize by
//   s=(t^2-2t+4)/(4t),  w=(t^2-4)/(2t).
QT<Avar2,t>:=PolynomialRing(Q,2,"grevlex");
KQT:=FieldOfFractions(QT);
sparam:=KQT!(t^2-2*t+4)/KQT!(4*t);
mapT:=hom<QA -> KQT | Avar2,sparam>;
h2T:=QT!Numerator(mapT(h2C1));
discT:=Discriminant(h2T,Avar2);
print "SQUARE_K_PARAM","s=(t^2-2t+4)/(4t)",
      "w=(t^2-4)/(2t)";
print "C1_SQUARE_K_H2_SHAPE",
      <TotalDegree(h2T),Degree(h2T,Avar2),Degree(h2T,t),#Terms(h2T)>;
print "C1_SQUARE_K_H2_FACTOR",Factorization(h2T);
print "C1_SQUARE_K_DISC_FACTOR",Factorization(discT);

// The remaining discriminant is square precisely when
// z^2=t^2+2t+4.  Parametrize this conic by
//   t=(3-u^2-2u)/(2u).
QU<Avar3,u>:=PolynomialRing(Q,2,"grevlex");
KQU:=FieldOfFractions(QU);
tparam:=KQU!(3-u^2-2*u)/KQU!(2*u);
mapU:=hom<QT -> KQU | KQU!Avar3,tparam>;
h2U:=QU!Numerator(mapU(h2T));
facU:=Factorization(h2U);
print "RATIONAL_A_PARAM","t=(3-u^2-2u)/(2u)";
print "C1_RATIONAL_A_H2_FACTOR",facU;

// Check smoothness of each rational a-branch by factoring the contact
// quintic discriminant over Q(u).
KUx<xx>:=PolynomialRing(KQU);
mapTU:=hom<QT -> KQU | KQU!0,tparam>;
sU:=mapTU(Numerator(sparam))/mapTU(Denominator(sparam));
rU:=sU*(2*sU+3);
kU:=(sU+1)^2*(2*sU-1)*(2*sU+3);
for fe in facU do
    pol:=fe[1];
    if Degree(pol,Avar3) ne 1 then continue; end if;
    c1:=Coefficient(pol,Avar3,1);
    c0:=Coefficient(pol,Avar3,0);
    aU:=-KQU!c0/KQU!c1;
    bU:=2*sU^2-3;
    q1:=xx*((bU+3)*xx^2+(aU-3)*xx+2)
          *(2*xx^2+(bU-3)*xx+(aU+3));
    discU:=QU!Numerator(Discriminant(q1));
    print "A_BRANCH",aU;
    print " S",sU,"R",rU,"K",kU;
    print " DISC_SHAPE",<Degree(discU,u),#Terms(discU)>;
    print " DISC_FACTOR",Factorization(discU);
end for;

print "CONTACT6_M612_TB_RECOVERY_BOUNDARY_DONE";
quit;
