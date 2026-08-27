// Exact check of the only small nontrivial row in the V1 N=10000 file,
// using the same seed/map/basis convention as the periodic run.
SetColumns(0); SetSeed(8);
Q:=Rationals(); Z:=Integers(); R<T>:=PolynomialRing(Q);
fixed:=[Q!-18,Q!20,Q!75]; d0:=Q!-1470/121;
Rx:=(fixed[1]+d0*T^2)/(fixed[1]+d0);
Ry:=(fixed[2]+d0*T^2)/(fixed[2]+d0);
C:=HyperellipticCurve(Rx*Ry); P0:=C![Q!1,Q!1,Q!1];
Eraw,phi:=EllipticCurve(C,P0); Einv:=Inverse(phi);
E,minmap:=MinimalModel(Eraw); minmapinv:=Inverse(minmap);
G1:=E![Q!3603,Q!216600,Q!1]; G2:=E![Q!93,Q!-2100,Q!1];
P:=-2*G1-G2; cp:=Einv(minmapinv(P)); tt:=Q!(cp[1]/cp[3]);
dd:=d0*tt^2; vals:=[(z+dd)/(z+d0):z in fixed];
sq:=[IsSquare(q):q in vals];
print "V1_SMALL_ROW",<-2,-1,1>,"P",P,"T",tt,
      "cover_values",vals,"cover_squares",sq;
if &and sq then
    den:=LCM([Denominator(z):z in fixed cat [dd]]);
    roots:=[Z!(den*z):z in fixed cat [dd]];
    g:=GCD([Abs(z):z in roots]); roots:=[z div g:z in roots];
    Pq<x>:=PolynomialRing(Q); f:=x*&*[x+(Q!z)^2:z in roots];
    TG,jmap:=TorsionSubgroup(Jacobian(HyperellipticCurve(f)));
    print "V1_SMALL_TARGET",roots,Invariants(TG),f;
end if;
quit;
