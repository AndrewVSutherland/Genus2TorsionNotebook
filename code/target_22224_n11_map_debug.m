SetColumns(0);Q:=Rationals();R<T>:=PolynomialRing(Q);fixed:=[Q!-960,Q!1800,Q!2535];d0:=Q!-1352;t0:=Q!1/13;
Rx:=(fixed[2]+d0*T^2)/(fixed[2]+d0*t0^2);Ry:=(fixed[3]+d0*T^2)/(fixed[3]+d0*t0^2);
C:=HyperellipticCurve(Rx*Ry);P0:=C![t0,Q!1,Q!1];Eraw,phi:=EllipticCurve(C,P0);Einv:=Inverse(phi);
E,minmap:=MinimalModel(Eraw);minmapinv:=Inverse(minmap);
print "C",C;print "Eraw",Eraw;print "E",E;
print "mininv",DefiningPolynomials(minmapinv);print "Einv",DefiningPolynomials(Einv);
print "gens",Generators(E);TG,tm:=TorsionSubgroup(E);print "tors",Invariants(TG),[tm(g):g in TG];quit;
