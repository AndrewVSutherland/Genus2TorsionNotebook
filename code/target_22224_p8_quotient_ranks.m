// Fast rank/map audit for the three elliptic quotients of the P8 fiber.
SetColumns(0); SetSeed(8);
Q:=Rationals(); R<T>:=PolynomialRing(Q);
fixed:=[Q!528,-726,-891]; d0:=Q!50; t1:=Q!19/5;
print "P8_RANKS_START";
for ij in [<1,2>,<1,3>,<2,3>] do
    i,j:=Explode(ij); x:=fixed[i]; y:=fixed[j];
    Rx:=(x+d0*T^2)/(x+d0); Ry:=(y+d0*T^2)/(y+d0);
    q:=Rx*Ry; C:=HyperellipticCurve(q); P0:=C![Q!1,Q!1,Q!1];
    vx:=Evaluate(Rx,t1); vy:=Evaluate(Ry,t1);
    ax,sx:=IsSquare(vx); ay,sy:=IsSquare(vy); assert ax and ay;
    P1:=C![t1,sx*sy,Q!1]; E,phi:=EllipticCurve(C,P0); Q1:=phi(P1);
    Emin,minmap:=MinimalModel(E); lo,hi:=RankBounds(Emin);
    TG,tmap:=TorsionSubgroup(Emin);
    print "P8_RANK",ij,"quartic",q,"rank_bounds",lo,hi,
          "torsion",Invariants(TG),"minimal",Emin,"known_point_E",Q1,
          "known_point_min",minmap(Q1),"inverse_quartic_map",DefiningPolynomials(Inverse(phi));
end for;
print "P8_RANKS_DONE"; quit;
