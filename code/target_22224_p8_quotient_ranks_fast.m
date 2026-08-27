//////////////////////////////////////////////////////////////////////
// Fast rank audit for the three elliptic quotients of the signed P8
// shared-triple fiber.  Deliberately avoids MordellWeilGroup on raw models.
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(8);
Q:=Rationals(); R<T>:=PolynomialRing(Q);
fixed:=[Q!528,Q!-726,Q!-891]; d0:=Q!50; t1:=Q!19/5;
SetLogFile("results/target_22224_p8_quotient_ranks_fast.log":Overwrite:=true);
print "P8_RANKS_FAST_START",fixed,d0,t1;
for ij in [<1,2>,<1,3>,<2,3>] do
    i,j:=Explode(ij); x:=fixed[i]; y:=fixed[j];
    Rx:=(x+d0*T^2)/(x+d0); Ry:=(y+d0*T^2)/(y+d0);
    C:=HyperellipticCurve(Rx*Ry); P0:=C![Q!1,Q!1,Q!1];
    sxok,sx:=IsSquare(Evaluate(Rx,t1)); syok,sy:=IsSquare(Evaluate(Ry,t1));
    assert sxok and syok;
    P1:=C![t1,sx*sy,Q!1];
    Eraw,phi:=EllipticCurve(C,P0); Emin,minmap:=MinimalModel(Eraw);
    Qmin:=minmap(phi(P1)); lo,hi:=RankBounds(Emin);
    TG,tmap:=TorsionSubgroup(Emin);
    print "P8_QUOTIENT_FAST",ij,"quartic",Rx*Ry,
          "minimal",Emin,"rank_bounds",lo,hi,
          "torsion",Invariants(TG),"known_point",Qmin,
          "known_order",Order(Qmin),"raw_to_min",minmap,
          "inverse_quartic_map",Inverse(phi);
end for;
print "P8_RANKS_FAST_DONE";
UnsetLogFile(); quit;
