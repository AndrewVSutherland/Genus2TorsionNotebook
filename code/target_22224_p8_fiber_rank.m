// Independent rank audit of the new signed shared-triple P8 fiber.
SetColumns(0);
SetLogFile("results/target_22224_p8_fiber_rank.log" : Overwrite := true);
Q:=Rationals(); R<T>:=PolynomialRing(Q);
fixed:=[Q!528,-726,-891]; d0:=Q!50; t1:=Q!19/5;
print "P8_FIBER_START",fixed,d0,t1;
for ij in [<1,2>,<1,3>,<2,3>] do
    i,j:=Explode(ij); a:=fixed[i]; b:=fixed[j];
    Ra:=(a+d0*T^2)/(a+d0); Rb:=(b+d0*T^2)/(b+d0);
    q:=Ra*Rb; C:=HyperellipticCurve(q);
    P0:=C![Q!1,Q!1,Q!1];
    va:=Evaluate(Ra,t1); vb:=Evaluate(Rb,t1);
    oka,sa:=IsSquare(va); okb,sb:=IsSquare(vb);
    assert oka and okb;
    P1:=C![t1,sa*sb,Q!1];
    E,phi:=EllipticCurve(C,P0); Q1:=phi(P1); Emin:=MinimalModel(E);
    lo,hi:=RankBounds(Emin); MW,mwmap:=MordellWeilGroup(E);
    print "P8_QUOTIENT",ij,"quartic",q,"known_order",Order(Q1),
          "rank_bounds",lo,hi,"minimal",Emin,"MW",Invariants(MW),
          "known_point",Q1;
end for;
print "P8_FIBER_DONE";
quit;
