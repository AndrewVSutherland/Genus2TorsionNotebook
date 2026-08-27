SetColumns(0);
Q:=Rationals();
E:=EllipticCurve([Q!0,0,0,1520244,202972176]);
print "P8_E",E;
print "P8_E_RANK_BOUNDS",RankBounds(E);
print "P8_E_GENERATORS",Generators(E);
print "P8_E_TORSION",Invariants(TorsionSubgroup(E));
quit;
