SetColumns(0);Q:=Rationals();Z:=Integers();E:=EllipticCurve([Q!0,0,0,-98064300,373711102000]);
print E;print "disc",Discriminant(E);print "vdisc",Valuation(Discriminant(E),13);
try print "local",LocalInformation(E,13);catch e print e;end try;
try print "kodaira",KodairaSymbol(E,13);catch e print e;end try;
print "ainv",aInvariants(E);print "cinv",cInvariants(E);print "j",jInvariant(E);
quit;
