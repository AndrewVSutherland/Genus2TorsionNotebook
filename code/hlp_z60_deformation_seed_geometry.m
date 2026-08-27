//////////////////////////////////////////////////////////////////////
// Lightweight geometric diagnostics for the split cyclic-[60] seed.
// Used to identify the reduced involution branch whose tangent is tested
// in hlp_z60_deformation_tangent.py.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB:=2; end if;
if Type(MemGB) eq MonStgElt then MemGB:=StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q:=Rationals(); P<x>:=PolynomialRing(Q);
f:=-185*(125*x^2-1728)*(2000*x^4-48525*x^2+294912);
C:=HyperellipticCurve(f);
print "HLP_Z60_DEFORMATION_SEED_GEOMETRY";
print "f",f;
try
    A:=AutomorphismGroup(C);
    print "AUTOMORPHISM_GROUP",A,"order",#A;
catch e
    print "AUTOMORPHISM_GROUP_UNAVAILABLE",e`Object;
end try;
try
    S:=Degree2Subcovers(C);
    print "DEGREE2_SUBCOVERS",#S;
    for i in [1..#S] do print "SUBCOVER",i,S[i]; end for;
catch e
    print "DEGREE2_SUBCOVERS_UNAVAILABLE",e`Object;
end try;
print "VISIBLE_REDUCED_INVOLUTION x |-> -x";
print "HLP_Z60_DEFORMATION_SEED_GEOMETRY_DONE";
quit;
