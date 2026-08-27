// Verify HLP's two isolated exact-torsion split curves ([7,7] and [63])
SetColumns(0);
SetMemoryLimit(8*10^9);
P<x> := PolynomialRing(Rationals());
C1 := HyperellipticCurve(x^6 + 3025*x^4 + 3232987*x^2 + 869675859);
T1 := TorsionSubgroup(Jacobian(C1));
printf "HLP_77 TORSION %o order %o\n", Invariants(T1), #T1;
C2 := HyperellipticCurve(897*x^6 - 197570*x^4 + 79136353*x^2 - 146398496);
T2 := TorsionSubgroup(Jacobian(C2));
printf "HLP_63 TORSION %o order %o\n", Invariants(T2), #T2;
printf "HLP_VERIFY_DONE\n";
quit;
