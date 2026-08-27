P<x> := PolynomialRing(Rationals());
A := [1,1,1,2,2];
s := 2208; m := -8303; n := -7200;
B := [2*s^2-s*n, 2*s^2+s*m-2*s*n-m*n, 2*s^2+s*m-s*n-m*n, -m*n, 4*s^2-4*s*n-m*n];
f := &*[A[i] + B[i]*x : i in [1..5]];
C := HyperellipticCurve(f);
printf "B = %o\n", B;
Cm := ReducedMinimalWeierstrassModel(C);
printf "reduced minimal model: %o\n", Cm;
printf "G2Invariants: %o\n", G2Invariants(Cm);
fm, hm := HyperellipticPolynomials(Cm);
fi := 4*fm + hm^2;   // y^2 = 4f+h^2 integral model
Ci := HyperellipticCurve(fi);
Ji := Jacobian(Ci);
printf "torsion of integral model 4f+h^2: %o\n", Invariants(TorsionSubgroup(Ji));
D := Integers()!Discriminant(Cm);
printf "disc = %o\nfactorization = %o\n", D, Factorization(D);
printf "conductor-ish bad primes: %o\n", PrimeDivisors(D);
quit;
