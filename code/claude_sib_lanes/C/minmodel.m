// minmodel.m — reduced minimal model + extra data for the NEW curve #2
// representative (u,r) = (-23/75, -9025/3519), (s,m,n) = (2208,-8303,-7200)
P<x> := PolynomialRing(Rationals());
A := [1,1,1,2,2];
si := 2208; mi := -8303; ni := -7200;
B := [2*si^2-si*ni, 2*si^2+si*mi-2*si*ni-mi*ni, 2*si^2+si*mi-si*ni-mi*ni,
      -mi*ni, 4*si^2-4*si*ni-mi*ni];
printf "B = %o\n", B;
f := &*[A[i] + B[i]*x : i in [1..5]];
C := HyperellipticCurve(f);
Cmin := ReducedMinimalWeierstrassModel(C);
printf "minimal model: %o\n", Cmin;
printf "G2Invariants: %o\n", G2Invariants(C);
printf "Igusa: %o\n", IgusaInvariants(Cmin);
D := Integers()!Discriminant(Cmin);
printf "disc(min model) = %o = %o\n", D, Factorization(D);
J := Jacobian(Cmin);
T := TorsionSubgroup(J);
printf "torsion of min model = %o\n", Invariants(T);
// conductor (may be slow; bounded attempt)
try
    N := Conductor(Cmin);
    printf "conductor = %o = %o\n", N, Factorization(N);
catch e
    printf "conductor computation failed/slow\n";
end try;
quit;
