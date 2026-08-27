// [2,2,2,6] upgrade: the M(2,2,2,6) chart (notes/m2226_order6_doubling.md)
// y^2 = x(x+2s^2-sn)(x+2s^2+sm-2sn-mn)(x+2s^2+sm-sn-mn)(2x-mn)(2x+4s^2-4sn-mn)
// marked [2,2,2,2] (six rational branch points) + order-6 class g.
SetColumns(0);
SetSeed(1);
SetMemoryLimit(8*10^9);
Q := Rationals();
P<x> := PolynomialRing(Q);
function StrictPrime(C)
    fC, hC := HyperellipticPolynomials(C);
    g := 4*fC + hC^2;
    dsc := Integers()!Numerator(Discriminant(g)/16);
    for p in PrimesInInterval(11, 600) do
        if dsc mod p eq 0 then continue; end if;
        Pp := P!Reverse(Coefficients(LPolynomial(ChangeRing(SimplifiedModel(C), GF(p)))));
        if not IsIrreducible(Pp) then continue; end if;
        K<pi> := NumberField(Pp);
        ok := true;
        for n in [1..12] do
            if Degree(MinimalPolynomial(pi^n)) ne 4 then ok := false; break; end if;
        end for;
        if ok then return p, true; end if;
    end for;
    return 0, false;
end function;
function M2226(s, m, n)
    return x*(x + 2*s^2 - s*n)*(x + 2*s^2 + s*m - 2*s*n - m*n)
            *(x + 2*s^2 + s*m - s*n - m*n)*(2*x - m*n)*(2*x + 4*s^2 - 4*s*n - m*n);
end function;
fibers := [ <25,-26,-15>, <7,2,1>, <3,-2,5>, <11,4,-3>, <5,7,2> ];
invs := {};
for fb in fibers do
    f := M2226(Q!fb[1], Q!fb[2], Q!fb[3]);
    if Degree(f) ne 6 or not IsSquarefree(f) then printf "  fiber %o degenerate\n", fb; continue; end if;
    C := HyperellipticCurve(f);
    Include(~invs, G2Invariants(C));
    T := TorsionSubgroup(Jacobian(SimplifiedModel(C)));
    printf "  fiber %o torsion %o\n", fb, Invariants(T);
    if Sprint(Invariants(T)) eq "[ 2, 2, 2, 6 ]" then
        sp, oks := StrictPrime(C);
        printf "  fiber %o strict prime %o (%o)\n", fb, sp, oks;
    end if;
end for;
printf "DISTINCT_INVS %o\n", #invs;
printf "DONE_C\n";
quit;
