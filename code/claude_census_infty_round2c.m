// [2,4] fixed: (gamma,delta)=(1,1) -> q2 = x^2+3x+1 (disc 5, irreducible);
// q1 = x^2+(2t+1)x+t^2 (disc 4t+1, keep nonsquare and != 5-class).
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
q2 := x^2 + 3*x + 1;
for t in [Q!3, Q!5, Q!7, Q!-2, Q!1/2, Q!-4, Q!11, Q!13] do
    q1 := x^2 + (2*t+1)*x + t^2;
    f := x*q1*q2;
    if not IsSquarefree(f) then printf "  t=%o degenerate\n", t; continue; end if;
    den := LCM([Denominator(c) : c in Coefficients(f)]);
    C := HyperellipticCurve(den^2*f);
    T := TorsionSubgroup(Jacobian(SimplifiedModel(C)));
    inv := Sprint(Invariants(T));
    printf "  t=%o torsion %o", t, inv;
    if inv eq "[ 2, 4 ]" then
        sp, oks := StrictPrime(C);
        printf "  STRICT p=%o(%o)", sp, oks;
    end if;
    printf "\n";
end for;
printf "DONE_2C\n";
quit;
