// Validator: find fresh random chart-U' tuples whose class order IS divisible by 31
// (positive-path test data for the B2 ladder). Seeded Magma RNG, independent of builders.
SetColumns(0);
SetMemoryLimit(3*10^9);
SetSeed(20260730);
primes := [101, 211, 499, 1009, 2003];
found := 0; tried := 0;
while found lt 4 and tried lt 20000 do
    tried +:= 1;
    P := primes[1 + (tried mod 5)];
    Fp := GF(P); R<x> := PolynomialRing(Fp);
    t := [ Random(Fp) : i in [1..9] ];
    u := x^2 + t[1]*x + t[2]; v := t[3]*x + t[4];
    w := t[5]*x^4 + t[6]*x^3 + t[7]*x^2 + t[8]*x + t[9];
    f := v^2 + u*w;
    if Degree(f) ne 6 or Discriminant(f) eq 0 then continue; end if;
    if t[1]^2 - 4*t[2] eq 0 then continue; end if;
    ok := false;
    DD := 0;
    try
        J := Jacobian(HyperellipticCurve(f));
        DD := J![u, v]; ok := true;
    catch e
        ok := false;
    end try;
    if not ok then continue; end if;
    N := Order(DD);
    if N mod 31 eq 0 then
        found +:= 1;
        printf "ROW %o %o %o %o %o %o %o %o %o %o %o 1\n",
            P, t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], N;
    end if;
end while;
printf "TRIED %o FOUND %o\n", tried, found;
print "ALLDONE";
quit;
