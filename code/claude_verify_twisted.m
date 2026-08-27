QQ := Rationals();
Px<x> := PolynomialRing(QQ);
for t in [[2,3,12,18],[1,2,4,50],[34,41,68,82],[22,28,55,70],[29,121,125,145],[43,57,129,171],[38,52,133,182],[118,132,177,198]] do
    a,b,c,d := Explode([QQ!v : v in t]);
    f := x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2);
    J := Jacobian(HyperellipticCurve(f));
    A := TorsionSubgroup(J);
    // simplicity certificate: irreducible 12th-power-transformed L-poly at a good prime
    simple := false;
    p := 13;
    dd := Integers()!(a*b*c*d*&*[Integers()!((x2-x1)) : x1,x2 in [0,a^2,b^2,c^2,d^2] | x2 gt x1]);
    for tries in [1..30] do
        p := NextPrime(p);
        while dd mod p eq 0 do p := NextPrime(p); end while;
        Lp := EulerFactor(BaseChange(J, GF(p)));
        LT := Polynomial(Reverse(Coefficients(Lp)));
        if IsIrreducible(LT) then
            K<pi> := NumberField(LT);
            g := MinimalPolynomial(pi^12);
            if Degree(g) eq 4 and IsIrreducible(g) then simple := true; break; end if;
        end if;
    end for;
    printf "tuple %o: torsion %o, simple-cert(12th power, up to 30 primes): %o\n", t, Invariants(A), simple;
end for;
quit;
