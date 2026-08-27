// probe: square-discriminant locus of the X1(9) Kubert family
// (square disc <=> cyclic 2-division field <=> self-2-gluable along C3
//  <=> [9,9] on the glued Jacobian)
SetColumns(0);
Qt<t> := FunctionField(Rationals());
P<x> := PolynomialRing(Rationals());

for N in [7, 9] do
    if N eq 7 then
        bN := t^3-t^2; cN := t^2-t;
    else
        bN := t^2*(t-1)*(t^2-t+1); cN := t^2*(t-1);
    end if;
    E := EllipticCurve([1-cN, -bN, -bN, 0, 0]);
    d := Discriminant(E);
    h := P!(Numerator(d)) * P!(Denominator(d));
    fac := Factorization(h);
    unit := h div &*[ f[1]^f[2] : f in fac ];
    assert unit * &*[ f[1]^f[2] : f in fac ] eq h;
    printf "N=%o disc factorization: unit=%o\n", N, unit;
    for f in fac do printf "   (%o)^%o\n", f[1], f[2]; end for;
    ker := &*[ P | f[1] : f in fac | IsOdd(f[2]) ];
    // squarefree rational part of unit
    uq := Rationals()!unit;
    un := Numerator(uq)*Denominator(uq);
    us := Sign(un) * &*[ Integers() | p[1] : p in Factorization(Integers()!AbsoluteValue(un)) | IsOdd(p[2]) ];
    D := us * ker;
    printf "N=%o square class D(t) = %o  (degree %o)\n", N, D, Degree(D);
    if Degree(D) ge 3 then
        C := HyperellipticCurve(D);
        printf "N=%o curve y^2 = D(t): genus %o\n", N, Genus(C);
        pts := Points(C : Bound := 10000);
        printf "N=%o points to height 1e4: %o\n", N, pts;
    elif Degree(D) ge 1 then
        printf "N=%o D is degree %o -- conic/rational, points dense\n", N, Degree(D);
    end if;
end for;
quit;
