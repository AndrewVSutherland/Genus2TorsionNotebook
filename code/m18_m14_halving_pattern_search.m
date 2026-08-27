//////////////////////////////////////////////////////////////////////
//  Pattern search for halving W_0 - infinity in the M_1(8,2^w) odd
//  model
//
//      y^2 = x*A(x)*B(x).
//
//  This compares exact divisibility with simple squareclass invariants
//  of A and B.
//////////////////////////////////////////////////////////////////////

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

if not assigned Bnd then
    Bnd := 8;
elif Type(Bnd) eq MonStgElt then
    Bnd := StringToInteger(Bnd);
end if;

function IsSquareQ(q)
    q := Q!q;
    if q lt 0 then
        return false;
    end if;
    return IsSquare(Numerator(q)) and IsSquare(Denominator(q));
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function Family(m, n, t)
    A := n^4*x^2
         + (m^3*n + 4*m^2*t + m*n^3 - 8*m*n*t + 4*n^2*t)*x
         + m^4;
    C := (m*n + 2*n^2 + 4*t)*x^2
         + (m^2 + 4*m*n + n^2 + 8*t)*x
         + (2*m^2 + m*n + 4*t);
    return A, C, x*A*C;
end function;

tested := 0;
hits := 0;
criterion := 0;
criterion_miss := 0;

for mi in [-Bnd..Bnd] do
    for ni in [-Bnd..Bnd] do
        for ti in [-Bnd..Bnd] do
            if mi eq 0 or ni eq 0 then
                continue;
            end if;
            m := Q!mi;
            n := Q!ni;
            t := Q!ti;
            A, C, f := Family(m,n,t);
            if Degree(f) ne 5 or Discriminant(f) eq 0 then
                continue;
            end if;
            fI, L := IntegralModelPolynomial(f);
            Curve := HyperellipticCurve(fI);
            Jac := Jacobian(Curve);
            D := Jac![x, Q!0];
            divisible := IsDivisibleBy(D, 2);
            tested +:= 1;

            A0 := Coefficient(A,0);
            A2 := Coefficient(A,2);
            C0 := Coefficient(C,0);
            C2 := Coefficient(C,2);
            ResAC := Resultant(A,C);
            h0 := Coefficient(A*C,0);
            h4 := Coefficient(A*C,4);

            crit := IsSquareQ(C0) and IsSquareQ(C2) and IsSquareQ(ResAC);
            if crit then
                criterion +:= 1;
                if not divisible then
                    criterion_miss +:= 1;
                    print "CRIT_MISS", [mi,ni,ti], "C0", C0, "C2", C2, "Res", ResAC;
                end if;
            end if;

            if divisible then
                hits +:= 1;
                print "HIT", [mi,ni,ti],
                      "A0", A0, "A2", A2, "C0", C0, "C2", C2,
                      "Res", ResAC,
                      "sq C0", IsSquareQ(C0),
                      "sq C2", IsSquareQ(C2),
                      "sq h0", IsSquareQ(h0),
                      "sq h0/h4", IsSquareQ(h0/h4),
                      "sq Res", IsSquareQ(ResAC);
            end if;
        end for;
    end for;
end for;

print "tested", tested, "hits", hits, "criterion", criterion, "criterion_miss", criterion_miss;
quit;
