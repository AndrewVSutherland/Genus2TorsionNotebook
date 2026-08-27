//////////////////////////////////////////////////////////////////////
//  Test whether the square condition for halving W_0 - infinity in
//  the M_1(8,2^w) odd model is sufficient.
//
//  The model is
//      y^2 = x*A(x)*B(x)
//  from NotesAndTodo, with parameters m,n,t.  The constant/leading
//  condition for the tangent equation is
//      (2*m^2 + m*n + 4*t)/(m*n + 2*n^2 + 4*t) = w^2.
//////////////////////////////////////////////////////////////////////

Q := Rationals();
P<x> := PolynomialRing(Q);

function FamilyPolynomial(m, n, w)
    t := (2*m^2 + (1-w^2)*m*n - 2*w^2*n^2)/(4*(w^2-1));
    A := n^4*x^2
         + (m^3*n + 4*m^2*t + m*n^3 - 8*m*n*t + 4*n^2*t)*x
         + m^4;
    B := (m*n + 2*n^2 + 4*t)*x^2
         + (m^2 + 4*m*n + n^2 + 8*t)*x
         + (2*m^2 + m*n + 4*t);
    return x*A*B, t, A, B;
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

tested := 0;
failed := 0;

for mi in [1..5] do
    for ni in [1..5] do
        for wi in [-5..5] do
            if wi in {-1,0,1} then
                continue;
            end if;
            m := Q!mi;
            n := Q!ni;
            w := Q!wi;
            f, t, A, B := FamilyPolynomial(m, n, w);
            if Degree(f) ne 5 or Discriminant(f) eq 0 then
                continue;
            end if;
            fI, L := IntegralModelPolynomial(f);
            C := HyperellipticCurve(fI);
            Jac := Jacobian(C);
            D := Jac![x, Q!0];
            divisible := IsDivisibleBy(D, 2);
            tested +:= 1;
            print "test", [mi,ni,wi], "t", t, "div", divisible;
            if not divisible then
                failed +:= 1;
                print "FAILED", [mi,ni,wi], f;
                break mi;
            end if;
        end for;
    end for;
end for;

print "tested", tested, "failed", failed;
quit;
