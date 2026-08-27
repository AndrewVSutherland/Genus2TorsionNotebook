// Validator: fully independent Magma check that ord(D_inf)=14 over Q for the
// B3 e2e "unplanned" rational candidate f = x^6+3/2x^4+5x^3+9/16x^2+31/4x+33/4
// (and the planted lambda-scaled f14). Method: Jacobian arithmetic, no CF.
SetColumns(0);
SetMemoryLimit(3*10^9);
Q := Rationals();
R<x> := PolynomialRing(Q);
fs := [ x^6 + (3/2)*x^4 + 5*x^3 + (9/16)*x^2 + (31/4)*x + 33/4,
        x^6 + (3/2)*x^4 + (1/2)*x^3 + (9/16)*x^2 + (1/8)*x + 1/16 ];
for f in fs do
    // integral model: x -> x/2, y -> y/8 clears denominators: 64*f(x/2)
    fi := 64*Evaluate(f, x/2);
    C := HyperellipticCurve(fi);
    Pinf := PointsAtInfinity(C);
    printf "npinf=%o disc_zero=%o\n", #Pinf, Discriminant(fi) eq 0;
    J := Jacobian(C);
    D := Pinf[1] - Pinf[2];
    N := Order(D);
    printf "DINF_ORDER f=%o N=%o\n", f, N;
end for;
print "ALLDONE";
quit;
