//////////////////////////////////////////////////////////////////////
// Verify open finite-field points on the iterated contact-7 norm system.
// Each row is <p,a,b,B0,B1,r>.  The coefficients of the monic quartic A
// are recovered triangularly from degrees 7,6,5,4.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(10^9);

rows := [
    <5, 2,1,1,3,0>,
    <5, 2,1,2,2,2>,
    <11,5,7,3,0,7>,
    <11,6,7,4,5,6>
];

print "Z49_STRUCTURAL_FINITE_VERIFY";
for row in rows do
    p,a0,b0,u0,v0,r0 := Explode(row);
    K := GF(p); P<x> := PolynomialRing(K);
    a:=K!a0; b:=K!b0; u:=K!u0; v:=K!v0; r:=K!r0;
    h := 1-(K!7/2)*x+a*x^2+b*x^3;
    numerator := h^2+(x-1)^7;
    assert numerator mod x^2 eq 0;
    f := ExactQuotient(numerator,x^2);
    B := u+v*x;
    rhs := (x-1)*(x-r)^7;
    A := x^4;
    for k in [3,2,1,0] do
        err := A^2-B^2*f-rhs;
        z := -Coefficient(err,4+k)/2;
        A +:= z*x^k;
    end for;
    identity := A^2-B^2*f eq rhs;
    open := r ne 1 and Evaluate(h,1) ne 0 and
            Evaluate(B,1) ne 0 and Evaluate(B,r) ne 0;
    smooth := Degree(f) eq 5 and Discriminant(f) ne 0;
    print "row",row,"identity",identity,"open",open,"smooth",smooth;
    if not (identity and open and smooth) then continue; end if;
    C := HyperellipticCurve(f); J := Jacobian(C);
    hp := Evaluate(h,1);
    sr := -Evaluate(A,r)/Evaluate(B,r);
    assert Evaluate(f,r) eq sr^2;
    D := J![x-1,hp];
    Q := J![x-r,sr];
    relplus := 7*Q eq D;
    relminus := 7*Q eq -D;
    print "  f",f;
    print "  A",A,"B",B,"R",<r,sr>;
    print "  order_D",Order(D),"order_Q",Order(Q),
          "7Q=D",relplus,"7Q=-D",relminus,"#J",#J;
end for;
print "Z49_STRUCTURAL_FINITE_VERIFY_DONE";
quit;
