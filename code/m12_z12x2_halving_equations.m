//////////////////////////////////////////////////////////////////////
//  Symbolic equations for halving the independent rational 2-torsion
//  class in the one-parameter M(12) family a=(1-r)/4.
//
//  In the odd quintic model obtained by moving x=2 to infinity, write
//      f(X) = (X-beta) g(X).
//
//  A half of [beta - infinity] has a tangent polynomial H of degree
//  at most 2 and a monic quadratic U such that
//      f(X) - H(X)^2 = c (X-beta) U(X)^2,
//  where c is the leading coefficient of g.  Since H(beta)=0, write
//      H(X) = (X-beta)(mX+n).
//  Dividing by X-beta gives the four coefficient equations
//      g(X) - (X-beta)(mX+n)^2 = c (X^2 + A X + B)^2.
//
//  This script prints these equations for the known divisible 2-torsion
//  root beta_div and the independent root beta_ind.
//////////////////////////////////////////////////////////////////////

Q := Rationals();
K<r> := RationalFunctionField(Q);
Px<x> := PolynomialRing(K);

a := (1-r)/4;
T := a*x^2 - x + r;
h := (x-r)*(T+1);
W := h^2 + 4*a*x^2*T*(T+1);

PX<X> := PolynomialRing(K);
f5 := PX!0;
for i in [0..Degree(W)] do
    for j in [0..i] do
        f5 +:= Coefficient(W, i)*Binomial(i,j)*2^(i-j)*X^(6-j);
    end for;
end for;

print "f5 factorization:";
print Factorization(f5);
print "beta_div =", (1-r)/(4*r);
print "beta_ind =", (2-r)/(4*(r-1));

procedure PrintHalvingEquations(beta, label)
    g := ExactQuotient(f5, X-beta);
    c := Coefficient(g, 4);

    R<A,B,m,n> := PolynomialRing(K, 4, "grevlex");
    PR<Y> := PolynomialRing(R);
    gR := PR!0;
    for i in [0..4] do
        gR +:= R!Coefficient(g, i)*Y^i;
    end for;

    U := Y^2 + A*Y + B;
    E := gR - (Y - R!beta)*(m*Y+n)^2 - R!c*U^2;

    print "";
    print label;
    print "degree(E) =", Degree(E);
    print "coefficient equations, degrees 0..3:";
    for i in [0..3] do
        print i, Coefficient(E, i);
    end for;

    I := ideal<R | [Coefficient(E, i) : i in [0..3]]>;
    G := GroebnerBasis(I);
    print "Groebner basis length", #G;
    print "last four basis elements:";
    for i in [Max(1, #G-3)..#G] do
        print G[i];
    end for;
end procedure;

PrintHalvingEquations((1-r)/(4*r), "known divisible root beta_div");
PrintHalvingEquations((2-r)/(4*(r-1)), "independent root beta_ind");
quit;
