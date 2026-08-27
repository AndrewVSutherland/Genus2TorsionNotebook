//////////////////////////////////////////////////////////////////////
//  Symbolic coefficient equations for halving Q in the odd M_1(8,2,2)
//  model.
//
//  A reduced degree-2 class H can double to the theta divisor point
//  Q=(-1,y_Q) only if there are
//
//      a(X)   = X^2 + A*X + B,
//      ell(X) = C*X^2 + D*X + E,
//
//  with
//
//      f_{u,v}(X) - ell(X)^2 = L*(X+1)*a(X)^2.
//
//  This script forms the coefficient equations and tries inexpensive
//  eliminations/linear solves that can be used for a rational search on
//  the halving surface.
//
//  Typical run from torsion_jac:
//      magma code/m3222_halving_surface_symbolic.m
//////////////////////////////////////////////////////////////////////

Q := Rationals();
R<u,v,A,B,C,D,E,L> := PolynomialRing(Q, 8, "grevlex");
PX<X> := PolynomialRing(R);

qtilde := -X^2
    + (u^2*v - u^2 + u*v^2 - u - v^2 - v - 2)*X
    - (u^2 + u*v + v^2 + u + v + 1);
f := ((1-u)*X + 1)*((1-v)*X + 1)*((u+v+2)*X + 1)*qtilde;
a := X^2 + A*X + B;
ell := C*X^2 + D*X + E;
expr := f - ell^2 - L*(X+1)*a^2;
eqs := [Coefficient(expr, i) : i in [0..5]];

print "Halving surface coefficient equations";
print "variables: u,v,A,B,C,D,E,L";
for i in [1..#eqs] do
    print "eq", i-1, eqs[i];
end for;

// The top coefficient is linear in L and gives L = leading coefficient.
lead := eqs[6];
print "top equation factors", Factorization(lead);

// Substitute L from the top equation and inspect the remaining system.
Lval := Coefficient(f, 5);
eqs5 := [Evaluate(eqs[i], [u,v,A,B,C,D,E,Lval]) : i in [1..5]];
print "After L = f_5, degrees of remaining equations",
      [TotalDegree(e) : e in eqs5];
for i in [1..#eqs5] do
    print "eqL", i-1, eqs5[i];
end for;

if assigned do_groebner and do_groebner then
    I := ideal<R | eqs cat []>;
    print "Computing Groebner basis for the raw coefficient ideal...";
    time G := GroebnerBasis(I);
    print "Groebner basis length", #G;
    for i in [1..Minimum(#G, 5)] do
        print "GB", i, G[i];
    end for;
else
    print "Skipping Groebner basis; pass -b do_groebner:=true to compute it.";
end if;

quit;
