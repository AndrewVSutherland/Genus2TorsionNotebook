
//////////////////////////////////////////////////////////////////////
//  Algebraic second-halving conditions in the M_1(8,4) family.
//
//  Base family:
//      C: y^2 = f(x) = x*A(x)*B(x)
//  with n=1, R=m/n, and
//      t = (2R^2 + (1-w^2)R - 2w^2)/(4(w^2-1)).
//
//  The [4,8] tangent cover is represented by variables U,V,M,N:
//      h=A*B,
//      h - x*(M*x+N)^2 = c4*(x^2+U*x+V)^2,
//  where c4=lc(h).  Its half of T_x=[x,0] is
//      H_x = [a0, v0],
//      a0=x^2+U*x+V,
//      v0=-(x*(M*x+N)) mod a0 = (M*U-N)*x + M*V.
//
//  This script prints compact coefficient systems for:
//    1. [4,16]: halve the explicit order-8 point P_R=(-R,Y_R), 4P_R=T_A.
//    2. [8,8]: halve H_x, giving an independent second order-8 chain.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Rng<R,w,U,V,M,N,a,b,c,d,e,rho,sigma> := PolynomialRing(Q, 13, "grevlex");
K := FieldOfFractions(Rng);
PX<x> := PolynomialRing(K);


do_print_full := assigned print_full;
do_print_factors := assigned print_factors;

function Clear(q)
    if q eq 0 then
        return Rng!0;
    end if;
    return Numerator(q);
end function;

procedure PrintEquationSummary(label, eqs)
    print label, "count", #eqs;
    for i in [1..#eqs] do
        e0 := Clear(eqs[i]);
        print label, "eq", i-1, "total_degree", TotalDegree(e0), "terms", #Terms(e0);
        if do_print_full then
            print e0;
        end if;
        if do_print_factors then
            print "factorization", Factorization(e0);
        end if;
    end for;
end procedure;

t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
A := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
B := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);
f := x*A*B;
h := A*B;
c4 := Coefficient(h, 4);

Qfac := R^2 - (K!1/2)*R*w^2 + (K!1/2)*R - w^2;
Y_R := -2*R*(R-1)^2*Qfac/(w^2-1);
print "Family data";
print "t =", t;
print "Y_R =", Y_R;
print "f(-R)-Y_R^2 is zero", Evaluate(f, -R) - Y_R^2 eq 0;
print "";

// First [4,8] tangent-cover equations.
a0 := x^2 + U*x + V;
ell0 := x*(M*x + N);
F0 := h - x*(M*x+N)^2 - c4*a0^2;
E0 := [Coefficient(F0, i) : i in [0..4]];
PrintEquationSummary("FIRST_COVER", E0);
print "";

// [4,16]: P_R is divisible by 2 iff f-l^2 = c4*(x+R)*q^2.
// Here q is monic quadratic and l is quadratic.  The sign of l(-R)=+-Y_R
// is irrelevant because P_R is divisible iff -P_R is divisible.
q416 := x^2 + a*x + b;
ell416 := c*x^2 + d*x + e;
F416 := f - ell416^2 - c4*(x+R)*q416^2;
E416 := [Coefficient(F416, i) : i in [0..4]];
PrintEquationSummary("TARGET_416", E416);
print "";

// [8,8]: H_x is divisible by 2.  If q is the Mumford u-polynomial of
// the new half, then the cubic interpolation polynomial ell must satisfy
// ell == -v0 mod a0.  Write ell = -v0 + a0*(rho*x+sigma).  Its leading
// coefficient is rho, so the degree-6 leading term forces the scalar
// in f-ell^2 = scalar*a0*q^2 to be -rho^2.
v0 := (M*U - N)*x + M*V;
q88 := x^2 + a*x + b;
ell88 := -v0 + a0*(rho*x + sigma);
F88 := f - ell88^2 + rho^2*a0*q88^2;
E88 := [Coefficient(F88, i) : i in [0..5]];
PrintEquationSummary("TARGET_88", E88);

quit;
