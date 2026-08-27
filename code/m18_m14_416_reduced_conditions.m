//////////////////////////////////////////////////////////////////////
//  Reduced [4,16] second-halving conditions for the M_1(8,4) family.
//
//  The [4,8] tangent cover is
//
//      h - x*(M*x+N)^2 = c4*(x^2 + U*x + V)^2,      h=A*B.
//
//  The order-8 point whose fourth multiple is T_A is
//
//      P_R = (-R, Y_R).
//
//  To get [4,16], halve P_R.  Write
//
//      q   = x^2 + aa*x + bb,
//      ell = cc*x^2 + dd*x + ee.
//
//  The raw halving identity is
//
//      f - ell^2 = c4*(x+R)*q^2.                    (416)
//
//  On the good open chart bb*c4*R != 0, the constant coefficient forces
//
//      ss^2 = -c4*R,      ee = eta*ss*bb, eta = +/-1.
//
//  The x^4 coefficient solves for aa, and the x coefficient solves for dd.
//  The remaining [4,16] conditions are the x^2 and x^3 coefficients of
//  (416), plus the square relation ss^2+c4*R=0.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned mode then
    mode := "summary";
end if;

Q := Rationals();
Rng<R,w,U,V,M,N,ss,cc,bb> := PolynomialRing(Q, 9, "grevlex");
K := FieldOfFractions(Rng);
PX<x> := PolynomialRing(K);

do_full := assigned print_full;
do_factors := assigned print_factors;

function Clear(q)
    if q eq 0 then
        return Rng!0;
    end if;
    return Numerator(q);
end function;

procedure PrintSystem(label, eqs)
    print label, "count", #eqs;
    for i in [1..#eqs] do
        e := Clear(eqs[i]);
        print label, "eq", i-1, "degree", TotalDegree(e), "terms", #Terms(e);
        if do_full then
            print e;
        end if;
        if do_factors and e ne 0 then
            print "factorization", Factorization(e);
        end if;
    end for;
end procedure;

tt := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
A := x^2 + (R^3 + 4*R^2*tt + R - 8*R*tt + 4*tt)*x + R^4;
B := (R + 2 + 4*tt)*x^2 + (R^2 + 4*R + 1 + 8*tt)*x
     + (2*R^2 + R + 4*tt);
h := A*B;
f := x*h;
c4 := Coefficient(h, 4);

print "Reduced [4,16] data";
print "t =", tt;
print "c4 =", c4;
print "square relation: ss^2 + c4*R = 0";
print "open exclusions include bb*ss*c4*R*(w^2-1) != 0";
print "";

// First [4,8] cover.
a0 := x^2 + U*x + V;
first := h - x*(M*x+N)^2 - c4*a0^2;
first_eqs := [Coefficient(first, i) : i in [1..4]];
PrintSystem("FIRST_COVER", first_eqs);
print "";

// Reduced [4,16] equations.
h0 := Coefficient(h, 0);
h1 := Coefficient(h, 1);
h2 := Coefficient(h, 2);
h3 := Coefficient(h, 3);
square_eq := ss^2 + c4*R;

procedure BuildCase(eta_sign)
    ee := (K!eta_sign)*ss*bb;
    aa := (h3 - cc^2 - c4*R)/(2*c4);
    dd := (h0 - c4*bb^2 - 2*c4*R*aa*bb)/(2*ee);

    q416 := x^2 + aa*x + bb;
    ell416 := cc*x^2 + dd*x + ee;
    raw416 := f - ell416^2 - c4*(x+R)*q416^2;

    print "CASE eta", eta_sign;
    reduced416 := [square_eq, Coefficient(raw416, 2), Coefficient(raw416, 3)];
    PrintSystem("SECOND_REDUCED_416", reduced416);
    print "";

    // Sanity: these coefficients are solved by construction, modulo the
    // square relation.
    sanity := [Coefficient(raw416, i) : i in [0,1,4]];
    PrintSystem("SOLVED_COEFFICIENTS_416", sanity);
    print "";
end procedure;

BuildCase(1);
BuildCase(-1);

quit;
