//////////////////////////////////////////////////////////////////////
//  Reduced [8,8] second-halving conditions for the M_1(8,4) family.
//
//  This uses the structural simplification from
//
//      f - ell1^2 = -rho^2*a0*q^2
//
//  on the first [4,8] cover.  On the good open chart the constant
//  coefficient forces V to be a square.  Since V = eps*R^2*w on the
//  two first-cover signs, this means w=s^2 for eps=+1 and w=-s^2
//  for eps=-1.
//
//  With z=1/rho and lambda=tau/z, the second-halving equation becomes
//
//      z^2*x*c4*a0 - 2*z*ell0*(x+tau)
//        - a0*(x+tau)^2 + q^2 = 0,
//
//  where a0=x^2+U*x+V, ell0=x*(M*x+N), q=x^2+a*x+b.
//  The x^3 coefficient solves for a; the constant coefficient is
//      b^2 = V*tau^2.
//
//  This script prints the first-cover equations plus the two remaining
//  reduced second-halving equations after setting b=eta*mu*tau, where
//  mu^2=V and eta=+-1.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned mode then
    mode := "summary";
end if;

Q := Rationals();
Rng<R,s,U,M,N,z,tau,eta> := PolynomialRing(Q, 8, "grevlex");
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

procedure BuildCase(delta, eps, eta_sign)
    // delta gives w = delta*s^2.  For the open square subcover,
    // eps=+1 pairs with delta=+1 and eps=-1 pairs with delta=-1.
    w := (K!delta)*s^2;
    tt := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
    A := x^2 + (R^3 + 4*R^2*tt + R - 8*R*tt + 4*tt)*x + R^4;
    B := (R + 2 + 4*tt)*x^2 + (R^2 + 4*R + 1 + 8*tt)*x
         + (2*R^2 + R + 4*tt);
    h := A*B;
    c4 := Coefficient(h, 4);

    V := (K!eps)*R^2*w;
    mu := R*s;
    a0 := x^2 + U*x + V;
    ell0 := x*(M*x+N);

    first := h - x*(M*x+N)^2 - c4*a0^2;
    first_eqs := [Coefficient(first, i) : i in [1..4]];

    a := U/2 + M*z + tau - (c4/2)*z^2;
    b := (K!eta_sign)*mu*tau;
    q := x^2 + a*x + b;
    second := z^2*x*c4*a0 - 2*z*ell0*(x+tau)
              - a0*(x+tau)^2 + q^2;
    second_eqs := [Coefficient(second, i) : i in [1..2]];

    print "";
    print "CASE", "delta", delta, "eps", eps, "eta", eta_sign,
          "meaning", "w = delta*s^2, V = R^2*s^2";
    PrintSystem("FIRST_REDUCED", first_eqs);
    PrintSystem("SECOND_REDUCED_88", second_eqs);
end procedure;

if mode eq "summary" then
    for eta_sign in [1,-1] do
        BuildCase(1, 1, eta_sign);
        BuildCase(-1, -1, eta_sign);
    end for;
else
    print "unknown mode", mode;
end if;

quit;
