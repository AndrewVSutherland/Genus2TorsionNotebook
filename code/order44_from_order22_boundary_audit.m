//////////////////////////////////////////////////////////////////////
// Projective/boundary audit for the order-44 scout at p=5 and p=17.
//
// For source_sign=e, the displayed parameter has a removable point at
// s=e and a genuine pole at s=-e.  The four base boundary charts are:
//
//   zero:     s=q,       X unchanged;
//   cancel:   s=e+q,     X unchanged (genuine smooth member);
//   pole:     s=-e+q,    X=q*Z;
//   infinity: s=1/q,     X=q^2*Z.
//
// In each chart h(q,Z)=g(s,lambda*Z) is regular at q=0.  For the three
// degenerate charts, this script checks whether h is a square through
// first order in q.  If
//
//     h = H0^2 + 2*q*H0*L1 + O(q^2),   ord_q(h5)>=2,
//
// then for q=p the exact square-quartic halving identity has a formal
// solution modulo p^2: ell=+/-(H0+q*L1), with u arbitrary to this
// order.  Such a seed is deliberately labelled FORMAL_DEGENERATE; it
// is not a Jacobian point on the square special fiber.
//
// No Groebner basis, resultant, or unbounded search is performed.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

if assigned primes then
    if Type(primes) eq MonStgElt then
        prime_list := [StringToInteger(a) : a in Split(primes, ",") | #a gt 0];
    else
        prime_list := primes;
    end if;
else
    prime_list := [5,17];
end if;

function IsDivisibleBy2Finite(D,G,phi)
    a := D @@ phi;
    coords := Eltseq(a);
    invs := Invariants(G);
    for i in [1..#coords] do
        if GCD(2,invs[i]) ne 1 and coords[i] mod 2 ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function ZeroOrder(poly)
    if poly eq 0 then return 10^6; end if;
    for i in [0..Degree(poly)] do
        if Coefficient(poly,i) ne 0 then return i; end if;
    end for;
    return 10^6;
end function;

function AtZero(c)
    num := Numerator(c);
    den := Denominator(c);
    assert Evaluate(den,0) ne 0;
    return Q!(Evaluate(num,0)/Evaluate(den,0));
end function;

function DerivativeAtZero(c)
    num := Numerator(c);
    den := Denominator(c);
    n0 := Evaluate(num,0);
    d0 := Evaluate(den,0);
    assert d0 ne 0;
    return Q!((Evaluate(Derivative(num),0)*d0
              - n0*Evaluate(Derivative(den),0))/d0^2);
end function;

function CoefficientsIntegralAtPrime(polys,p)
    for f in polys do
        for c in Coefficients(f) do
            if Denominator(Q!c) mod p eq 0 then return false; end if;
        end for;
    end for;
    return true;
end function;

function ReducePolynomial(f,F,P)
    X := P.1;
    return &+[ (F!Numerator(Q!Coefficient(f,i))/F!Denominator(Q!Coefficient(f,i)))*X^i
               : i in [0..Degree(f)] ];
end function;

function FamilyPolynomial(family,K,PX,par)
    X := PX.1;
    if family eq "Flynn" then
        return X^6 + 2*X^5 + (2*par+3)*X^4 + 2*X^3
               + (par^2+1)*X^2 + 2*par*(1-par)*X + par^2;
    end if;
    return X^6 - 4*X^5 + 8*(1+par)*X^4 - (10+32*par)*X^3
           + 8*(1+6*par+2*par^2)*X^2
           - 4*(1+6*par+16*par^2)*X + 64*par^2+1;
end function;

function ChartModel(family,e,chart)
    Kq<q> := FunctionField(Q);
    PX<X> := PolynomialRing(Kq);

    if chart eq "zero" then
        s := q; lam := q^0;
    elif chart eq "cancel" then
        s := Kq!e+q; lam := q^0;
    elif chart eq "pole" then
        s := -Kq!e+q; lam := q;
    else
        assert chart eq "infinity";
        s := 1/q; lam := q^2;
    end if;

    num := -s^2*(s^2+1)*(s^4-s^2+1) + 2*e*s^5;
    den := (s^2-1)^2;
    if family eq "Flynn" then
        par := num/den;
        r := s^2;
    else
        par := num/(4*den);
        r := 1+s^2;
    end if;
    f := FamilyPolynomial(family,Kq,PX,par);
    g := PX!(&+[Coefficient(f,i)*(r*X+1)^i*X^(6-i) : i in [0..6]]);
    h := PX!(&+[Coefficient(g,i)*lam^i*X^i : i in [0..Degree(g)]]);

    // The selected lambda makes every coefficient regular at q=0.
    assert &and [ZeroOrder(Numerator(Coefficient(h,i)))
                 ge ZeroOrder(Denominator(Coefficient(h,i)))
                 : i in [0..Degree(h)]];

    h0 := PolynomialRing(Q)!
          (&+[AtZero(Coefficient(h,i))*(PolynomialRing(Q).1)^i
              : i in [0..Degree(h)]]);
    h1 := PolynomialRing(Q)!
          (&+[DerivativeAtZero(Coefficient(h,i))*(PolynomialRing(Q).1)^i
              : i in [0..Degree(h)]]);
    h5 := Coefficient(h,5);
    ordh5 := ZeroOrder(Numerator(h5))-ZeroOrder(Denominator(h5));
    return h0,h1,ordh5,par,r;
end function;

procedure AuditChart(family,e,chart)
    h0,h1,ordh5,par,r := ChartModel(family,e,chart);
    print "CHART",family,"source_sign",e,"chart",chart;
    print "  special_h",h0;

    if chart eq "cancel" then
        print "  classification GENUINE_FINITE_MEMBER";
        print "  source_parameter_at_boundary",AtZero(par),
              "marked_root_at_boundary",AtZero(r),
              "disc_zero_over_Q",Discriminant(h0) eq 0;
        for p in prime_list do
            F := GF(p);
            PF<Xf> := PolynomialRing(F);
            if not CoefficientsIntegralAtPrime([h0],p) then
                print "  p",p,"bad_denominator_kept";
                continue;
            end if;
            hp := ReducePolynomial(h0,F,PF);
            if Degree(hp) ne 5 or Discriminant(hp) eq 0 then
                print "  p",p,"bad_reduction_kept";
                continue;
            end if;
            J := Jacobian(HyperellipticCurve(hp));
            D := J![Xf,F!1];
            G,phi := AbelianGroup(J);
            div2 := IsDivisibleBy2Finite(D,G,phi);
            print "  p",p,"good",true,"D_order",Order(D),
                  "J_invariants",Invariants(G),"D_divisible_by_2",div2,
                  "p2_exact_seeds",(div2 select "requires_lift_check" else "0");
        end for;
        return;
    end if;

    ok,H0 := IsSquare(h0);
    print "  classification DEGENERATE_SPECIAL_FIBER";
    print "  square_special_fiber",ok,"H0",H0,"ord_q_h5",ordh5;
    if not ok then
        print "  FORMAL_P2 status DEAD_not_square_at_order_0";
        return;
    end if;
    if Coefficient(H0,0) eq -1 then H0 := -H0; end if;
    divisible := IsDivisibleBy(h1,2*H0);
    if divisible then
        L1 := ExactQuotient(h1,2*H0);
    else
        L1 := Parent(h1)!0;
    end if;
    print "  first_derivative_h",h1;
    print "  square_through_first_order",divisible,"L1",L1;
    for p in prime_list do
        goodden := CoefficientsIntegralAtPrime([H0,L1],p);
        live := divisible and ordh5 ge 2 and goodden;
        print "  p",p,"FORMAL_DEGENERATE_P2_LIVE",live,
              "seed q=p ell=half_sign*(H0+p*L1) u=arbitrary_mod_p";
    end for;
end procedure;

print "ORDER44_FROM_ORDER22_BOUNDARY_AUDIT";
print "primes",prime_list;
for family in ["Flynn","DaowsudSchmidt"] do
    for e in [-1,1] do
        for chart in ["zero","cancel","pole","infinity"] do
            AuditChart(family,e,chart);
        end for;
    end for;
end for;

quit;
