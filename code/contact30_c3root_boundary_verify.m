//////////////////////////////////////////////////////////////////////
// Exact verification of the global survivors of the projective C3-root
// sieve.  The sieve intentionally retains every bad reduction disk.  Its
// only height-10000 survivors are R=1,2,3,7/3 on both parameter branches;
// this script proves that all eight are the global c=0 boundary.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
P<x> := PolynomialRing(Q);

function FamilyAtR(R, eps)
    denR := R^2-5;
    t := (5*R^2-20*R+19)/denR;
    Y := -2*(5*R^2-22*R+25)/denR;
    u := t^3;
    s := t^5+t^4+(Q!5/2)*t^3+(Q!1/2)*t
       + eps*t*(t-Q!1/2)*(t+1)*Y;
    C := (u^2+1)/(2*u);
    c := (u^2-1)/(2*u);
    denq := u^6+6*u^4*s-2*u^4+15*u^3*s-u*s^3+u^2;
    if denq eq 0 then
        return t,u,c,denq,P!0,false;
    end if;
    numq := 15*u^5+90*u^4+20*u^3*s-6*u^2*s^2+231*u^3
          + 2*u^2*s-15*u*s^2+90*u^2-20*u*s+15*u-2*s;
    q := numq/denq;
    A := (s+q)/2;
    B := (15-s*q)/2;
    C3 := 2*x^3+(A-3)*x^2+(B+3)*x+(C-1);
    return t,u,c,denq,C3,true;
end function;

boundaries := [Q!1,Q!2,Q!3,Q!7/3];
for R in boundaries do
    for eps in [-1,1] do
        t,u,c,denq,C3,defined := FamilyAtR(R,Q!eps);
        assert t in {Q!-1,Q!1};
        assert u in {Q!-1,Q!1};
        assert c eq 0;
        roots := defined select Roots(C3) else [];
        print "BOUNDARY", "R", R, "eps", eps, "t", t, "u", u,
              "c", c, "q_den", denq, "C3_defined", defined,
              "C3_factorization", defined select Factorization(C3) else [],
              "rational_roots", roots;
    end for;
end for;

// Exact global identity behind the list.
S<R> := PolynomialRing(Q);
t := (5*R^2-20*R+19)/(R^2-5);
print "NUMERATOR(t-1) =", Factorization(Numerator(t-1));
print "NUMERATOR(t+1) =", Factorization(Numerator(t+1));
assert {R-2,R-3} eq {fe[1] : fe in Factorization(Numerator(t-1))};
assert {R-1,R-Q!7/3} eq {fe[1] : fe in Factorization(Numerator(t+1))};

print "EXACT_BOUNDARY_VERIFICATION_OK";
quit;
