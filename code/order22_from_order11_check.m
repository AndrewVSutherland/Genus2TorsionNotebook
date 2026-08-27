//////////////////////////////////////////////////////////////////////
// Order-22 families by forcing a rational branch point in order-11
// infinity-torsion families.
//
// If D_inf = infinity_+ - infinity_- has order 11 on an even sextic and
// W=(r,0) is a rational branch point, then 2(W-infinity_+) = -D_inf,
// so W-infinity_+ has order 22.
//////////////////////////////////////////////////////////////////////

Q := Rationals();
P<x> := PolynomialRing(Q);

function Flynn11(t)
    t := Q!t;
    return x^6 + 2*x^5 + (2*t + 3)*x^4 + 2*x^3
           + (t^2 + 1)*x^2 + 2*t*(1 - t)*x + t^2;
end function;

function DaowsudSchmidt11(u)
    u := Q!u;
    return x^6 - 4*x^5 + 8*(1 + u)*x^4 - (10 + 32*u)*x^3
           + 8*(1 + 6*u + 2*u^2)*x^2
           - 4*(1 + 6*u + 16*u^2)*x + 64*u^2 + 1;
end function;

function Flynn22Parameter(s, eps)
    s := Q!s;
    if s eq 0 or s^2 eq 1 then error "s must avoid 0 and +/-1"; end if;
    return (-s^2*(s^2 + 1)*(s^4 - s^2 + 1) + 2*eps*s^5)/(s^2 - 1)^2;
end function;

function DaowsudSchmidt22Parameter(s, eps)
    s := Q!s;
    if s eq 0 or s^2 eq 1 then error "s must avoid 0 and +/-1"; end if;
    return (-s^2*(s^2 + 1)*(s^4 - s^2 + 1) + 2*eps*s^5)/(4*(s^2 - 1)^2);
end function;

function PrimitiveIntegralModel(f)
    den := LCM([ Denominator(Coefficient(f, i)) : i in [0..Degree(f)] ]);
    F := den*f;
    cont := GCD([ Integers()!Coefficient(F, i) : i in [0..Degree(F)] ]);
    return F/cont;
end function;

procedure IrreducibleFrobeniusCertificate(C, f)
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47] do
        fp := ChangeRing(f, GF(p));
        if Discriminant(fp) eq 0 then continue; end if;
        Lp := LPolynomial(ChangeRing(C, GF(p)));
        fac := Factorization(Lp);
        if #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4 then
            print "irreducible Lp", p, Lp;
            return;
        end if;
    end for;
    print "no irreducible Lp in tested primes";
end procedure;

procedure Check(label, f, root)
    F := PrimitiveIntegralModel(f);
    print label;
    print "root", root, "f(root)", Evaluate(f, root);
    print "primitive", F;
    print "discriminant_zero", Discriminant(F) eq 0;
    C := HyperellipticCurve(F);
    J := Jacobian(C);
    T, mp := TorsionSubgroup(J);
    print "TorsionSubgroup", Invariants(T);
    IrreducibleFrobeniusCertificate(C, F);
end procedure;

for eps in [-1,1] do
    s := Q!2;
    t := Flynn22Parameter(s, eps);
    Check(Sprintf("Flynn22 s=%o eps=%o t=%o", s, eps, t), Flynn11(t), s^2);
end for;

for eps in [-1,1] do
    s := Q!2;
    u := DaowsudSchmidt22Parameter(s, eps);
    Check(Sprintf("DaowsudSchmidt22 s=%o eps=%o u=%o", s, eps, u),
          DaowsudSchmidt11(u), 1 + s^2);
end for;
