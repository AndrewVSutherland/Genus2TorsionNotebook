//////////////////////////////////////////////////////////////////////
// Published order-11 genus-2 torsion families.
// Sources:
//   Flynn, Large Rational Torsion on Abelian Varieties.
//   Daowsud-Schmidt, Continued fractions for rational torsion.
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

procedure Check(label, f)
    print label;
    print "f", f;
    print "Discriminant zero", Discriminant(f) eq 0;
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    T, mp := TorsionSubgroup(J);
    print "TorsionSubgroup", Invariants(T);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43] do
        fp := ChangeRing(f, GF(p));
        if Discriminant(fp) eq 0 then continue; end if;
        Lp := LPolynomial(ChangeRing(C, GF(p)));
        fac := Factorization(Lp);
        if #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4 then
            print "irreducible Lp", p, Lp;
            break;
        end if;
    end for;
end procedure;

Check("Flynn t=1", Flynn11(1));
Check("Daowsud-Schmidt u=1", DaowsudSchmidt11(1));
