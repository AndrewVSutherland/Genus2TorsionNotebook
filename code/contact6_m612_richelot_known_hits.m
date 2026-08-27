//////////////////////////////////////////////////////////////////////
// Distinguished Richelot dual of the known contact-6 [6,6] hits.
//
// The contact sextic has the rational quadratic-pair factorization
//
//   f = x * ((b+3)x^2+(a-3)x+2)
//         * (2x^2+(b-3)x+(a+3)).
//
// Pairing {0,infinity} and the roots of the two quadratics gives the
// pointwise-rational (2,2) Richelot kernel.  A degree-4 isogeny preserves
// the full 3-primary [3,3].  We test whether the dual gains 2-exponent 4,
// which would give target torsion [6,12].
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
P<x> := PolynomialRing(Q);

function IntegralSquareScale(f)
    d := LCM([ Denominator(Coefficient(f,i)) : i in [0..Degree(f)] ]);
    return P!(d^2*f), d;
end function;

function TorsionData(C)
    f, h := HyperellipticPolynomials(C);
    if Degree(h) ge 0 then
        error "expected y^2=f";
    end if;
    fi, d := IntegralSquareScale(P!f);
    Ji := Jacobian(HyperellipticCurve(fi));
    G, mp := TorsionSubgroup(Ji);
    return Invariants(G), fi;
end function;

function InvOrder(inv)
    ans := 1;
    for n in inv do ans *:= n; end for;
    return ans;
end function;

hits := [
    <"core_simple", 133/39, -7/13>,
    <"core_nonsimple_1", -19/9, 3/2>,
    <"core_nonsimple_2", -43/25, 1/8>,
    <"core_nonsimple_3", -15/8, 5/9>
];

for hit in hits do
    label, a, b := Explode(hit);
    A := x;
    B := (b+3)*x^2 + (a-3)*x + 2;
    C := 2*x^2 + (b-3)*x + (a+3);
    f := A*B*C;
    J := Jacobian(HyperellipticCurve(f));
    source_inv, source_fi := TorsionData(Curve(J));
    print "SOURCE", label, "a", a, "b", b;
    print "SOURCE_FACTORS", A, B, C;
    print "SOURCE_TORSION", source_inv, "order", InvOrder(source_inv);
    print "SOURCE_CURVE", source_fi;

    Rs := RichelotIsogenousSurfaces(J);
    print "RICHELOT_COUNT", #Rs;
    for i in [1..#Rs] do
        print "RICHELOT_TYPE", i, Type(Rs[i]);
        if Type(Rs[i]) eq JacHyp then
            inv, fi := TorsionData(Curve(Rs[i]));
            print "RICHELOT", i, "torsion", inv, "order", InvOrder(inv);
            print "RICHELOT", i, "factorization", Factorization(fi);
            print "RICHELOT", i, "curve", fi;
        else
            print "RICHELOT", i, Rs[i];
        end if;
    end for;
    print "END_SOURCE", label;
end for;

quit;
