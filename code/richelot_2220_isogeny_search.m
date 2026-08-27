//////////////////////////////////////////////////////////////////////
//  Richelot / 2-power isogeny-class search around the known
//  geometrically simple [2,2,20] example.
//
//  The seed is the contact-5 specialization
//
//      h = 1 + t*x + ((t^2 - 1)/2)*x^2,
//      f = h^2 - ((t+1)^4/4)*x^5,
//      t = -8233/7225.
//
//  Magma's Richelot routines return rational-coefficient y^2=f(x)
//  models.  TorsionSubgroup wants integral coefficients, so for exact
//  torsion we replace f by d^2*f, where d clears denominators.  This is
//  the Q-isomorphic change of variables Y=d*y.
//////////////////////////////////////////////////////////////////////

Q := Rationals();
P<x> := PolynomialRing(Q);

function IntegralSquareScale(f)
    denoms := [ Denominator(Coefficient(f, i)) : i in [0..Degree(f)] ];
    d := LCM(denoms);
    return P!(d^2*f), d;
end function;

function TorsionInvariantsForCurve(C)
    f, h := HyperellipticPolynomials(C);
    if Degree(h) ge 0 then
        error "expected a model y^2=f(x)";
    end if;
    f_int, d := IntegralSquareScale(P!f);
    C_int := HyperellipticCurve(f_int);
    J_int := Jacobian(C_int);
    A, phi := TorsionSubgroup(J_int);
    return Invariants(A), f_int;
end function;

function InvOrder(inv)
    n := 1;
    for m in inv do
        n *:= m;
    end for;
    return n;
end function;

procedure PrintJacobianReport(label, J)
    C := Curve(J);
    inv, f_int := TorsionInvariantsForCurve(C);
    print label, "torsion", inv, "order", InvOrder(inv), "degree", Degree(f_int);
    print label, "factor_degrees", [ Degree(ff[1]) : ff in Factorization(f_int) ];
    print label, "curve", f_int;
    print label, "G2", G2Invariants(HyperellipticCurve(f_int));
end procedure;

t := -8233/7225;
h := 1 + t*x + ((t^2 - 1)/2)*x^2;
f_seed := h^2 - ((t+1)^4/4)*x^5;
C_seed := HyperellipticCurve(f_seed);
J_seed := Jacobian(C_seed);

print "Richelot search around known [2,2,20] seed";
print "seed_t", t;
PrintJacobianReport("seed", J_seed);

print "";
print "Immediate rational Richelot neighbors";
richelots := RichelotIsogenousSurfaces(J_seed);
print "immediate_count", #richelots;
for i in [1..#richelots] do
    label := Sprintf("richelot_%o", i);
    print label, "type", Type(richelots[i]);
    if Type(richelots[i]) eq JacHyp then
        PrintJacobianReport(label, richelots[i]);
    else
        print label, richelots[i];
    end if;
end for;

print "";
print "TwoPowerIsogenies traversal";
Js, products, weil_restrictions := TwoPowerIsogenies(J_seed);
print "twopower_jacobians", #Js,
      "twopower_products", #products,
      "twopower_weil_restrictions", #weil_restrictions;

best_order := 0;
best_labels := [];
for i in [1..#Js] do
    label := Sprintf("twopower_%o", i);
    C := Curve(Js[i]);
    inv, f_int := TorsionInvariantsForCurve(C);
    ord := InvOrder(inv);
    if ord gt best_order then
        best_order := ord;
        best_labels := [ label ];
    elif ord eq best_order then
        Append(~best_labels, label);
    end if;
    print label, "torsion", inv, "order", ord, "degree", Degree(f_int);
    print label, "factor_degrees", [ Degree(ff[1]) : ff in Factorization(f_int) ];
    print label, "curve", f_int;
    print label, "G2", G2Invariants(HyperellipticCurve(f_int));
end for;

for i in [1..#products] do
    print Sprintf("product_%o", i), products[i];
end for;

for i in [1..#weil_restrictions] do
    print Sprintf("weil_restriction_%o", i), weil_restrictions[i];
end for;

print "best_twopower_order", best_order, "best_labels", best_labels;
print "done";

quit;
