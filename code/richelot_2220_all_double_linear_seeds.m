//////////////////////////////////////////////////////////////////////
//  Find all [2,2,20] seeds on the contact-5 double-linear locus and
//  run the Richelot / 2-power isogeny-class test on each.
//
//  The double-linear locus is controlled by the auxiliary genus-2 curve
//
//      Y^2 = (r+1)(r^2+2r+2)(r^3-r^2-4r+2).
//
//  A previous height search found only r=-2,-1,0,1/3, with only r=1/3
//  nondegenerate.  This script uses Chabauty on the auxiliary curve:
//  RankBounds(J)=[1,1], so Chabauty from an infinite-order point gives
//  all rational points on the auxiliary curve.  It then maps every
//  nondegenerate point to the contact-5 family and runs Magma's
//  TwoPowerIsogenies traversal.
//////////////////////////////////////////////////////////////////////

Q := Rationals();
P<x> := PolynomialRing(Q);

function IntegralSquareScale(f)
    denoms := [ Denominator(Coefficient(f, i)) : i in [0..Degree(f)] ];
    d := LCM(denoms);
    return P!(d^2*f), d;
end function;

function InvOrder(inv)
    n := 1;
    for m in inv do
        n *:= m;
    end for;
    return n;
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

function TLinear(z)
    den := z^4 + 4*z^3 + 8*z^2 + 8*z + 4;
    if den eq 0 then
        return false, Q!0;
    end if;
    return true, -(z^4 + 4*z + 4)/den;
end function;

function FamilyCurve(t)
    h := 1 + t*x + ((t^2 - 1)/2)*x^2;
    f := h^2 - ((t+1)^4/4)*x^5;
    return HyperellipticCurve(f), f;
end function;

function AffineCoordinates(pt)
    coords := Eltseq(pt);
    X := Q!coords[1];
    Y := Q!coords[2];
    Z := Q!coords[3];
    if Z eq 0 then
        return false, Q!0, Q!0;
    end if;
    return true, X/Z, Y/Z^3;
end function;

procedure PrintTwoPowerReport(label, C)
    J := Jacobian(C);
    seed_inv, seed_f := TorsionInvariantsForCurve(C);
    print label, "seed_torsion", seed_inv, "seed_order", InvOrder(seed_inv);
    print label, "seed_curve", seed_f;

    richelots := RichelotIsogenousSurfaces(J);
    print label, "immediate_richelot_count", #richelots;
    for i in [1..#richelots] do
        rlabel := Sprintf("%o_richelot_%o", label, i);
        print rlabel, "type", Type(richelots[i]);
        if Type(richelots[i]) eq JacHyp then
            inv, f_int := TorsionInvariantsForCurve(Curve(richelots[i]));
            print rlabel, "torsion", inv, "order", InvOrder(inv);
            print rlabel, "factor_degrees",
                  [ Degree(ff[1]) : ff in Factorization(f_int) ];
        end if;
    end for;

    Js, products, weil_restrictions := TwoPowerIsogenies(J);
    print label, "twopower_jacobians", #Js,
          "twopower_products", #products,
          "twopower_weil_restrictions", #weil_restrictions;

    best_order := InvOrder(seed_inv);
    best_labels := [ label cat "_seed" ];
    for i in [1..#Js] do
        tlabel := Sprintf("%o_twopower_%o", label, i);
        inv, f_int := TorsionInvariantsForCurve(Curve(Js[i]));
        ord := InvOrder(inv);
        if ord gt best_order then
            best_order := ord;
            best_labels := [ tlabel ];
        elif ord eq best_order then
            Append(~best_labels, tlabel);
        end if;
        print tlabel, "torsion", inv, "order", ord;
        print tlabel, "factor_degrees",
              [ Degree(ff[1]) : ff in Factorization(f_int) ];
    end for;

    print label, "best_order_including_seed", best_order,
          "best_labels", best_labels;
end procedure;

f_aux := (x+1)*(x^2+2*x+2)*(x^3-x^2-4*x+2);
C_aux := HyperellipticCurve(f_aux);
J_aux := Jacobian(C_aux);

print "All double-linear [2,2,20] seeds plus Richelot search";
print "auxiliary_curve", f_aux;
print "rank_bounds", RankBounds(J_aux);
T_aux, phi_aux := TorsionSubgroup(J_aux);
print "auxiliary_torsion", Invariants(T_aux);

// Infinite-order divisor from (0,2) and (1/3,40/27).
u := x*(x - 1/3);
v := Interpolation([Q!0, Q!1/3], [Q!2, Q!40/27]);
D := J_aux![u, v];
print "chabauty_generator", D, "order", Order(D);

points := Chabauty(D : ptC := C_aux![0,2,1]);
print "chabauty_points_count", #points;
print "chabauty_points", points;

seed_count := 0;
seen_t := {};
for pt in points do
    affine, r, Y := AffineCoordinates(pt);
    if not affine then
        print "POINT_INFINITY", pt;
        continue;
    end if;

    print "POINT_AFFINE", "pt", pt, "r", r, "Y", Y;
    if r eq -2 then
        print "  skip boundary r=-2";
        continue;
    end if;

    s := -(r^3 + r^2 + 2)/(r + 2);
    delta_root := Y/(r + 2);
    z := (s + delta_root)/2;
    w := (s - delta_root)/2;
    if z eq w then
        print "  skip collision z=w", "z", z;
        continue;
    end if;

    okz, t := TLinear(z);
    okw, tw := TLinear(w);
    if not okz or not okw or t ne tw then
        print "  skip inconsistent t", "z", z, "w", w, "okz", okz, "okw", okw;
        continue;
    end if;
    if t eq -1 then
        print "  skip degenerate t=-1", "z", z, "w", w;
        continue;
    end if;

    tkey := Sprint(t);
    if tkey in seen_t then
        print "  skip duplicate t", t;
        continue;
    end if;
    Include(~seen_t, tkey);

    C, f := FamilyCurve(t);
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        print "  skip singular family point", "t", t;
        continue;
    end if;

    inv, f_int := TorsionInvariantsForCurve(C);
    seed_count +:= 1;
    label := Sprintf("seed_%o", seed_count);
    print "SEED", label, "r", r, "Y", Y, "z", z, "w", w,
          "t", t, "torsion", inv, "order", InvOrder(inv);
    print label, "curve", f_int;
    PrintTwoPowerReport(label, C);
end for;

print "double_linear_seed_count", seed_count;
print "unique_t_count", #seen_t;
print "done";

quit;
