//////////////////////////////////////////////////////////////////////
//  Search the double-linear residual-quartic locus in the contact-5
//  order-20 family.
//
//  In the scaled residual quartic
//
//      u*(y^4 + 4*y^3 + 8*y^2 + 8*y + 4) - 4*y*(y+1)^2,
//
//  two rational roots z,w lead to the symmetric condition
//
//      (4-p)*s^2 + (-2*p^2 + 4*p + 8)*s
//          - p^3 + p^2 + 4*p + 4 = 0,
//
//  where s=z+w and p=z*w.  The discriminant in s is
//  4*p*(p-2)^2, so put p=r^2.  The final condition that z,w are rational
//  is the genus-2 curve
//
//      Y^2 = (r+1)*(r^2+2*r+2)*(r^3-r^2-4*r+2).
//
//  The non-boundary rational points on this curve give residual factor type
//  1+1+2 and are candidates for J(Q)_tors = [2,2,20].
//
//  Typical run:
//      magma -b height:=10000 code/contact5_order80_double_linear_search.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 1000;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

Q := Rationals();
P<x> := PolynomialRing(Q);

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            r := Q!num/den;
            key := Sprint(r);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, r);
            end if;
        end for;
    end for;
    return vals;
end function;

function TLinear(z)
    den := z^4 + 4*z^3 + 8*z^2 + 8*z + 4;
    if den eq 0 then
        return false, Q!0;
    end if;
    return true, -(z^4 + 4*z + 4)/den;
end function;

function FamilyPolynomial(t)
    b := (t^2 - 1)/2;
    h := 1 + t*x + b*x^2;
    f := h^2 - ((t + 1)^4/4)*x^5;
    return f, b;
end function;

function IntegralModel(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function FactorTypeString(fac)
    degs := Sort([ Degree(ff[1]) : ff in fac ]);
    return Join([ IntegerToString(d) : d in degs ], "+");
end function;

params := RationalParametersOfHeight(height);
checked := 0;
curve_points := [];
hits := [];
seen_t := {};

print "Contact-5/order-80 double-linear locus search";
print "height", height, "parameters", #params;

for r in params do
    checked +:= 1;
    F := (r+1)*(r^2+2*r+2)*(r^3-r^2-4*r+2);
    ok, Y := IsSquare(F);
    if not ok then
        continue;
    end if;

    Append(~curve_points, <r,Y>);
    print "CURVE_POINT", "r", r, "Y", Y;

    if r eq -2 then
        continue;
    end if;
    s := -(r^3 + r^2 + 2)/(r + 2);
    delta_root := Y/(r + 2);
    z := (s + delta_root)/2;
    w := (s - delta_root)/2;
    if z eq w then
        continue;
    end if;

    okz, t := TLinear(z);
    okw, tw := TLinear(w);
    if not okz or not okw or t ne tw or t eq -1 then
        continue;
    end if;

    tkey := Sprint(t);
    if tkey in seen_t then
        continue;
    end if;
    Include(~seen_t, tkey);

    f, b := FamilyPolynomial(t);
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        continue;
    end if;
    q := ExactQuotient(f, x - 1);
    facq := Factorization(q);
    ftype := FactorTypeString(facq);

    fI, L := IntegralModel(f);
    C := HyperellipticCurve(fI);
    J := Jacobian(C);
    G, phi := TorsionSubgroup(J);
    invs := Invariants(G);
    Append(~hits, <r,Y,z,w,t,b,invs,ftype,fI,facq>);

    print "HIT", "r", r, "Y", Y, "z", z, "w", w,
          "t", t, "b", b, "torsion", invs, "factor_type", ftype;
    print "  f =", fI;
    print "  residual_factorization =", facq;

    if Y ne 0 then
        Append(~curve_points, <r,-Y>);
    end if;
end for;

print "DONE height", height;
print "checked", checked, "curve_points_with_positive_Y_only", #curve_points,
      "unique_hits", #hits;

quit;
