//////////////////////////////////////////////////////////////////////
//  Richelot / 2-power isogeny sweep over many [2,20] sources in the
//  contact-5/order-20 family with an extra rational 2-torsion point.
//
//  The two extra-2 parametrized loci are:
//
//      linear:
//          t = -(z^4 + 4*z + 4)/(z^4 + 4*z^3 + 8*z^2 + 8*z + 4)
//
//      quadratic-quadratic:
//          t = -(r^6 - 2*r^5 + 2*r^4 - 4*r^3 + 4*r^2 - 8*r + 8)
//               /((r^2 - 2)^2*(r^2 - 2*r + 2)).
//
//  For each smooth source that exact-checks to torsion order at least
//  min_source_order, run immediate rational Richelot neighbors and
//  Magma's TwoPowerIsogenies traversal.  This tests whether a [2,20]
//  source can land on a Richelot neighbor with [2,2,20], [2,4,20],
//  [2,40], or larger.
//
//  Typical runs:
//      magma -b height:=8 max_sources_per_label:=10 \
//          code/richelot_contact5_extra2_sweep.m
//
//      magma -b height:=30 max_sources_per_label:=50 \
//          code/richelot_contact5_extra2_sweep.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 20;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

if not assigned max_sources_per_label then
    max_sources_per_label := 25;
elif Type(max_sources_per_label) eq MonStgElt then
    max_sources_per_label := StringToInteger(max_sources_per_label);
end if;

if not assigned min_source_order then
    min_source_order := 40;
elif Type(min_source_order) eq MonStgElt then
    min_source_order := StringToInteger(min_source_order);
end if;

if not assigned target_order then
    target_order := 80;
elif Type(target_order) eq MonStgElt then
    target_order := StringToInteger(target_order);
end if;

Q := Rationals();
P<x> := PolynomialRing(Q);

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for H in [1..B] do
        for den in [1..H] do
            for num in [-H..H] do
                if Maximum(Abs(num), den) ne H then
                    continue;
                end if;
                if GCD(num, den) ne 1 then
                    continue;
                end if;
                z := Q!num/den;
                key := Sprint(z);
                if key notin seen then
                    Include(~seen, key);
                    Append(~vals, z);
                end if;
            end for;
        end for;
    end for;
    return vals;
end function;

function IntegralSquareScale(f)
    denoms := [ Denominator(Coefficient(f, i)) : i in [0..Degree(f)] ];
    d := LCM(denoms);
    return P!(d^2*f), d;
end function;

function InvOrder(inv)
    if #inv eq 0 then
        return 1;
    end if;
    n := 1;
    for m in inv do
        n *:= m;
    end for;
    return n;
end function;

function InvExponent(inv)
    if #inv eq 0 then
        return 1;
    end if;
    e := 1;
    for m in inv do
        e := LCM(e, m);
    end for;
    return e;
end function;

function FactorTypeString(fac)
    degs := Sort([ Degree(ff[1]) : ff in fac ]);
    return Join([ IntegerToString(d) : d in degs ], "+");
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

function FamilyPolynomial(t)
    b := (t^2 - 1)/2;
    h := 1 + t*x + b*x^2;
    return h^2 - ((t + 1)^4/4)*x^5;
end function;

function TLinear(z)
    den := z^4 + 4*z^3 + 8*z^2 + 8*z + 4;
    if den eq 0 then
        return false, Q!0;
    end if;
    return true, -(z^4 + 4*z + 4)/den;
end function;

function TQuadraticQuadratic(r)
    den := (r^2 - 2)^2*(r^2 - 2*r + 2);
    if den eq 0 then
        return false, Q!0;
    end if;
    num := r^6 - 2*r^5 + 2*r^4 - 4*r^3 + 4*r^2 - 8*r + 8;
    return true, -num/den;
end function;

function TFromParameter(label, z)
    if label eq "linear" then
        return TLinear(z);
    end if;
    return TQuadraticQuadratic(z);
end function;

procedure ConsiderObject(~global_best_order, ~global_best_exponent,
                         ~global_best_records, ~interesting_records,
                         source_label, obj_label, inv, f_int)
    ord := InvOrder(inv);
    exp := InvExponent(inv);
    if ord gt global_best_order then
        global_best_order := ord;
        global_best_exponent := exp;
        global_best_records := [ <source_label, obj_label, inv, ord, exp> ];
    elif ord eq global_best_order then
        if exp gt global_best_exponent then
            global_best_exponent := exp;
        end if;
        Append(~global_best_records, <source_label, obj_label, inv, ord, exp>);
    end if;

    if ord ge target_order or exp gt 20 then
        Append(~interesting_records, <source_label, obj_label, inv, ord, exp, f_int>);
        print "INTERESTING", "source", source_label, "object", obj_label,
              "torsion", inv, "order", ord, "exponent", exp;
        print "INTERESTING_CURVE", f_int;
    end if;
end procedure;

procedure SweepSource(label, param, t, source_index,
                      ~global_best_order, ~global_best_exponent,
                      ~global_best_records, ~interesting_records)
    f := FamilyPolynomial(t);
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        print "SOURCE_SKIP_SINGULAR", "label", label, "param", param, "t", t;
        return;
    end if;

    q := ExactQuotient(f, x - 1);
    facq := Factorization(q);
    ftype := FactorTypeString(facq);

    C := HyperellipticCurve(f);
    source_inv, source_f_int := TorsionInvariantsForCurve(C);
    source_order := InvOrder(source_inv);
    source_exp := InvExponent(source_inv);
    source_label := Sprintf("%o_%o", label, source_index);

    print "SOURCE", source_label, "family", label, "param", param, "t", t,
          "torsion", source_inv, "order", source_order, "exponent", source_exp,
          "factor_type", ftype;
    print "SOURCE_CURVE", source_f_int;

    if source_order lt min_source_order then
        print "SOURCE_SKIP_SMALL", source_label;
        return;
    end if;

    ConsiderObject(~global_best_order, ~global_best_exponent,
                   ~global_best_records, ~interesting_records,
                   source_label, "source", source_inv, source_f_int);

    J := Jacobian(C);

    try
        richelots := RichelotIsogenousSurfaces(J);
        print "RICHELOT_COUNT", source_label, #richelots;
        for i in [1..#richelots] do
            obj_label := Sprintf("richelot_%o", i);
            print "RICHELOT_TYPE", source_label, obj_label, Type(richelots[i]);
            if Type(richelots[i]) eq JacHyp then
                inv, f_int := TorsionInvariantsForCurve(Curve(richelots[i]));
                print "RICHELOT", source_label, obj_label, "torsion", inv,
                      "order", InvOrder(inv), "exponent", InvExponent(inv),
                      "factor_type", FactorTypeString(Factorization(f_int));
                ConsiderObject(~global_best_order, ~global_best_exponent,
                               ~global_best_records, ~interesting_records,
                               source_label, obj_label, inv, f_int);
            end if;
        end for;
    catch e
        print "RICHELOT_ERROR", source_label, e`Object;
    end try;

    try
        Js, products, weil_restrictions := TwoPowerIsogenies(J);
        print "TWOPOWER_COUNT", source_label, "jacobians", #Js,
              "products", #products, "weil_restrictions", #weil_restrictions;
        for i in [1..#Js] do
            obj_label := Sprintf("twopower_%o", i);
            inv, f_int := TorsionInvariantsForCurve(Curve(Js[i]));
            print "TWOPOWER", source_label, obj_label, "torsion", inv,
                  "order", InvOrder(inv), "exponent", InvExponent(inv),
                  "factor_type", FactorTypeString(Factorization(f_int));
            ConsiderObject(~global_best_order, ~global_best_exponent,
                           ~global_best_records, ~interesting_records,
                           source_label, obj_label, inv, f_int);
        end for;
        for i in [1..#products] do
            print "TWOPOWER_PRODUCT", source_label, i, products[i];
        end for;
        for i in [1..#weil_restrictions] do
            print "TWOPOWER_WEIL_RESTRICTION", source_label, i, weil_restrictions[i];
        end for;
    catch e
        print "TWOPOWER_ERROR", source_label, e`Object;
    end try;
end procedure;

params := RationalParametersOfHeight(height);
seen_t := {};
checked := AssociativeArray();
tested_sources := AssociativeArray();
accepted_sources := AssociativeArray();
global_best_order := 0;
global_best_exponent := 0;
global_best_records := [];
interesting_records := [];

for label in ["linear", "qq"] do
    checked[label] := 0;
    tested_sources[label] := 0;
    accepted_sources[label] := 0;
end for;

print "Richelot sweep over contact-5 extra-2 [2,20] sources";
print "height", height, "params", #params,
      "max_sources_per_label", max_sources_per_label,
      "min_source_order", min_source_order, "target_order", target_order;

for label in ["linear", "qq"] do
    for z in params do
        if tested_sources[label] ge max_sources_per_label then
            continue;
        end if;

        checked[label] +:= 1;
        ok, t := TFromParameter(label, z);
        if not ok or t eq -Q!1 then
            continue;
        end if;

        tkey := Sprint(t);
        if tkey in seen_t then
            continue;
        end if;
        Include(~seen_t, tkey);

        tested_sources[label] +:= 1;
        before := #interesting_records;
        SweepSource(label, z, t, tested_sources[label],
                    ~global_best_order, ~global_best_exponent,
                    ~global_best_records, ~interesting_records);
        accepted_sources[label] +:= 1;
        if #interesting_records gt before then
            print "NEW_INTERESTING_AFTER_SOURCE", label, tested_sources[label];
        end if;
    end for;
end for;

print "DONE";
print "checked_linear", checked["linear"], "checked_qq", checked["qq"],
      "tested_linear", tested_sources["linear"], "tested_qq", tested_sources["qq"],
      "unique_t", #seen_t;
print "best_order", global_best_order, "best_exponent", global_best_exponent;
print "best_records", global_best_records;
print "interesting_count", #interesting_records;
for rec in interesting_records do
    print "INTERESTING_RECORD", rec[1], rec[2], "torsion", rec[3],
          "order", rec[4], "exponent", rec[5];
end for;

quit;
