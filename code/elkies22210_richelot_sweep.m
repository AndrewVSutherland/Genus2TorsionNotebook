//////////////////////////////////////////////////////////////////////
//  Elkies [2,2,2,10] family: source and Richelot/2-power sweep.
//
//  From Elkies, "Families of genus-2 curves with 5-torsion",
//  Proposition 2.11: full level-2 atypical 5-torsion is parametrized by
//  points (r1:...:r5) on the Clebsch-Klein cubic
//
//      sum_i r_i = 0,     sum_i r_i^3 = 0.
//
//  For a primitive integral representative with nonzero distinct squares,
//  the corresponding split model is
//
//      C: y^2 = x * prod_i (x - r_i^2).
//
//  This gives full rational 2-torsion and a rational 5-torsion point,
//  so the generic exact torsion is [2,2,2,10].  We search bounded
//  integral projective points and test whether any source or rational
//  Richelot / 2-power isogenous Jacobian has larger rational torsion.
//
//  Typical runs:
//      magma -b height:=12 max_sources:=20 code/elkies22210_richelot_sweep.m
//      magma -b height:=25 max_sources:=100 run_isogenies:=false \
//          code/elkies22210_richelot_sweep.m
//      magma -b height:=30 max_sources:=100000 run_isogenies:=false \
//          run_source_halving:=true code/elkies22210_richelot_sweep.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 12;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

if not assigned max_sources then
    max_sources := 25;
elif Type(max_sources) eq MonStgElt then
    max_sources := StringToInteger(max_sources);
end if;

if not assigned run_isogenies then
    run_isogenies := true;
elif Type(run_isogenies) eq MonStgElt then
    run_isogenies := run_isogenies in {"true", "True", "1", "yes"};
end if;

if not assigned max_twopower_sources then
    max_twopower_sources := max_sources;
elif Type(max_twopower_sources) eq MonStgElt then
    max_twopower_sources := StringToInteger(max_twopower_sources);
end if;

// Independently audit whether any of the 15 nonzero rational J[2]
// classes on the full-split source is divisible by 2 over Q.  This is
// opt-in so that the historical source/Richelot sweep has unchanged
// runtime and output by default.
if not assigned run_source_halving then
    run_source_halving := false;
elif Type(run_source_halving) eq MonStgElt then
    run_source_halving := run_source_halving in {"true", "True", "1", "yes"};
end if;

if not assigned source_halving_verbose then
    source_halving_verbose := true;
elif Type(source_halving_verbose) eq MonStgElt then
    source_halving_verbose := source_halving_verbose in
        {"true", "True", "1", "yes"};
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

function AbsGCD(vals)
    nz := [ Abs(Z!v) : v in vals | v ne 0 ];
    if #nz eq 0 then
        return 0;
    end if;
    return GCD(nz);
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

function FactorTypeString(f)
    fac := Factorization(f);
    degs := Sort([ Degree(ff[1]) : ff in fac ]);
    return Join([ IntegerToString(d) : d in degs ], "+");
end function;

function TorsionInvariantsForCurve(C)
    J := Jacobian(C);
    A, phi := TorsionSubgroup(J);
    return Invariants(A);
end function;

function TorsionInvariantsForJacobian(J)
    A, phi := TorsionSubgroup(J);
    return Invariants(A);
end function;

function SourcePolynomial(rs)
    return x * &*[ x - (Q!r)^2 : r in rs ];
end function;

function SourceTwoTorsionClasses(J, rs)
    // The even-degree model has the six finite Weierstrass roots
    // 0,r1^2,...,r5^2.  Its 15 nonzero J[2] classes are indexed by
    // unordered pairs of these roots and represented by [u,0], where
    // u=(x-beta_i)(x-beta_j).
    betas := [ Q!0 ] cat [ (Q!r)^2 : r in rs ];
    out := [];
    for i in [1..5] do
        for j in [i+1..6] do
            u := (x-betas[i])*(x-betas[j]);
            T := J![u, Q!0];
            Append(~out, <i-1, j-1, betas[i], betas[j], u, T>);
        end for;
    end for;
    assert #out eq 15;
    for i in [1..14] do
        for j in [i+1..15] do
            assert out[i][6] ne out[j][6];
        end for;
    end for;
    return out;
end function;

procedure AuditSourceHalving(J, rs, source_label, source_inv,
                              ~sources_audited, ~classes_audited,
                              ~divisible_classes, ~halving_records)
    classes := SourceTwoTorsionClasses(J, rs);
    ZJ := J!0;
    sources_audited +:= 1;
    print "SOURCE_HALVING_START", source_label,
          "classes", #classes,
          "source_torsion", source_inv;

    for rec in classes do
        i, j, beta_i, beta_j, u, T := Explode(rec);
        assert T ne ZJ and 2*T eq ZJ;
        classes_audited +:= 1;

        divisible, half := IsDivisibleBy(T, 2);
        if source_halving_verbose then
            print "SOURCE_HALVING_CLASS", source_label,
                  "pair", [i,j], "roots", [beta_i,beta_j],
                  "u", u, "divisible", divisible;
        end if;
        if not divisible then
            continue;
        end if;

        assert 2*half eq T;
        half_order := Order(half);
        assert half_order eq 4;
        divisible_classes +:= 1;

        // TorsionSubgroup was already computed exactly for the source.
        // Full rational J[2], the marked 5-part, and this order-4 half
        // force order at least 2^5*5=160.
        assert InvOrder(source_inv) ge 160;
        Append(~halving_records,
               <source_label, i, j, beta_i, beta_j, half_order, source_inv>);
        print "SOURCE_HALVING_HIT", source_label,
              "pair", [i,j], "roots", [beta_i,beta_j],
              "half_order", half_order,
              "certified_torsion", source_inv,
              "certified_order", InvOrder(source_inv);
    end for;

    print "SOURCE_HALVING_END", source_label,
          "classes", #classes,
          "divisible", #[ rec : rec in halving_records |
              rec[1] eq source_label ];
end procedure;

function SquareKey(rs)
    return Join([ IntegerToString(s) : s in Sort([ Z!(r^2) : r in rs ]) ], ",");
end function;

function IsCanonicalTuple(rs)
    if AbsGCD(rs) ne 1 then
        return false;
    end if;
    if &or [ r eq 0 : r in rs ] then
        return false;
    end if;
    sqs := [ r^2 : r in rs ];
    if #Seqset(sqs) ne #sqs then
        return false;
    end if;
    // Fix the global sign ambiguity.
    for r in rs do
        if r ne 0 then
            return r gt 0;
        end if;
    end for;
    return false;
end function;

function IrreducibleLPolynomialCertificate(f)
    C := HyperellipticCurve(f);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79] do
        try
            fp := ChangeRing(f, GF(p));
            if Degree(fp) notin {5,6} or Discriminant(fp) eq 0 then
                continue;
            end if;
            Cp := ChangeRing(C, GF(p));
            Lp := LPolynomial(Cp);
            fac := Factorization(Lp);
            if #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4 then
                return true, p, Lp;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false, 0, P!0;
end function;

procedure Consider(~best_order, ~best_exponent, ~best_records,
                   ~interesting_records, source_label, object_label, inv, f)
    ord := InvOrder(inv);
    exp := InvExponent(inv);
    if ord gt best_order then
        best_order := ord;
        best_exponent := exp;
        best_records := [ <source_label, object_label, inv, ord, exp> ];
    elif ord eq best_order then
        if exp gt best_exponent then
            best_exponent := exp;
        end if;
        Append(~best_records, <source_label, object_label, inv, ord, exp>);
    end if;

    if ord gt 80 or exp gt 10 then
        Append(~interesting_records, <source_label, object_label, inv, ord, exp, f>);
        print "INTERESTING", source_label, object_label, "torsion", inv,
              "order", ord, "exponent", exp;
        print "INTERESTING_CURVE", f;
    end if;
end procedure;

procedure SweepSource(rs, source_index, ~twopower_sources_done,
                      ~best_order, ~best_exponent, ~best_records,
                      ~interesting_records,
                      ~source_halving_sources,
                      ~source_halving_classes,
                      ~source_halving_divisible,
                      ~source_halving_records)
    label := Sprintf("source_%o", source_index);
    f := SourcePolynomial(rs);
    C := HyperellipticCurve(f);
    inv := TorsionInvariantsForCurve(C);
    simple, pcert, Lp := IrreducibleLPolynomialCertificate(f);

    print "SOURCE", label, "rs", rs, "squares", Sort([Z!(r^2) : r in rs]),
          "torsion", inv, "order", InvOrder(inv), "exponent", InvExponent(inv),
          "factor_type", FactorTypeString(f),
          "simple_cert", simple, "pcert", pcert;
    if simple then
        print "SOURCE_LPOLY", label, Lp;
    end if;
    print "SOURCE_CURVE", label, f;

    Consider(~best_order, ~best_exponent, ~best_records, ~interesting_records,
             label, "source", inv, f);

    if run_source_halving then
        AuditSourceHalving(Jacobian(C), rs, label, inv,
                           ~source_halving_sources,
                           ~source_halving_classes,
                           ~source_halving_divisible,
                           ~source_halving_records);
    end if;

    if not run_isogenies then
        return;
    end if;

    J := Jacobian(C);

    try
        richelots := RichelotIsogenousSurfaces(J);
        print "RICHELOT_COUNT", label, #richelots;
        for i in [1..#richelots] do
            obj := Sprintf("richelot_%o", i);
            print "RICHELOT_TYPE", label, obj, Type(richelots[i]);
            if Type(richelots[i]) eq JacHyp then
                D := Curve(richelots[i]);
                invR := TorsionInvariantsForCurve(D);
                fR, hR := HyperellipticPolynomials(D);
                print "RICHELOT", label, obj, "torsion", invR,
                      "order", InvOrder(invR), "exponent", InvExponent(invR);
                Consider(~best_order, ~best_exponent, ~best_records,
                         ~interesting_records, label, obj, invR, P!fR);
            end if;
        end for;
    catch e
        print "RICHELOT_ERROR", label, e`Object;
    end try;

    if twopower_sources_done ge max_twopower_sources then
        print "TWOPOWER_SKIP_LIMIT", label;
        return;
    end if;
    twopower_sources_done +:= 1;

    try
        Js, products, weil_restrictions := TwoPowerIsogenies(J);
        print "TWOPOWER_COUNT", label, "jacobians", #Js,
              "products", #products, "weil_restrictions", #weil_restrictions;
        for i in [1..#Js] do
            obj := Sprintf("twopower_%o", i);
            invT := TorsionInvariantsForJacobian(Js[i]);
            fT, hT := HyperellipticPolynomials(Curve(Js[i]));
            print "TWOPOWER", label, obj, "torsion", invT,
                  "order", InvOrder(invT), "exponent", InvExponent(invT);
            Consider(~best_order, ~best_exponent, ~best_records,
                     ~interesting_records, label, obj, invT, P!fT);
        end for;
        for i in [1..#products] do
            print "TWOPOWER_PRODUCT", label, i, products[i];
        end for;
        for i in [1..#weil_restrictions] do
            print "TWOPOWER_WEIL_RESTRICTION", label, i, weil_restrictions[i];
        end for;
    catch e
        print "TWOPOWER_ERROR", label, e`Object;
    end try;
end procedure;

print "Elkies [2,2,2,10] Clebsch-Klein source/Richelot sweep";
print "height", height, "max_sources", max_sources,
      "run_isogenies", run_isogenies, "max_twopower_sources", max_twopower_sources;
if run_source_halving then
    print "SOURCE_HALVING_CONFIG", "enabled", true,
          "verbose", source_halving_verbose;
end if;

seen_square_keys := {};
source_count := 0;
tuple_checked := 0;
clebsch_points := 0;
best_order := 0;
best_exponent := 0;
best_records := [];
interesting_records := [];
twopower_sources_done := 0;
triple_checked := 0;
source_halving_sources := 0;
source_halving_classes := 0;
source_halving_divisible := 0;
source_halving_records := [];

for r1 in [-height..height] do
    for r2 in [-height..height] do
        for r3 in [-height..height] do
            triple_checked +:= 1;
            s := r1 + r2 + r3;
            if s eq 0 then
                continue;
            end if;
            A := r1^3 + r2^3 + r3^3;

            // With r5 = -(r1+r2+r3+r4), the cubic condition is
            // 3*s*r4^2 + 3*s^2*r4 + (s^3-A) = 0.
            disc := 3*s*(4*A - s^3);
            is_square, root := IsSquare(disc);
            if not is_square then
                continue;
            end if;

            for root_signed in Sort(Setseq({-root, root})) do
                q4 := (Q!(-3*s^2 + root_signed))/(Q!(6*s));
                if Denominator(q4) ne 1 then
                    continue;
                end if;
                r4 := Z!q4;
                r5 := -(r1+r2+r3+r4);
                if Abs(r4) gt height or Abs(r5) gt height then
                    continue;
                end if;
                rs := [r1,r2,r3,r4,r5];
                tuple_checked +:= 1;
                if not IsCanonicalTuple(rs) then
                    continue;
                end if;
                clebsch_points +:= 1;
                key := SquareKey(rs);
                if key in seen_square_keys then
                    continue;
                end if;
                Include(~seen_square_keys, key);

                source_count +:= 1;
                SweepSource(rs, source_count, ~twopower_sources_done,
                            ~best_order, ~best_exponent, ~best_records,
                            ~interesting_records,
                            ~source_halving_sources,
                            ~source_halving_classes,
                            ~source_halving_divisible,
                            ~source_halving_records);
                if source_count ge max_sources then
                    print "SOURCE_LIMIT_REACHED";
                    print "DONE_PARTIAL";
                    print "triple_checked", triple_checked,
                          "tuple_checked", tuple_checked,
                          "clebsch_points", clebsch_points,
                          "unique_sources", source_count;
                    print "best_order", best_order, "best_exponent", best_exponent;
                    print "best_records", best_records;
                    print "interesting_count", #interesting_records;
                    if run_source_halving then
                        print "SOURCE_HALVING_SUMMARY",
                              "sources", source_halving_sources,
                              "classes", source_halving_classes,
                              "divisible", source_halving_divisible,
                              "records", source_halving_records;
                    end if;
                    for rec in interesting_records do
                        print "INTERESTING_RECORD", rec[1], rec[2],
                              "torsion", rec[3], "order", rec[4],
                              "exponent", rec[5];
                    end for;
                    quit;
                end if;
            end for;
        end for;
    end for;
end for;

print "DONE";
print "triple_checked", triple_checked,
      "tuple_checked", tuple_checked,
      "clebsch_points", clebsch_points,
      "unique_sources", source_count;
print "best_order", best_order, "best_exponent", best_exponent;
print "best_records", best_records;
print "interesting_count", #interesting_records;
if run_source_halving then
    print "SOURCE_HALVING_SUMMARY",
          "sources", source_halving_sources,
          "classes", source_halving_classes,
          "divisible", source_halving_divisible,
          "records", source_halving_records;
end if;
for rec in interesting_records do
    print "INTERESTING_RECORD", rec[1], rec[2],
          "torsion", rec[3], "order", rec[4], "exponent", rec[5];
end for;

quit;
