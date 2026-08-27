//////////////////////////////////////////////////////////////////////
//  M(2,2,2,8) plus rational 3-torsion.
//
//  For C: y^2 = f(x), deg(f)=5, a nonzero rational 3-torsion
//  class represented by a degree-2 Mumford divisor is equivalent to
//  a cubic h having triple contact with C along that divisor:
//
//      h(x)^2 - f(x) = m^2 q(x)^3,
//
//  where
//
//      q = x^2 + U*x + V,       h = m*x^3 + N*x^2 + R*x + S.
//
//  For the M(2,2,2,8) model
//
//      f = x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2),
//
//  this gives explicit algebraic equations over the K3 surface
//  s2(a,b,c,d)^2 = 4abcd.  The search below uses this condition as
//  motivation, but first applies the cheap necessary condition that
//  3 divides #J(F_p) at every good prime p.
//
//  Typical runs from torsion_jac:
//
//      magma -b height:=25 mode:="named" code/m2228_three_torsion_search.m
//
//      magma -b height:=20 max_multiple:=4 mode:="multiples" \
//          code/m2228_three_torsion_search.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 25;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

if not assigned max_multiple then
    max_multiple := 4;
elif Type(max_multiple) eq MonStgElt then
    max_multiple := StringToInteger(max_multiple);
end if;

if not assigned mode then
    mode := "named";
end if;

if not assigned progress_interval then
    progress_interval := 0;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

load "code/m2248_section_multiples.m";

Q := Rationals();
PX<X> := PolynomialRing(Q);

prime_list := [
    5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
    73,79,83,89,97,101,103,107,109,113,127,131,137,139,
    149,151,157,163,167,173,179,181,191,193,197,199
];

procedure PrintThreeTorsionContactCondition()
    print "3-torsion cubic-contact condition";
    print "  f = x^5 + e1*x^4 + e2*x^3 + e3*x^2 + e4*x";
    print "  q = x^2 + U*x + V";
    print "  h = m*x^3 + N*x^2 + R*x + S";
    print "  h^2 - f = m^2*q^3";
    print "Coefficient equations:";
    print "  2*m*N - 1 = 3*m^2*U";
    print "  N^2 + 2*m*R - e1 = 3*m^2*(U^2 + V)";
    print "  2*m*S + 2*N*R - e2 = m^2*(U^3 + 6*U*V)";
    print "  R^2 + 2*N*S - e3 = 3*m^2*(U^2*V + V^2)";
    print "  2*R*S - e4 = 3*m^2*U*V^2";
    print "  S^2 = m^2*V^3";
    print "Nondegenerate conditions: m != 0, U^2 - 4V != 0, gcd(q,f)=1.";
end procedure;

function CurvePolynomialFromTuple(vals)
    f := X;
    for v in vals do
        f *:= X + (Q!v)^2;
    end for;
    return f;
end function;

function TorsionPrime3Filter(f, primes)
    C := HyperellipticCurve(f);
    used := [];
    gcd_bound := 0;

    for p in primes do
        try
            fp := ChangeRing(f, GF(p));
        catch e
            continue;
        end try;

        if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
            continue;
        end if;

        try
            Cp := ChangeRing(C, GF(p));
            n := Integers()!#Jacobian(Cp);
        catch e
            continue;
        end try;

        if gcd_bound eq 0 then
            gcd_bound := n;
        else
            gcd_bound := GCD(gcd_bound, n);
        end if;
        Append(~used, <p, n, gcd_bound>);

        if n mod 3 ne 0 then
            return false, p, n, gcd_bound, used;
        end if;
    end for;

    return true, 0, 0, gcd_bound, used;
end function;

procedure RecordCandidate(~candidates, source, param, tuple, gcd_bound, used)
    Append(~candidates, <source, param, tuple, gcd_bound, used>);
    print "  candidate", source, "param", param, "gcd_bound", gcd_bound;
    print "    tuple", tuple;
    print "    used", used;
end procedure;

procedure SearchNamedSectionFamilies(height)
    params := RationalParametersOfHeight(height);
    total_all := 0;
    usable_all := 0;
    candidates := [];

    print "Searching named M(2,2,2,8) section families";
    print "height", height, "parameters", #params;
    print "primes", prime_list;

    for family in M2248SectionFamilyNames do
        total := 0;
        usable := 0;
        obstructed := 0;
        bad_or_unusable := 0;
        candidates_before := #candidates;

        for idx in [1..#params] do
            t := params[idx];
            total +:= 1;
            total_all +:= 1;

            vals := [];
            ok := true;
            try
                vals := SectionTuple(family, t);
            catch e
                ok := false;
            end try;

            if not ok or not IsUsableTuple(vals) then
                bad_or_unusable +:= 1;
                continue;
            end if;

            tuple := PrimitiveIntegerTuple(vals);
            f := CurvePolynomialFromTuple(tuple);
            if Degree(f) ne 5 or Discriminant(f) eq 0 then
                bad_or_unusable +:= 1;
                continue;
            end if;

            usable +:= 1;
            usable_all +:= 1;
            survives, pbad, nbad, gcd_bound, used := TorsionPrime3Filter(f, prime_list);
            if survives then
                RecordCandidate(~candidates, family, t, tuple, gcd_bound, used);
            else
                obstructed +:= 1;
            end if;

            if progress_interval gt 0 and total_all mod progress_interval eq 0 then
                print "progress", total_all, "usable", usable_all,
                      "candidates", #candidates;
            end if;
        end for;

        print family, "total", total, "usable", usable,
              "bad_or_unusable", bad_or_unusable,
              "obstructed", obstructed,
              "survived_3_filter", #candidates - candidates_before;
    end for;

    print "Named-family search done";
    print "total", total_all, "usable", usable_all, "candidates", #candidates;
end procedure;

procedure SearchSectionMultiples(height, max_multiple)
    params := RationalParametersOfHeight(height);
    E, P0, torsion, T := ATEllipticData();
    total_all := 0;
    usable_all := 0;
    candidates := [];

    print "Searching M(2,2,2,8) section multiples";
    print "height", height, "parameters", #params, "max_multiple", max_multiple;
    print "primes", prime_list;

    for n in [1..max_multiple] do
        for tors_data in torsion do
            label := tors_data[1];
            ok_section, vals_generic := ATSectionTuple(n, label);
            if not ok_section then
                print "section", n, label, "skipped generically";
                continue;
            end if;

            source := Sprintf("n=%o,torsion=%o", n, label);
            total := 0;
            usable := 0;
            obstructed := 0;
            bad_or_unusable := 0;

            for idx in [1..#params] do
                t := params[idx];
                total +:= 1;
                total_all +:= 1;

                ok_spec, vals := SpecializeTupleAtQ(vals_generic, t);
                if not ok_spec or not IsUsableTuple(vals) then
                    bad_or_unusable +:= 1;
                    continue;
                end if;

                tuple := PrimitiveIntegerTuple(vals);
                f := CurvePolynomialFromTuple(tuple);
                if Degree(f) ne 5 or Discriminant(f) eq 0 then
                    bad_or_unusable +:= 1;
                    continue;
                end if;

                usable +:= 1;
                usable_all +:= 1;
                survives, pbad, nbad, gcd_bound, used := TorsionPrime3Filter(f, prime_list);
                if survives then
                    RecordCandidate(~candidates, source, t, tuple, gcd_bound, used);
                else
                    obstructed +:= 1;
                end if;

                if progress_interval gt 0 and total_all mod progress_interval eq 0 then
                    print "progress", total_all, "usable", usable_all,
                          "candidates", #candidates;
                end if;
            end for;

            print "section", n, label, "total", total, "usable", usable,
                  "bad_or_unusable", bad_or_unusable,
                  "obstructed", obstructed,
                  "survived_3_filter_so_far", #candidates;
        end for;
    end for;

    print "Section-multiple search done";
    print "total", total_all, "usable", usable_all, "candidates", #candidates;
end procedure;

PrintThreeTorsionContactCondition();

if mode eq "named" then
    SearchNamedSectionFamilies(height);
elif mode eq "multiples" then
    SearchSectionMultiples(height, max_multiple);
else
    error "Unknown mode. Use mode:=\"named\" or mode:=\"multiples\".";
end if;

quit;
