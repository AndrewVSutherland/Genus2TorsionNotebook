//////////////////////////////////////////////////////////////////////
//  Full rational 2-torsion, i.e. M(2,2,2,2), plus possible 7-torsion.
//
//  We use the normalized full-split model
//
//      C: y^2 = x*(x-1)*(x-a)*(x-b)*(x-c),
//
//  with distinct a,b,c not in {0,1}.  This gives full rational
//  2-torsion on the Jacobian.  A necessary condition for rational
//  7-torsion is that 7 divides #J(F_p) for every good prime p != 7.
//
//  Modes:
//      finite      finite-field density check on the full-split chart
//      search      rational height search with reduction filter
//      prime_diag  prime-by-prime diagnostic for the rational search box
//
//  Typical runs:
//      magma -b mode:="finite" code/m2222_plus7_search.m
//      magma -b mode:="search" height:=8 code/m2222_plus7_search.m
//      magma -b mode:="prime_diag" height:=8 code/m2222_plus7_search.m
//////////////////////////////////////////////////////////////////////

if not assigned mode then
    mode := "finite";
end if;
if not assigned height then
    height := 8;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;
if not assigned max_tests then
    max_tests := 200;
elif Type(max_tests) eq MonStgElt then
    max_tests := StringToInteger(max_tests);
end if;

Q := Rationals();
Z := Integers();
Qx<x> := PolynomialRing(Q);

finite_primes := [5,11,13,17,19,23,29,31,37,41,43,47];
filter_primes := [5,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71];

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            r := Q!num/Q!den;
            key := Sprint(r);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, r);
            end if;
        end for;
    end for;
    return vals;
end function;

function FullSplitPolynomial(a, b, c)
    return x*(x-1)*(x-a)*(x-b)*(x-c);
end function;

function GoodHyperellipticPolynomial(f)
    return Degree(f) in {5,6} and Discriminant(f) ne 0;
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return Qx!(L^2*f), L;
end function;

function HasSevenInvariants(invs)
    return &or [ n mod 7 eq 0 : n in invs ];
end function;

function PassesSevenReduction(f, primes)
    for p in primes do
        try
            fp := ChangeRing(f, GF(p));
        catch e
            continue;
        end try;
        if not GoodHyperellipticPolynomial(fp) then
            continue;
        end if;
        C := HyperellipticCurve(fp);
        if (#Jacobian(C) mod 7) ne 0 then
            return false, p;
        end if;
    end for;
    return true, 0;
end function;

function TupleKey(vals)
    return Sprint(Sort([Q!v : v in vals]));
end function;

procedure FiniteM2222Plus7()
    print "FINITE M(2,2,2,2) plus 7 necessary condition";
    print "finite_primes", finite_primes;
    for p in finite_primes do
        K := GF(p);
        P<X> := PolynomialRing(K);
        vals := [u : u in K | u ne 0 and u ne 1];
        total := 0;
        good := 0;
        seven := 0;
        samples := [];
        for i in [1..#vals] do
            for j in [i+1..#vals] do
                for k in [j+1..#vals] do
                    a := vals[i]; b := vals[j]; c := vals[k];
                    total +:= 1;
                    f := X*(X-1)*(X-a)*(X-b)*(X-c);
                    if not GoodHyperellipticPolynomial(f) then
                        continue;
                    end if;
                    good +:= 1;
                    C := HyperellipticCurve(f);
                    if (#Jacobian(C) mod 7) eq 0 then
                        seven +:= 1;
                        if #samples lt 8 then
                            Append(~samples, <Z!a,Z!b,Z!c>);
                        end if;
                    end if;
                end for;
            end for;
        end for;
        print "p", p, "total", total, "good", good,
              "seven_possible", seven, "samples", samples;
    end for;
end procedure;

procedure SearchM2222Plus7()
    print "SEARCH M(2,2,2,2) plus 7";
    print "height", height, "filter_primes", filter_primes;
    params0 := RationalParametersOfHeight(height);
    params := [u : u in params0 | u ne 0 and u ne 1];
    print "params", #params;

    checked := 0;
    smooth := 0;
    reduction_survivors := 0;
    exact_tests := 0;
    hits := [];

    for i in [1..#params] do
        if #hits ge max_hits then break; end if;
        for j in [i+1..#params] do
            if #hits ge max_hits then break; end if;
            for k in [j+1..#params] do
                a := params[i]; b := params[j]; c := params[k];
                if #(Set([a,b,c,Q!0,Q!1])) ne 5 then
                    continue;
                end if;
                checked +:= 1;
                f := FullSplitPolynomial(a,b,c);
                if not GoodHyperellipticPolynomial(f) then
                    continue;
                end if;
                smooth +:= 1;
                ok, badp := PassesSevenReduction(f, filter_primes);
                if not ok then
                    continue;
                end if;
                reduction_survivors +:= 1;
                fI, L := IntegralModelPolynomial(f);
                C := HyperellipticCurve(fI);
                J := Jacobian(C);
                G, phi := TorsionSubgroup(J);
                invs := Invariants(G);
                exact_tests +:= 1;
                print "REDUCTION_SURVIVOR", [a,b,c], "torsion", invs, "f", fI;
                if HasSevenInvariants(invs) then
                    Append(~hits, <a,b,c,invs,fI>);
                    print "HIT", [a,b,c], "torsion", invs, "f", fI;
                    if #hits ge max_hits then
                        break;
                    end if;
                end if;
                if exact_tests ge max_tests then
                    break i;
                end if;
            end for;
        end for;
    end for;

    print "DONE m2222_plus7 height", height,
          "checked", checked,
          "smooth", smooth,
          "reduction_survivors", reduction_survivors,
          "exact_tests", exact_tests,
          "hits", #hits;
end procedure;

procedure PrimeDiagnostic()
    print "PRIME_DIAG M(2,2,2,2) plus 7";
    print "height", height, "filter_primes", filter_primes;
    params0 := RationalParametersOfHeight(height);
    params := [u : u in params0 | u ne 0 and u ne 1];

    checked := 0;
    smooth := 0;
    pass_all := 0;

    good_counts := AssociativeArray(Integers());
    bad_counts := AssociativeArray(Integers());
    pass_counts := AssociativeArray(Integers());
    kill_counts := AssociativeArray(Integers());
    first_kill_counts := AssociativeArray(Integers());
    for p in filter_primes do
        good_counts[p] := 0;
        bad_counts[p] := 0;
        pass_counts[p] := 0;
        kill_counts[p] := 0;
        first_kill_counts[p] := 0;
    end for;

    for i in [1..#params] do
        for j in [i+1..#params] do
            for k in [j+1..#params] do
                a := params[i]; b := params[j]; c := params[k];
                if #(Set([a,b,c,Q!0,Q!1])) ne 5 then
                    continue;
                end if;
                checked +:= 1;
                f := FullSplitPolynomial(a,b,c);
                if not GoodHyperellipticPolynomial(f) then
                    continue;
                end if;
                smooth +:= 1;
                first_fail := 0;
                for p in filter_primes do
                    try
                        fp := ChangeRing(f, GF(p));
                    catch e
                        bad_counts[p] +:= 1;
                        continue;
                    end try;
                    if not GoodHyperellipticPolynomial(fp) then
                        bad_counts[p] +:= 1;
                        continue;
                    end if;
                    good_counts[p] +:= 1;
                    C := HyperellipticCurve(fp);
                    if (#Jacobian(C) mod 7) eq 0 then
                        pass_counts[p] +:= 1;
                    else
                        kill_counts[p] +:= 1;
                        if first_fail eq 0 then
                            first_fail := p;
                        end if;
                    end if;
                end for;
                if first_fail eq 0 then
                    pass_all +:= 1;
                else
                    first_kill_counts[first_fail] +:= 1;
                end if;
            end for;
        end for;
    end for;

    print "checked", checked;
    print "smooth", smooth;
    print "pass_all", pass_all;
    for p in filter_primes do
        print "PRIME", p,
              "good", good_counts[p],
              "bad_or_boundary", bad_counts[p],
              "pass7", pass_counts[p],
              "kill7", kill_counts[p],
              "first_kill", first_kill_counts[p];
    end for;
end procedure;

if mode eq "finite" then
    FiniteM2222Plus7();
elif mode eq "search" then
    SearchM2222Plus7();
elif mode eq "prime_diag" then
    PrimeDiagnostic();
else
    error "Unknown mode";
end if;
