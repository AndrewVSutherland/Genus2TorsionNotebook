//////////////////////////////////////////////////////////////////////
//  M(2,2,2,4) plus possible rational 3-torsion.
//
//  Base family:
//      C: y^2 = x (x+a^2) (x+b^2) (x+c^2) (x+d^2).
//
//  The visible full rational 2-torsion and the half of (0,0)-infinity
//  give torsion containing [2,2,2,4] on the nonsingular open chart.
//  A rational 3-torsion class is tested using the same cubic-contact
//  equations as in m2228_three_torsion_surface_sieve.m.
//
//  Modes:
//    finite:
//      finite-field good-chart diagnostic for the cubic-contact cover.
//
//    search:
//      scan primitive positive sorted tuples 0 < a < b < c < d <= height.
//      By default only tuples bad at p=13 are kept, since the good chart
//      has no 3-contact points at p=13.
//
//    exact_file:
//      read tuple_file containing rows [a,b,c,d] and apply the same filters.
//
//  Typical runs:
//      magma -b mode:="finite" code/m2224_plus3_search.m
//
//      magma -b mode:="search" height:=80 prime_bound:=37 \
//          point_prime_bound:=73 max_exact:=20 code/m2224_plus3_search.m
//
//      magma -b mode:="exact_file" tuple_file:="data/file.txt" \
//          code/m2224_plus3_search.m
//////////////////////////////////////////////////////////////////////

if not assigned mode then
    mode := "finite";
end if;

if not assigned height then
    height := 50;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

if not assigned prime_bound then
    prime_bound := 37;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;

if not assigned point_prime_bound then
    point_prime_bound := 73;
elif Type(point_prime_bound) eq MonStgElt then
    point_prime_bound := StringToInteger(point_prime_bound);
end if;

if not assigned max_exact then
    max_exact := 50;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;

if not assigned max_print then
    max_print := 20;
elif Type(max_print) eq MonStgElt then
    max_print := StringToInteger(max_print);
end if;

if not assigned progress_interval then
    progress_interval := 100000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

if not assigned boundary13_only then
    boundary13_only := true;
elif Type(boundary13_only) eq MonStgElt then
    boundary13_only := boundary13_only in {"true", "True", "1", "yes"};
end if;

if not assigned boundary13_live_only then
    boundary13_live_only := false;
elif Type(boundary13_live_only) eq MonStgElt then
    boundary13_live_only := boundary13_live_only in {"true", "True", "1", "yes"};
end if;

if assigned tuple_file then
    exact_tuple_file := tuple_file;
else
    exact_tuple_file := "data/m2224_plus3_candidates.txt";
end if;

Q := Rationals();
Z := Integers();
PX<X> := PolynomialRing(Q);

finite_primes := [5,7,11,13,17,19,23,29,31,37];
filter_primes := [p : p in PrimesUpTo(prime_bound) | p notin {2,3,5,7,13}];
point_primes := [
    p : p in PrimesUpTo(point_prime_bound)
    | p notin {2,3,5,7,13} and p notin filter_primes
];
bound_primes := [
    11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,
    79,83,89,97,101,103,107,109,113,127,131,137,139,149,
    151,157,163,167,173,179,181,191,193,197,199
];

function ContactCoefficientKeys(p)
    F := GF(p);
    PF<x> := PolynomialRing(F);
    keys := {};
    witnesses := AssociativeArray();
    raw_count := 0;

    for L in F do
        if L eq 0 then
            continue;
        end if;
        M := L^2;

        for U in F do
            for v in F do
                if v eq 0 or U^2 - 4*v^2 eq 0 then
                    continue;
                end if;

                q := x^2 + U*x + v^2;
                for e1 in F do
                    A := 2*M*e1 + 6*(U^2+v^2) - (M+3*U)^2;
                    e2 := ((M+3*U)*A + 8*v^3 - 4*U^3 - 24*U*v^2)/(4*M);
                    e3 := (A^2 + 16*(M+3*U)*v^3
                           - 48*(U^2*v^2+v^4))/(16*M);
                    e4 := (A*v^3 - 6*U*v^4)/(2*M);
                    if e4 eq 0 then
                        continue;
                    end if;

                    f := x^5 + e1*x^4 + e2*x^3 + e3*x^2 + e4*x;
                    if Degree(GCD(q, f)) gt 0 then
                        continue;
                    end if;

                    key := <Z!e1, Z!e2, Z!e3, Z!e4>;
                    Include(~keys, key);
                    raw_count +:= 1;
                    if not IsDefined(witnesses, key) then
                        witnesses[key] := <Z!L, Z!U, Z!v>;
                    end if;
                end for;
            end for;
        end for;
    end for;

    return keys, witnesses, raw_count;
end function;

function CurveKeyFromSquares(sq)
    A := sq[1]; B := sq[2]; C := sq[3]; D := sq[4];
    e1 := A+B+C+D;
    e2 := A*B + A*C + A*D + B*C + B*D + C*D;
    e3 := A*B*C + A*B*D + A*C*D + B*C*D;
    e4 := A*B*C*D;
    return <Z!e1, Z!e2, Z!e3, Z!e4>;
end function;

function CurveKeyModp(tup, p)
    F := GF(p);
    squares := [(F!a)^2 : a in tup];
    if &*squares eq 0 or #(Set(squares)) lt #squares then
        return false, <0,0,0,0>;
    end if;
    return true, CurveKeyFromSquares(squares);
end function;

function IsBadAtPrime(tup, p)
    F := GF(p);
    squares := [];
    for a in tup do
        if F!a eq 0 then
            return true;
        end if;
        Append(~squares, (F!a)^2);
    end for;
    return #(Set(squares)) lt #squares;
end function;

function HasSquareCollisionAtPrime(tup, p)
    F := GF(p);
    squares := [(F!a)^2 : a in tup];
    return #(Set(squares)) lt #squares;
end function;

function CurvePolynomialFromTuple(tup)
    f := X;
    for a in tup do
        f *:= X + (Q!a)^2;
    end for;
    return f;
end function;

function ReadTupleFile(filename)
    S := Read(filename);
    rows := Split(S, "\n");
    tuples := [];

    for raw in rows do
        row := raw;
        if #row gt 0 and row[#row] eq "\r" then
            row := row[1..#row-1];
        end if;
        if #row ge 2 and row[1] eq "[" and row[#row] eq "]" then
            body := row[2..#row-1];
            parts := Split(body, ",");
            tup := [StringToInteger(part) : part in parts | #part gt 0];
            if #tup eq 4 then
                Append(~tuples, tup);
            end if;
        end if;
    end for;

    return tuples;
end function;

function TorsionGcdBound(f, primes)
    C := HyperellipticCurve(f);
    gcd_bound := 0;
    used := [];

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
            n := Z!#Jacobian(ChangeRing(C, GF(p)));
        catch e
            continue;
        end try;

        if gcd_bound eq 0 then
            gcd_bound := n;
        else
            gcd_bound := GCD(gcd_bound, n);
        end if;
        Append(~used, <p, n, gcd_bound>);

        if gcd_bound mod 3 ne 0 then
            return gcd_bound, used;
        end if;
    end for;

    return gcd_bound, used;
end function;

function IrreducibleFrobeniusCertificate(f)
    C := HyperellipticCurve(f);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79] do
        try
            fp := ChangeRing(f, GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
                continue;
            end if;
            Lp := LPolynomial(ChangeRing(C, GF(p)));
            fac := Factorization(Lp);
            if #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4 then
                return true, p, Lp;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false, 0, PX!0;
end function;

function PrecomputeContactData(primes)
    data := AssociativeArray();
    for p in primes do
        keys, witnesses, raw_count := ContactCoefficientKeys(p);
        data[p] := keys;
        print "prime", p, "contact_keys", #keys, "raw_contact_param", raw_count;
    end for;
    return data;
end function;

function PointCountAllowedKeys(p)
    F := GF(p);
    PF<x> := PolynomialRing(F);
    squares := Setseq({a^2 : a in F | a ne 0});
    allowed := {};
    total := 0;

    for i in [1..#squares-3] do
        for j in [i+1..#squares-2] do
            for k in [j+1..#squares-1] do
                for l in [k+1..#squares] do
                    branch := [squares[i], squares[j], squares[k], squares[l]];
                    key := CurveKeyFromSquares(branch);
                    f := x;
                    for a in branch do
                        f *:= x + a;
                    end for;
                    C := HyperellipticCurve(f);
                    n := Z!#Jacobian(C);
                    total +:= 1;
                    if n mod 3 eq 0 then
                        Include(~allowed, key);
                    end if;
                end for;
            end for;
        end for;
    end for;

    return allowed, total;
end function;

function PrecomputePointCountData(primes)
    data := AssociativeArray();
    for p in primes do
        allowed, total := PointCountAllowedKeys(p);
        data[p] := allowed;
        print "prime", p, "point_keys", total, "allowed_3div", #allowed;
    end for;
    return data;
end function;

function PassesContactFilters(tup, contact_data)
    for p in Keys(contact_data) do
        good, key := CurveKeyModp(tup, p);
        if good and not (key in contact_data[p]) then
            return false, p;
        end if;
    end for;
    return true, 0;
end function;

function PassesPointCountFilters(tup, point_count_data)
    for p in Keys(point_count_data) do
        good, key := CurveKeyModp(tup, p);
        if good and not (key in point_count_data[p]) then
            return false, p;
        end if;
    end for;
    return true, 0;
end function;

procedure ExactCheckTuple(tup, ~exact_tests, ~hits, ~printed)
    f := CurvePolynomialFromTuple(tup);
    if Discriminant(f) eq 0 then
        return;
    end if;

    C := HyperellipticCurve(f);
    J := Jacobian(C);
    exact_tests +:= 1;
    try
        G, phi := TorsionSubgroup(J);
        inv := Invariants(G);
    catch e
        print "TORSION_ERROR", tup, e`Object;
        return;
    end try;

    has3 := #inv gt 0 and &or [n mod 3 eq 0 : n in inv];
    if has3 then
        simple, pcert, Lp := IrreducibleFrobeniusCertificate(f);
        Append(~hits, <tup, inv, simple, pcert, Lp>);
        print "HIT", tup, "torsion", inv,
              "simple", simple, "prime", pcert, "L", Lp;
    elif printed lt max_print then
        print "exact_no3", tup, "torsion", inv;
        printed +:= 1;
    end if;
end procedure;

procedure ProcessTuple(tup, contact_data, point_count_data,
                       ~boundary13, ~contact_survivors, ~point_survivors,
                       ~no3_bound, ~survived_bound,
                       ~exact_tests, ~hits, ~printed)
    if boundary13_only and not IsBadAtPrime(tup, 13) then
        return;
    end if;
    if boundary13_live_only and not HasSquareCollisionAtPrime(tup, 13) then
        return;
    end if;
    boundary13 +:= 1;

    pass, killp := PassesContactFilters(tup, contact_data);
    if not pass then
        return;
    end if;
    contact_survivors +:= 1;

    point_pass, point_killp := PassesPointCountFilters(tup, point_count_data);
    if not point_pass then
        return;
    end if;
    point_survivors +:= 1;

    f := CurvePolynomialFromTuple(tup);
    if Discriminant(f) eq 0 then
        return;
    end if;

    gcd_bound, used := TorsionGcdBound(f, bound_primes);
    if gcd_bound mod 3 ne 0 then
        no3_bound +:= 1;
        return;
    end if;

    survived_bound +:= 1;
    if printed lt max_print then
        print "BOUND_SURVIVOR", tup, "gcd_bound", gcd_bound, "used", used;
        printed +:= 1;
    end if;

    if exact_tests lt max_exact then
        ExactCheckTuple(tup, ~exact_tests, ~hits, ~printed);
    end if;
end procedure;

procedure RunFinite()
    print "M(2,2,2,4)+3 finite good-reduction diagnostic";
    for p in finite_primes do
        keys, witnesses, raw_count := ContactCoefficientKeys(p);
        F := GF(p);
        good := 0;
        contact := 0;
        curve_keys := {};
        contact_curve_keys := {};
        samples := [];

        for a in F do
            if a eq 0 then continue; end if;
            for b in F do
                if b eq 0 then continue; end if;
                for c in F do
                    if c eq 0 then continue; end if;
                    for d in F do
                        if d eq 0 then continue; end if;
                        squares := [a^2,b^2,c^2,d^2];
                        if #(Set(squares)) ne 4 then
                            continue;
                        end if;

                        key := CurveKeyFromSquares(squares);
                        Include(~curve_keys, key);
                        good +:= 1;
                        if key in keys then
                            contact +:= 1;
                            Include(~contact_curve_keys, key);
                            if #samples lt 5 then
                                Append(~samples, <Z!a,Z!b,Z!c,Z!d,
                                                  key, witnesses[key]>);
                            end if;
                        end if;
                    end for;
                end for;
            end for;
        end for;

        print "p", p,
              "contact_keys", #keys,
              "raw_contact_param", raw_count,
              "good_points", good,
              "contact_points", contact,
              "curvekeys", #curve_keys,
              "contact_curvekeys", #contact_curve_keys;
        print " samples", samples;
    end for;
end procedure;

procedure RunSearch()
    print "M(2,2,2,4)+3 integer tuple search";
    print "height", height,
          "prime_bound", prime_bound,
          "filter_primes", filter_primes,
          "point_prime_bound", point_prime_bound,
          "point_primes", point_primes,
          "boundary13_only", boundary13_only,
          "boundary13_live_only", boundary13_live_only,
          "max_exact", max_exact;
    contact_data := PrecomputeContactData(filter_primes);
    point_count_data := PrecomputePointCountData(point_primes);

    checked := 0;
    boundary13 := 0;
    contact_survivors := 0;
    point_survivors := 0;
    no3_bound := 0;
    survived_bound := 0;
    exact_tests := 0;
    hits := [];
    printed := 0;

    for a in [1..height-3] do
        for b in [a+1..height-2] do
            for c in [b+1..height-1] do
                for d in [c+1..height] do
                    tup := [a,b,c,d];
                    if GCD(tup) ne 1 then
                        continue;
                    end if;
                    checked +:= 1;

                    if progress_interval gt 0 and checked mod progress_interval eq 0 then
                        print "progress", checked,
                              "boundary13", boundary13,
                              "contact_survivors", contact_survivors,
                              "point_survivors", point_survivors,
                              "survived_bound", survived_bound,
                              "exact", exact_tests,
                              "hits", #hits;
                    end if;

                    ProcessTuple(tup, contact_data, point_count_data,
                                 ~boundary13, ~contact_survivors, ~point_survivors,
                                 ~no3_bound, ~survived_bound,
                                 ~exact_tests, ~hits, ~printed);
                end for;
            end for;
        end for;
    end for;

    print "SUMMARY checked", checked,
          "boundary13", boundary13,
          "contact_survivors", contact_survivors,
          "point_survivors", point_survivors,
          "no3_bound", no3_bound,
          "survived_bound", survived_bound,
          "exact", exact_tests,
          "hits", #hits;
    for H in hits do
        print H;
    end for;
end procedure;

procedure RunExactFile()
    tuples := ReadTupleFile(exact_tuple_file);
    print "M(2,2,2,4)+3 tuple-file search";
    print "tuple_file", exact_tuple_file,
          "tuples", #tuples,
          "prime_bound", prime_bound,
          "filter_primes", filter_primes,
          "point_prime_bound", point_prime_bound,
          "point_primes", point_primes,
          "boundary13_only", boundary13_only,
          "boundary13_live_only", boundary13_live_only,
          "max_exact", max_exact;
    contact_data := PrecomputeContactData(filter_primes);
    point_count_data := PrecomputePointCountData(point_primes);

    boundary13 := 0;
    contact_survivors := 0;
    point_survivors := 0;
    no3_bound := 0;
    survived_bound := 0;
    exact_tests := 0;
    hits := [];
    printed := 0;

    for tup in tuples do
        ProcessTuple(tup, contact_data, point_count_data,
                     ~boundary13, ~contact_survivors, ~point_survivors,
                     ~no3_bound, ~survived_bound,
                     ~exact_tests, ~hits, ~printed);
    end for;

    print "SUMMARY tuples", #tuples,
          "boundary13", boundary13,
          "contact_survivors", contact_survivors,
          "point_survivors", point_survivors,
          "no3_bound", no3_bound,
          "survived_bound", survived_bound,
          "exact", exact_tests,
          "hits", #hits;
    for H in hits do
        print H;
    end for;
end procedure;

if mode eq "finite" then
    RunFinite();
elif mode eq "search" then
    RunSearch();
elif mode eq "exact_file" then
    RunExactFile();
else
    error "Unknown mode. Use mode:=\"finite\", \"search\", or \"exact_file\".";
end if;

quit;
