//////////////////////////////////////////////////////////////////////
//  M(12) plus possible rational 7-torsion.
//
//  A necessary condition for rational 7-torsion is that 7 divides
//  #J(F_p) for every good reduction prime p != 7.  This script applies
//  that filter to the M(12) family used in m12_simple_search.m.
//
//  Modes:
//      finite      scan the finite affine (z,r) chart
//      search      exact rational search after the reduction filter
//      prime_diag  prime-by-prime kill diagnostic for rational candidates
//////////////////////////////////////////////////////////////////////

if not assigned mode then
    mode := "finite";
end if;
if not assigned height then
    height := 30;
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
Qx<x> := PolynomialRing(Q);
QX<X> := PolynomialRing(Q);

finite_primes := [3,5,11,13,17,19,23,29,31,37,41,43];
filter_primes := [3,5,11,13,17,19,23,29,31,37,41,43,47,53,59,61];

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            q := Q!num/Q!den;
            key := Sprint(q);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, q);
            end if;
        end for;
    end for;
    return vals;
end function;

function M12Polynomial(a, r)
    T := a*x^2 - x + r;
    h := (x-r)*(T+1);
    f := a*x^2*T*(T+1);
    W := h^2 + 4*f;
    return W, T, h;
end function;

function M12PolynomialFinite(K, z, r)
    P<XK> := PolynomialRing(K);
    a := (1-z^2)/(4*(r+1));
    T := a*XK^2 - XK + r;
    h := (XK-r)*(T+1);
    f := a*XK^2*T*(T+1);
    W := h^2 + 4*f;
    return W, T, h, a;
end function;

function OddQuinticAtRoot(W, w)
    out := QX!0;
    for i in [0..Degree(W)] do
        for j in [0..i] do
            out +:= Coefficient(W, i)*Binomial(i,j)*w^(i-j)*X^(6-j);
        end for;
    end for;
    return out;
end function;

function GoodHyperellipticPolynomial(f)
    return Degree(f) in {5,6} and Discriminant(f) ne 0;
end function;

function HasSevenInvariants(invs)
    return &or [ n mod 7 eq 0 : n in invs ];
end function;

function PassesSevenReduction(W, primes)
    for p in primes do
        try
            Wp := ChangeRing(W, GF(p));
        catch e
            continue;
        end try;
        if not GoodHyperellipticPolynomial(Wp) then
            continue;
        end if;
        C := HyperellipticCurve(Wp);
        if (#Jacobian(C) mod 7) ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

procedure FiniteM12Plus7()
    print "FINITE M(12) plus 7 necessary condition";
    print "finite_primes", finite_primes;
    for p in finite_primes do
        K := GF(p);
        total := 0;
        good := 0;
        seven := 0;
        samples := [];
        for z in K do
            for r in K do
                if r eq -1 or z^2 eq 1 then
                    continue;
                end if;
                total +:= 1;
                W, T, h, a := M12PolynomialFinite(K, z, r);
                if a eq 0 or not GoodHyperellipticPolynomial(W) then
                    continue;
                end if;
                good +:= 1;
                C := HyperellipticCurve(W);
                if (#Jacobian(C) mod 7) eq 0 then
                    seven +:= 1;
                    if #samples lt 8 then
                        Append(~samples, <Integers()!z, Integers()!r>);
                    end if;
                end if;
            end for;
        end for;
        print "p", p, "total", total, "good", good,
              "seven_possible", seven, "samples", samples;
    end for;
end procedure;

procedure SearchM12Plus7()
    print "SEARCH M(12) plus 7";
    print "height", height, "filter_primes", filter_primes;
    params := RationalParametersOfHeight(height);
    checked := 0;
    smooth := 0;
    split_T := 0;
    seven_survivors := 0;
    torsion_tests := 0;
    hits := [];
    stop := false;

    for z in params do
        if stop then break; end if;
        for r in params do
            if r eq -1 or z^2 eq 1 then
                continue;
            end if;
            a := (1-z^2)/(4*(r+1));
            if a eq 0 then
                continue;
            end if;
            checked +:= 1;
            W, T, h := M12Polynomial(a, r);
            if not GoodHyperellipticPolynomial(W) then
                continue;
            end if;
            smooth +:= 1;
            roots_quad := Roots(T+1);
            if #roots_quad eq 0 then
                continue;
            end if;
            split_T +:= 1;
            if not PassesSevenReduction(W, filter_primes) then
                continue;
            end if;
            seven_survivors +:= 1;

            for root_data in roots_quad do
                if root_data[2] ne 1 then
                    continue;
                end if;
                w := root_data[1];
                if w eq 0 then
                    continue;
                end if;
                f5 := OddQuinticAtRoot(W, w);
                if Degree(f5) ne 5 or Discriminant(f5) eq 0 then
                    continue;
                end if;
                Y0 := Evaluate(h, 0);
                Xp := -1/w;
                Yp := Y0*Xp^3;
                if Evaluate(f5, Xp) ne Yp^2 then
                    continue;
                end if;
                C5 := HyperellipticCurve(f5);
                J := Jacobian(C5);
                D := J![X-Xp, Yp];
                if Order(D) ne 12 then
                    continue;
                end if;
                torsion_tests +:= 1;
                G, phi := TorsionSubgroup(J);
                invs := Invariants(G);
                print "SURVIVOR", "a", a, "r", r, "z", z, "w", w,
                      "torsion", invs, "f5", f5;
                if HasSevenInvariants(invs) then
                    Append(~hits, <a,r,z,w,invs,f5>);
                    print "HIT", "a", a, "r", r, "z", z, "w", w,
                          "torsion", invs, "f5", f5;
                    if #hits ge max_hits then
                        stop := true;
                        break;
                    end if;
                end if;
                if torsion_tests ge max_tests then
                    stop := true;
                    break;
                end if;
            end for;
            if stop then break; end if;
        end for;
    end for;

    print "DONE m12_plus7 height", height,
          "checked", checked,
          "smooth", smooth,
          "split_T", split_T,
          "seven_survivors", seven_survivors,
          "torsion_tests", torsion_tests,
          "hits", #hits;
end procedure;

procedure PrimeDiagnostic()
    print "PRIME_DIAG M(12) plus 7";
    print "height", height, "filter_primes", filter_primes;
    params := RationalParametersOfHeight(height);
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

    for z in params do
        for r in params do
            if r eq -1 or z^2 eq 1 then
                continue;
            end if;
            a := (1-z^2)/(4*(r+1));
            if a eq 0 then
                continue;
            end if;
            checked +:= 1;
            W, T, h := M12Polynomial(a, r);
            if not GoodHyperellipticPolynomial(W) then
                continue;
            end if;
            smooth +:= 1;
            first_fail := 0;
            for p in filter_primes do
                try
                    Wp := ChangeRing(W, GF(p));
                catch e
                    bad_counts[p] +:= 1;
                    continue;
                end try;
                if not GoodHyperellipticPolynomial(Wp) then
                    bad_counts[p] +:= 1;
                    continue;
                end if;
                good_counts[p] +:= 1;
                C := HyperellipticCurve(Wp);
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
    FiniteM12Plus7();
elif mode eq "search" then
    SearchM12Plus7();
elif mode eq "prime_diag" then
    PrimeDiagnostic();
else
    error "Unknown mode";
end if;
