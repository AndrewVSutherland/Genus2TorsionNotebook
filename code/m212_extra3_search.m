//////////////////////////////////////////////////////////////////////
//  M(2,12) plus an independent rational 3-torsion point.
//
//  The M(2,12) chart already has a marked order-12 class D, hence the
//  built-in rational 3-torsion class 4D.  An additional independent
//  rational 3-torsion class forces, at every good prime p != 3, the
//  finite group J(F_p)[3] to have rank at least 2.
//
//  Modes:
//    finite: scan the finite (z,r) chart and count good residues with
//            3-rank at least 2 in J(F_p).
//    search: rational-height search using the same 3-rank filter before
//            exact torsion.
//
//  Typical runs:
//      magma -b mode:=finite code/m212_extra3_search.m
//      magma -b mode:=search height:=20 prime_bound:=43 code/m212_extra3_search.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned mode then
    mode := "finite";
end if;
if not assigned height then
    height := 20;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned prime_bound then
    prime_bound := 43;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;
if not assigned max_exact then
    max_exact := 100;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned max_print then
    max_print := 20;
elif Type(max_print) eq MonStgElt then
    max_print := StringToInteger(max_print);
end if;
if not assigned progress_interval then
    progress_interval := 10000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Q := Rationals();
Z := Integers();
Qx<x> := PolynomialRing(Q);
QX<X> := PolynomialRing(Q);

finite_primes := [5,7,11,13,17,19,23,29,31,37,41,43];

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            q := Q!num / Q!den;
            key := Sprint(q);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, q);
            end if;
        end for;
    end for;
    return vals;
end function;

function M212Polynomial(a, r)
    T := a*x^2 - x + r;
    h := (x-r)*(T+1);
    W := h^2 + 4*a*x^2*T*(T+1);
    return W, T, h;
end function;

function M212PolynomialFromZR(z, r)
    if r eq -1 or z^2 eq 1 then
        return false, Qx!0, Qx!0, Qx!0, Q!0;
    end if;
    a := (1-z^2)/(4*(r+1));
    if a eq 0 then
        return false, Qx!0, Qx!0, Qx!0, Q!0;
    end if;
    W, T, h := M212Polynomial(a, r);
    return true, W, T, h, a;
end function;

function M212PolynomialFinite(K, z, r)
    P<t> := PolynomialRing(K);
    if r eq -K!1 or z^2 eq K!1 then
        return false, P!0, P!0, K!0;
    end if;
    a := (1-z^2)/(4*(r+1));
    if a eq 0 then
        return false, P!0, P!0, K!0;
    end if;
    T := a*t^2 - t + r;
    h := (t-r)*(T+1);
    W := h^2 + 4*a*t^2*T*(T+1);
    return true, W, T, a;
end function;

function GoodHyperellipticPolynomial(f)
    return Degree(f) in {5,6} and Discriminant(f) ne 0;
end function;

function ThreeRankFromInvariants(invs)
    return #[n : n in invs | (Z!n) mod 3 eq 0];
end function;

function JacobianInvariantsFinite(C)
    J := Jacobian(C);
    A, phi := AbelianGroup(J);
    return Invariants(A);
end function;

function ThreeRankFinitePolynomial(f)
    C := HyperellipticCurve(f);
    invs := JacobianInvariantsFinite(C);
    return ThreeRankFromInvariants(invs), invs;
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return Qx!(L^2*f), L;
end function;

function TorsionOrder(invs)
    if #invs eq 0 then
        return 1;
    end if;
    return &*invs;
end function;

function HasIndependentThreeInvariants(invs)
    return ThreeRankFromInvariants(invs) ge 2;
end function;

procedure FiniteM212Extra3()
    print "FINITE M(2,12) plus independent 3-torsion necessary condition";
    print "Need 3-rank(J(F_p)[3]) >= 2 at every good p != 3.";
    for p in finite_primes do
        K := GF(p);
        total := 0;
        good := 0;
        rank0 := 0;
        rank1 := 0;
        rank2 := 0;
        rank3plus := 0;
        bad := 0;
        samples := [];
        for zz in K do
            for rr in K do
                total +:= 1;
                ok, W, T, a := M212PolynomialFinite(K, zz, rr);
                if not ok or not GoodHyperellipticPolynomial(W) then
                    bad +:= 1;
                    continue;
                end if;
                good +:= 1;
                trank, invs := ThreeRankFinitePolynomial(W);
                if trank eq 0 then
                    rank0 +:= 1;
                elif trank eq 1 then
                    rank1 +:= 1;
                elif trank eq 2 then
                    rank2 +:= 1;
                    if #samples lt 8 then
                        Append(~samples, <Z!zz, Z!rr, invs>);
                    end if;
                else
                    rank3plus +:= 1;
                    if #samples lt 8 then
                        Append(~samples, <Z!zz, Z!rr, invs>);
                    end if;
                end if;
            end for;
        end for;
        print "p", p, "total", total, "good", good, "bad", bad,
              "rank0", rank0, "rank1", rank1,
              "rank2", rank2, "rank3plus", rank3plus,
              "samples", samples;
    end for;
end procedure;

function PassesExtra3Reduction(W, primes)
    bad_primes := [];
    used := [];
    for p in primes do
        try
            Wp := ChangeRing(W, GF(p));
        catch e
            Append(~bad_primes, p);
            continue;
        end try;
        if not GoodHyperellipticPolynomial(Wp) then
            Append(~bad_primes, p);
            continue;
        end if;
        trank, invs := ThreeRankFinitePolynomial(Wp);
        Append(~used, <p, trank, invs>);
        if trank lt 2 then
            return false, p, trank, invs, bad_primes, used;
        end if;
    end for;
    return true, 0, 0, [], bad_primes, used;
end function;

procedure SearchM212Extra3()
    primes := [p : p in PrimesUpTo(prime_bound) | p notin {2,3}];
    params := RationalParametersOfHeight(height);
    checked := 0;
    smooth := 0;
    reduction_survivors := 0;
    exact_tests := 0;
    hits := [];
    printed := 0;
    kill_counts := AssociativeArray();

    print "SEARCH M(2,12) plus independent 3-torsion";
    print "height", height, "parameters", #params,
          "prime_bound", prime_bound, "primes", primes,
          "max_exact", max_exact;

    for z in params do
        for r in params do
            checked +:= 1;
            if progress_interval gt 0 and checked mod progress_interval eq 0 then
                print "progress", checked, "smooth", smooth,
                      "reduction_survivors", reduction_survivors,
                      "exact", exact_tests, "hits", #hits;
            end if;
            ok, W, T, h, a := M212PolynomialFromZR(z, r);
            if not ok or not GoodHyperellipticPolynomial(W) then
                continue;
            end if;
            smooth +:= 1;
            pass, killp, trank, invs_kill, bad_primes, used := PassesExtra3Reduction(W, primes);
            if not pass then
                if IsDefined(kill_counts, killp) then
                    kill_counts[killp] +:= 1;
                else
                    kill_counts[killp] := 1;
                end if;
                continue;
            end if;
            reduction_survivors +:= 1;
            if printed lt max_print then
                print "SURVIVOR", "z", z, "r", r, "a", a,
                      "bad_primes", bad_primes, "used", used;
                printed +:= 1;
            end if;
            if exact_tests ge max_exact then
                continue;
            end if;
            WI, L := IntegralModelPolynomial(W);
            C := HyperellipticCurve(WI);
            J := Jacobian(C);
            G, phi := TorsionSubgroup(J);
            invs := Invariants(G);
            exact_tests +:= 1;
            print "EXACT", "z", z, "r", r, "a", a,
                  "torsion", invs, "order", TorsionOrder(invs),
                  "three_rank", ThreeRankFromInvariants(invs);
            if HasIndependentThreeInvariants(invs) then
                Append(~hits, <z,r,a,invs,WI>);
                print "EXTRA3_HIT", "z", z, "r", r, "a", a,
                      "torsion", invs, "W", WI;
            end if;
        end for;
    end for;

    print "DONE M(2,12)+extra3 search";
    print "checked", checked, "smooth", smooth,
          "reduction_survivors", reduction_survivors,
          "exact_tests", exact_tests, "hits", #hits;
    print "kill_counts";
    for p in Sort([k : k in Keys(kill_counts)]) do
        print p, kill_counts[p];
    end for;
end procedure;

if mode eq "finite" then
    FiniteM212Extra3();
elif mode eq "search" then
    SearchM212Extra3();
else
    error "Unknown mode. Use mode:=finite or mode:=search.";
end if;

quit;
