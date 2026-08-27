//////////////////////////////////////////////////////////////////////
//  M_1(8,4) tangent-cover family plus possible rational 3-torsion.
//
//  The M_1(8,4) tangent-cover family already gives many simple
//  Jacobians with torsion [4,8].  This script asks whether the same
//  chart can also have rational 3-torsion.  A necessary condition is
//  that 3 divides #J(F_p) for every good reduction prime p != 3.
//
//  Modes:
//      finite      finite-field density check on the tangent-cover chart
//      search      exact rational search with reduction filter
//      prime_diag  prime-by-prime kill diagnostic for rational candidates
//
//  Typical runs:
//      magma -b mode:="finite" code/elkies3_m14_plus3_search.m
//      magma -b mode:="search" height:=50 code/elkies3_m14_plus3_search.m
//      magma -b mode:="prime_diag" height:=50 code/elkies3_m14_plus3_search.m
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
Z := Integers();
Qx<x> := PolynomialRing(Q);

filter_primes := [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61];
finite_primes := [5,7,11,13,17,19,23,29,31,37,41,43];

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

function IsSquareQ(q)
    q := Q!q;
    if q lt 0 then
        return false;
    end if;
    return IsSquare(Numerator(q)) and IsSquare(Denominator(q));
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return Qx!(L^2*f), L;
end function;

function GoodHyperellipticPolynomial(f)
    return Degree(f) in {5,6} and Discriminant(f) ne 0;
end function;

function HasThreeInvariants(invs)
    return &or [ n mod 3 eq 0 : n in invs ];
end function;

function M14FamilyPolynomial(R, w)
    m := R;
    n := Q!1;
    t := (2*m^2 + (1-w^2)*m*n - 2*w^2*n^2)/(4*(w^2-1));
    A := n^4*x^2
         + (m^3*n + 4*m^2*t + m*n^3 - 8*m*n*t + 4*n^2*t)*x
         + m^4;
    B := (m*n + 2*n^2 + 4*t)*x^2
         + (m^2 + 4*m*n + n^2 + 8*t)*x
         + (2*m^2 + m*n + 4*t);
    return x*A*B, t, A, B;
end function;

function M14FamilyPolynomialFinite(K, R, w)
    P<X> := PolynomialRing(K);
    m := R;
    n := K!1;
    t := (2*m^2 + (1-w^2)*m*n - 2*w^2*n^2)/(4*(w^2-1));
    A := n^4*X^2
         + (m^3*n + 4*m^2*t + m*n^3 - 8*m*n*t + 4*n^2*t)*X
         + m^4;
    B := (m*n + 2*n^2 + 4*t)*X^2
         + (m^2 + 4*m*n + n^2 + 8*t)*X
         + (2*m^2 + m*n + 4*t);
    return X*A*B, t;
end function;

function PlusDisc(R, w)
    return -4*(w-R)*(R-1)^2*(R+1)*(R+w)*(w+1)
            *(R*w - 3*R + 3*w - 1);
end function;

function MinusDisc(R, w)
    return 4*(w-1)*(R+1)*(R*w + 3*R + 3*w + 1)
           *(R^4 - 2*R^3 + R^2*w^2 - R^2 + 2*R*w^2 - w^2);
end function;

function M14CoverPossible(R, w)
    return IsSquareQ(PlusDisc(R,w)) or IsSquareQ(MinusDisc(R,w));
end function;

function PassesThreeReduction(f, primes)
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
        if (#Jacobian(C) mod 3) ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

procedure FiniteM14Plus3()
    print "FINITE M_1(8,4) tangent-cover chart: need 3 | #J(Fp)";
    print "finite_primes", finite_primes;
    for p in finite_primes do
        K := GF(p);
        checked := 0;
        cover := 0;
        good := 0;
        three := 0;
        samples := [];
        for R in K do
            for w in K do
                if R eq 0 or w eq 0 or w eq 1 or w eq -1 then
                    continue;
                end if;
                checked +:= 1;
                if not (IsSquare(K!PlusDisc(R,w)) or IsSquare(K!MinusDisc(R,w))) then
                    continue;
                end if;
                cover +:= 1;
                f, t := M14FamilyPolynomialFinite(K, R, w);
                if not GoodHyperellipticPolynomial(f) then
                    continue;
                end if;
                good +:= 1;
                C := HyperellipticCurve(f);
                if (#Jacobian(C) mod 3) eq 0 then
                    three +:= 1;
                    if #samples lt 8 then
                        Append(~samples, <Z!R,Z!w>);
                    end if;
                end if;
            end for;
        end for;
        print "p", p, "checked", checked, "cover", cover, "good", good,
              "three_possible", three, "samples", samples;
    end for;
end procedure;

procedure SearchM14Plus3()
    print "SEARCH M_1(8,4) tangent cover plus 3";
    print "height", height, "filter_primes", filter_primes;
    params := RationalParametersOfHeight(height);
    checked := 0;
    cover := 0;
    smooth := 0;
    three_survivors := 0;
    verified := 0;
    torsion_tests := 0;
    hits := [];
    stop := false;

    for R in params do
        if stop then break; end if;
        for w in params do
            if R eq 0 or w in {Q!-1, Q!0, Q!1} then
                continue;
            end if;
            checked +:= 1;
            if not M14CoverPossible(R,w) then
                continue;
            end if;
            cover +:= 1;
            f, t, A, B := M14FamilyPolynomial(R,w);
            if not GoodHyperellipticPolynomial(f) then
                continue;
            end if;
            smooth +:= 1;
            if not PassesThreeReduction(f, filter_primes) then
                continue;
            end if;
            three_survivors +:= 1;

            fI, L := IntegralModelPolynomial(f);
            C := HyperellipticCurve(fI);
            J := Jacobian(C);
            D := J![x, Q!0];
            ok, half := IsDivisibleBy(D, 2);
            if not ok then
                continue;
            end if;
            verified +:= 1;
            torsion_tests +:= 1;
            G, phi := TorsionSubgroup(J);
            invs := Invariants(G);
            print "SURVIVOR", "R", R, "w", w, "t", t,
                  "half_order", Order(half), "torsion", invs;
            if HasThreeInvariants(invs) then
                Append(~hits, <R,w,t,invs,fI>);
                print "HIT", "R", R, "w", w, "torsion", invs, "f", fI;
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
    end for;

    print "DONE plus3 height", height, "checked", checked, "cover", cover,
          "smooth", smooth, "three_survivors", three_survivors,
          "verified", verified, "torsion_tests", torsion_tests, "hits", #hits;
end procedure;

procedure PrimeDiagnostic()
    print "PRIME_DIAG M_1(8,4) tangent cover plus 3";
    print "height", height, "filter_primes", filter_primes;
    params := RationalParametersOfHeight(height);

    checked := 0;
    cover := 0;
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

    for R in params do
        for w in params do
            if R eq 0 or w in {Q!-1, Q!0, Q!1} then
                continue;
            end if;
            checked +:= 1;
            if not M14CoverPossible(R,w) then
                continue;
            end if;
            cover +:= 1;
            f, t, A, B := M14FamilyPolynomial(R,w);
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
                if (#Jacobian(C) mod 3) eq 0 then
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
    print "cover", cover;
    print "smooth", smooth;
    print "pass_all", pass_all;
    for p in filter_primes do
        print "PRIME", p,
              "good", good_counts[p],
              "bad_or_boundary", bad_counts[p],
              "pass3", pass_counts[p],
              "kill3", kill_counts[p],
              "first_kill", first_kill_counts[p];
    end for;
end procedure;

if mode eq "finite" then
    FiniteM14Plus3();
elif mode eq "search" then
    SearchM14Plus3();
elif mode eq "prime_diag" then
    PrimeDiagnostic();
else
    error "Unknown mode";
end if;
