//////////////////////////////////////////////////////////////////////
//  Simple 35-torsion attempt: contact-5 family plus rational 7-torsion.
//
//  The contact-5 family is
//
//      h = 1 + a*x + b*x^2,
//      f = h^2 - (1+a+b)^2*x^5.
//
//  For smooth degree-5 specializations with 1+a+b != 0, the class
//      J![x, 1]
//  has order 5.  This script searches for an additional rational
//  7-torsion class using the necessary condition that, at every good
//  reduction prime p != 7, 7 divides #J(F_p).  Exact torsion is computed
//  only for survivors.
//
//  Modes:
//      finite    enumerate (a,b) over F_p and count 7-divisible #J(F_p)
//      search    rational parameter height search
//
//  Typical runs:
//      magma -b mode:="finite" code/contact5_plus7_search.m
//      magma -b mode:="search" height:=12 code/contact5_plus7_search.m
//////////////////////////////////////////////////////////////////////

if not assigned mode then
    mode := "finite";
end if;
if not assigned height then
    height := 12;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned max_exact then
    max_exact := 500;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned progress_interval then
    progress_interval := 20000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Q := Rationals();
P<x> := PolynomialRing(Q);

finite_primes := [3,5,11,13,17,19,23,29,31,37,41,43];
filter_primes := [3,5,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97];

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

function IntegralModel(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function Contact5Polynomial(a, b)
    h := 1 + a*x + b*x^2;
    k := (1+a+b)^2;
    return h^2 - k*x^5, h, k;
end function;

function GoodHyperellipticPolynomial(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function TorsionOrder(invs)
    if #invs eq 0 then
        return 1;
    end if;
    return &*invs;
end function;

function HasInvariantDivisibleBy(invs, n)
    return &or [ m mod n eq 0 : m in invs ];
end function;

function SimpleCertificate(f)
    C := HyperellipticCurve(f);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149] do
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
    return false, 0, P!0;
end function;

function PassesSevenReduction(f)
    used := [];
    for p in filter_primes do
        if p eq 7 then
            continue;
        end if;
        try
            fp := ChangeRing(f, GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
                continue;
            end if;
            C := HyperellipticCurve(fp);
            n := Integers()!#Jacobian(C);
            Append(~used, <p,n>);
            if n mod 7 ne 0 then
                return false, p, n, used;
            end if;
        catch e
            continue;
        end try;
    end for;
    return true, 0, 0, used;
end function;

procedure RunFinite()
    print "FINITE contact-5 plus 7 necessary condition";
    print "finite_primes", finite_primes;
    for p in finite_primes do
        F := GF(p);
        PF<X> := PolynomialRing(F);
        total := 0;
        good := 0;
        pass7 := 0;
        samples := [];

        for a in F do
            for b in F do
                total +:= 1;
                h := 1 + a*X + b*X^2;
                k := (1+a+b)^2;
                if k eq 0 then
                    continue;
                end if;
                f := h^2 - k*X^5;
                if Degree(f) ne 5 or Discriminant(f) eq 0 then
                    continue;
                end if;
                good +:= 1;
                C := HyperellipticCurve(f);
                n := Integers()!#Jacobian(C);
                if n mod 7 eq 0 then
                    pass7 +:= 1;
                    if #samples lt 10 then
                        Append(~samples, <Integers()!a, Integers()!b, n>);
                    end if;
                end if;
            end for;
        end for;

        print "p", p, "total", total, "good", good,
              "pass7", pass7, "samples", samples;
    end for;
end procedure;

procedure RunSearch()
    params := RationalParametersOfHeight(height);
    checked := 0;
    smooth := 0;
    seven_survivors := 0;
    exact_tests := 0;
    hits := [];
    first_kill := AssociativeArray();
    torsion_counts := AssociativeArray();

    print "SEARCH contact-5 plus 7";
    print "height", height, "parameters", #params,
          "max_exact", max_exact, "filter_primes", filter_primes;

    for a in params do
        for b in params do
            checked +:= 1;
            if progress_interval gt 0 and checked mod progress_interval eq 0 then
                print "progress", "checked", checked, "smooth", smooth,
                      "seven_survivors", seven_survivors,
                      "exact_tests", exact_tests, "hits", #hits;
            end if;

            f, h, k := Contact5Polynomial(a, b);
            if k eq 0 or not GoodHyperellipticPolynomial(f) then
                continue;
            end if;
            smooth +:= 1;

            pass, pbad, nbad, used := PassesSevenReduction(f);
            if not pass then
                if IsDefined(first_kill, pbad) then
                    first_kill[pbad] +:= 1;
                else
                    first_kill[pbad] := 1;
                end if;
                continue;
            end if;
            seven_survivors +:= 1;
            if exact_tests ge max_exact then
                continue;
            end if;

            fI, L := IntegralModel(f);
            C := HyperellipticCurve(fI);
            J := Jacobian(C);
            G, phi := TorsionSubgroup(J);
            invs := Invariants(G);
            exact_tests +:= 1;

            invkey := Sprint(invs);
            if IsDefined(torsion_counts, invkey) then
                torsion_counts[invkey] +:= 1;
            else
                torsion_counts[invkey] := 1;
            end if;

            if HasInvariantDivisibleBy(invs, 35) then
                simple, pcert, Lp := SimpleCertificate(fI);
                Append(~hits, <a,b,invs,simple,pcert,fI,used>);
                print "HIT35", "a", a, "b", b, "torsion", invs,
                      "order", TorsionOrder(invs),
                      "simple", simple, "pcert", pcert;
                print "  f =", fI;
                print "  used", used;
            else
                print "SURVIVOR_NOT_35", "a", a, "b", b,
                      "torsion", invs, "used", used;
            end if;
        end for;
    end for;

    print "DONE contact-5 plus 7";
    print "checked", checked, "smooth", smooth,
          "seven_survivors", seven_survivors,
          "exact_tests", exact_tests, "hits", #hits;
    print "FIRST_KILL_7";
    for p in Sort([ k : k in Keys(first_kill) ]) do
        print " ", p, first_kill[p];
    end for;
    print "Torsion counts among exact survivors";
    for key in Sort([ k : k in Keys(torsion_counts) ]) do
        print " ", key, torsion_counts[key];
    end for;
end procedure;

if mode eq "finite" then
    RunFinite();
elif mode eq "search" then
    RunSearch();
else
    print "unknown mode", mode;
end if;

quit;
