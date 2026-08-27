//////////////////////////////////////////////////////////////////////
//  Contact-7 rational-root family plus possible rational 3-torsion.
//
//  In the contact-7 family
//
//      h = 1 - (7/2)*x + a*x^2 + b*x^3,
//      f = (h^2 + (x - 1)^7)/x^2,
//
//  the point P=(1,h(1)) gives a rational 7-torsion class.
//  Force a rational Weierstrass root by writing
//
//      r = 1 - s^2,        h(r) = eps*s^7, eps = +/-1,
//
//  which determines
//
//      a = (eps*s^7 - 1 + (7/2)*r - b*r^3)/r^2.
//
//  Then the curve has rational 2-torsion as well.  This script searches
//  for rational 3-torsion by applying the necessary condition that the
//  prime-to-p part of 42 divides #J(F_p) at every good prime p.
//
//  Modes:
//      finite    enumerate the (s,b,eps) chart over small finite fields
//      search    rational height search
//
//  Typical runs:
//      magma -b mode:="finite" code/contact7_root_plus3_search.m
//      magma -b mode:="search" height:=12 code/contact7_root_plus3_search.m
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
    max_exact := 1000;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned progress_interval then
    progress_interval := 20000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;
if not assigned verify_all then
    verify_all := false;
elif Type(verify_all) eq MonStgElt then
    verify_all := verify_all in {"true", "True", "1", "yes"};
end if;

Q := Rationals();
P<x> := PolynomialRing(Q);

target := 42;
finite_primes := [3,5,7,11,13,17,19,23,29,31,37,41,43];
filter_primes := [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97];

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

function RequiredPrimeToP(n, p)
    required := n;
    while required mod p eq 0 do
        required div:= p;
    end while;
    return required;
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

function Contact7RootPolynomial(s, b, eps)
    r := 1 - s^2;
    if r eq 0 then
        return false, P!0, P!0, Q!0, Q!0;
    end if;
    a := (eps*s^7 - 1 + (Q!7/2)*r - b*r^3)/r^2;
    h := 1 - (Q!7/2)*x + a*x^2 + b*x^3;
    num := h^2 + (x - 1)^7;
    if Coefficient(num, 0) ne 0 or Coefficient(num, 1) ne 0 then
        return false, P!0, h, a, r;
    end if;
    f := ExactQuotient(num, x^2);
    if Evaluate(f, r) ne 0 then
        return false, f, h, a, r;
    end if;
    return true, f, h, a, r;
end function;

function Contact7RootPolynomialFinite(F, s, b, eps)
    PF<X> := PolynomialRing(F);
    r := 1 - s^2;
    if r eq 0 then
        return false, PF!0, PF!0, F!0, F!0;
    end if;
    a := (eps*s^7 - 1 + (F!7/F!2)*r - b*r^3)/r^2;
    h := 1 - (F!7/F!2)*X + a*X^2 + b*X^3;
    f := ExactQuotient(h^2 + (X - 1)^7, X^2);
    if Evaluate(f, r) ne 0 then
        return false, f, h, a, r;
    end if;
    return true, f, h, a, r;
end function;

function PassesTargetReduction(f)
    used := [];
    for p in filter_primes do
        try
            fp := ChangeRing(f, GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
                continue;
            end if;
            C := HyperellipticCurve(fp);
            n := Integers()!#Jacobian(C);
            Append(~used, <p,n>);
            required := RequiredPrimeToP(target, p);
            if required gt 1 and n mod required ne 0 then
                return false, p, n, used;
            end if;
        catch e
            continue;
        end try;
    end for;
    return true, 0, 0, used;
end function;

procedure RunFinite()
    print "FINITE contact-7 rational-root plus 3 necessary condition";
    print "target", target, "finite_primes", finite_primes;

    for p in finite_primes do
        F := GF(p);
        total := 0;
        good := 0;
        pass3 := 0;
        pass_target := 0;
        samples := [];

        for s in F do
            for b in F do
                for eps in [F!-1, F!1] do
                    total +:= 1;
                    ok, f, h, a, r := Contact7RootPolynomialFinite(F, s, b, eps);
                    if not ok or Degree(f) ne 5 or Discriminant(f) eq 0 then
                        continue;
                    end if;
                    if Evaluate(h, F!1) eq 0 then
                        continue;
                    end if;
                    good +:= 1;
                    C := HyperellipticCurve(f);
                    n := Integers()!#Jacobian(C);
                    if p eq 3 or n mod 3 eq 0 then
                        pass3 +:= 1;
                    end if;
                    required := RequiredPrimeToP(target, p);
                    if required eq 1 or n mod required eq 0 then
                        pass_target +:= 1;
                        if #samples lt 10 then
                            Append(~samples, <Integers()!s, Integers()!b, Integers()!eps, Integers()!a, Integers()!r, n>);
                        end if;
                    end if;
                end for;
            end for;
        end for;

        print "p", p, "total", total, "good", good,
              "pass3", pass3, "pass_target", pass_target,
              "samples", samples;
    end for;
end procedure;

procedure RunSearch()
    params := RationalParametersOfHeight(height);
    checked := 0;
    smooth := 0;
    verified7 := 0;
    target_survivors := 0;
    exact_tests := 0;
    hits := [];
    first_kill := AssociativeArray();
    torsion_counts := AssociativeArray();

    print "SEARCH contact-7 rational-root plus 3";
    print "height", height, "parameters", #params, "target", target,
          "max_exact", max_exact, "verify_all", verify_all;
    print "filter_primes", filter_primes;

    for s in params do
        for b in params do
            for eps in [Q!-1, Q!1] do
                checked +:= 1;
                if progress_interval gt 0 and checked mod progress_interval eq 0 then
                    print "progress", "checked", checked, "smooth", smooth,
                          "verified7", verified7, "target_survivors", target_survivors,
                          "exact_tests", exact_tests, "hits", #hits;
                end if;

                ok, f, h, a, r := Contact7RootPolynomial(s, b, eps);
                if not ok or not GoodHyperellipticPolynomial(f) then
                    continue;
                end if;
                if Evaluate(h, Q!1) eq 0 then
                    continue;
                end if;
                smooth +:= 1;

                if verify_all then
                    C0 := HyperellipticCurve(f);
                    J0 := Jacobian(C0);
                    if Order(J0![x - 1, Evaluate(h, Q!1)]) ne 7 then
                        continue;
                    end if;
                end if;
                verified7 +:= 1;

                pass, pbad, nbad, used := PassesTargetReduction(f);
                if not pass then
                    if IsDefined(first_kill, pbad) then
                        first_kill[pbad] +:= 1;
                    else
                        first_kill[pbad] := 1;
                    end if;
                    continue;
                end if;

                target_survivors +:= 1;
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

                if HasInvariantDivisibleBy(invs, target) then
                    simple, pcert, Lp := SimpleCertificate(fI);
                    Append(~hits, <s,b,eps,a,r,invs,simple,pcert,fI,used>);
                    print "HIT42", "s", s, "b", b, "eps", eps,
                          "a", a, "r", r, "torsion", invs,
                          "order", TorsionOrder(invs),
                          "simple", simple, "pcert", pcert;
                    print "  f =", fI;
                    print "  factorization =", Factorization(f);
                    print "  used", used;
                else
                    print "SURVIVOR_NOT_42", "s", s, "b", b, "eps", eps,
                          "a", a, "r", r, "torsion", invs, "used", used;
                end if;
            end for;
        end for;
    end for;

    print "DONE contact-7 rational-root plus 3";
    print "checked", checked, "smooth", smooth, "verified7", verified7,
          "target_survivors", target_survivors, "exact_tests", exact_tests,
          "hits", #hits;
    print "FIRST_KILL";
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
