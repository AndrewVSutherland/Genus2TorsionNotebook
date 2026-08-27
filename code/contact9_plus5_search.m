//////////////////////////////////////////////////////////////////////
//  Contact-9 direct family plus possible rational 5-torsion.
//
//  The contact-9 family
//
//      h = 1 - (9/2)*x + (63/8)*x^2 - (105/16)*x^3 + a*x^4,
//      f = (h^2 + (x - 1)^9)/x^4
//
//  has a rational marked 9-torsion class for smooth curves with h(1) != 0.
//  This script searches for extra 5-torsion by requiring the prime-to-p part
//  of target 45 to divide #J(F_p) at good primes.
//
//  Typical runs:
//      magma -b mode:=finite code/contact9_plus5_search.m
//      magma -b mode:=search height:=200 code/contact9_plus5_search.m
//////////////////////////////////////////////////////////////////////

if not assigned mode then
    mode := "finite";
end if;
if not assigned height then
    height := 100;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned max_exact then
    max_exact := 1000;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned progress_interval then
    progress_interval := 2000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Q := Rationals();
P<x> := PolynomialRing(Q);

target := 45;
finite_primes := [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97];
filter_primes := [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149];

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

function Contact9Polynomial(a)
    h := 1 - (Q!9/2)*x + (Q!63/8)*x^2 - (Q!105/16)*x^3 + a*x^4;
    f := ExactQuotient(h^2 + (x - 1)^9, x^4);
    return f, h;
end function;

function Contact9PolynomialFinite(F, a)
    PF<X> := PolynomialRing(F);
    h := 1 - (F!9/F!2)*X + (F!63/F!8)*X^2 - (F!105/F!16)*X^3 + a*X^4;
    f := ExactQuotient(h^2 + (X - 1)^9, X^4);
    return f, h;
end function;

function GoodContactCurve(f, h)
    return Degree(f) eq 5 and Discriminant(f) ne 0 and Evaluate(h, Q!1) ne 0;
end function;

function SimpleCertificate(f)
    C := HyperellipticCurve(f);
    for p in filter_primes do
        try
            fp := ChangeRing(f, GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
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

function PassesTargetReduction(f)
    used := [];
    for p in filter_primes do
        try
            fp := ChangeRing(f, GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
                continue;
            end if;
            C := HyperellipticCurve(fp);
            jp := Integers()!#Jacobian(C);
            Append(~used, <p,jp>);
            required := RequiredPrimeToP(target, p);
            if required gt 1 and jp mod required ne 0 then
                return false, p, jp, used;
            end if;
        catch e
            continue;
        end try;
    end for;
    return true, 0, 0, used;
end function;

procedure RunFinite()
    print "CONTACT9 plus 5 finite diagnostic";
    print "target", target, "finite_primes", finite_primes;
    for p in finite_primes do
        F := GF(p);
        total := 0;
        good := 0;
        pass_target := 0;
        samples := [];
        for a in F do
            total +:= 1;
            f, h := Contact9PolynomialFinite(F, a);
            if Degree(f) ne 5 or Discriminant(f) eq 0 or Evaluate(h, F!1) eq 0 then
                continue;
            end if;
            good +:= 1;
            C := HyperellipticCurve(f);
            jp := Integers()!#Jacobian(C);
            required := RequiredPrimeToP(target, p);
            if required eq 1 or jp mod required eq 0 then
                pass_target +:= 1;
                if #samples lt 12 then
                    Append(~samples, <Integers()!a, jp>);
                end if;
            end if;
        end for;
        print "p", p, "total", total, "good", good,
              "pass_target", pass_target, "samples", samples;
    end for;
end procedure;

procedure RunSearch()
    vals := RationalParametersOfHeight(height);
    checked := 0;
    smooth := 0;
    marked9 := 0;
    target_survivors := 0;
    exact_tests := 0;
    hits := 0;
    first_kill := AssociativeArray();
    torsion_counts := AssociativeArray();

    print "CONTACT9 plus 5 rational search";
    print "height", height, "params", #vals, "target", target,
          "max_exact", max_exact;
    print "filter_primes", filter_primes;

    for a in vals do
        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "progress", "checked", checked, "smooth", smooth,
                  "marked9", marked9, "target_survivors", target_survivors,
                  "exact_tests", exact_tests, "hits", hits;
        end if;

        f, h := Contact9Polynomial(a);
        if not GoodContactCurve(f, h) then
            continue;
        end if;
        smooth +:= 1;

        C0 := HyperellipticCurve(f);
        J0 := Jacobian(C0);
        D9 := J0![x - 1, Evaluate(h, Q!1)];
        if Order(D9) ne 9 then
            continue;
        end if;
        marked9 +:= 1;

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
        G, phi := TorsionSubgroup(Jacobian(HyperellipticCurve(fI)));
        invs := Invariants(G);
        exact_tests +:= 1;
        key := Sprint(invs);
        if IsDefined(torsion_counts, key) then
            torsion_counts[key] +:= 1;
        else
            torsion_counts[key] := 1;
        end if;

        if HasInvariantDivisibleBy(invs, target) or TorsionOrder(invs) mod target eq 0 then
            hits +:= 1;
            simple, pcert, Lp := SimpleCertificate(fI);
            print "HIT45", "a", a, "torsion", invs,
                  "order", TorsionOrder(invs),
                  "simple", simple, "pcert", pcert;
            print "  f =", fI;
            print "  used", used;
        else
            print "SURVIVOR_NOT_45", "a", a, "torsion", invs, "used", used;
        end if;
    end for;

    print "DONE contact9 plus 5";
    print "checked", checked, "smooth", smooth, "marked9", marked9,
          "target_survivors", target_survivors, "exact_tests", exact_tests,
          "hits", hits;
    print "FIRST_KILL";
    for p in Sort([ k : k in Keys(first_kill) ]) do
        print " ", p, first_kill[p];
    end for;
    print "TORSION_COUNTS";
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
