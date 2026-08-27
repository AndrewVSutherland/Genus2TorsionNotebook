//////////////////////////////////////////////////////////////////////
//  Direct contact-9 family and rational-root subfamily.
//
//  On an odd genus-2 model y^2 = f(x), the function x^2*y - h(x)
//  has pole order 9 at infinity.  Force
//
//      h(x)^2 - x^4*f(x) = -(x - 1)^9.
//
//  Taking h congruent to (1-x)^(9/2) modulo x^4 gives
//
//      h = 1 - (9/2)*x + (63/8)*x^2 - (105/16)*x^3 + a*x^4,
//      f = (h^2 + (x - 1)^9)/x^4.
//
//  Then P=(1,h(1)) gives a rational divisor class of order 9, provided
//  the curve is smooth and h(1) != 0.
//
//  The rational-root subfamily sets r = 1-s^2 and h(r)=eps*s^9, eps=+-1,
//  which determines a and forces f(r)=0.  This gives rational 2-torsion
//  together with the marked 9-torsion class.
//
//  Typical runs:
//      magma -b mode:=verify code/contact9_family_search.m
//      magma -b mode:=family height:=20 code/contact9_family_search.m
//      magma -b mode:=root_finite target:=72 code/contact9_family_search.m
//      magma -b mode:=root height:=40 target:=72 code/contact9_family_search.m
//////////////////////////////////////////////////////////////////////

if not assigned mode then
    mode := "verify";
end if;
if not assigned height then
    height := 20;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned target then
    target := 72;
elif Type(target) eq MonStgElt then
    target := StringToInteger(target);
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

Q := Rationals();
P<x> := PolynomialRing(Q);

odd_primes := [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97];

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

function Contact9RootPolynomial(s, eps)
    r := 1 - s^2;
    if r eq 0 then
        return false, P!0, P!0, Q!0, Q!0;
    end if;
    h0 := 1 - (Q!9/2)*r + (Q!63/8)*r^2 - (Q!105/16)*r^3;
    a := (eps*s^9 - h0)/r^4;
    f, h := Contact9Polynomial(a);
    if Evaluate(f, r) ne 0 then
        return false, f, h, a, r;
    end if;
    return true, f, h, a, r;
end function;

function GoodContactCurve(f, h)
    return Degree(f) eq 5 and Discriminant(f) ne 0 and Evaluate(h, Q!1) ne 0;
end function;

function SimpleCertificate(f)
    C := HyperellipticCurve(f);
    for p in odd_primes do
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

function PassesTargetReduction(f, n)
    used := [];
    for p in odd_primes do
        try
            fp := ChangeRing(f, GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
                continue;
            end if;
            C := HyperellipticCurve(fp);
            jp := Integers()!#Jacobian(C);
            Append(~used, <p,jp>);
            required := RequiredPrimeToP(n, p);
            if required gt 1 and jp mod required ne 0 then
                return false, p, jp, used;
            end if;
        catch e
            continue;
        end try;
    end for;
    return true, 0, 0, used;
end function;

function Contact9PolynomialFinite(F, a)
    PF<X> := PolynomialRing(F);
    h := 1 - (F!9/F!2)*X + (F!63/F!8)*X^2 - (F!105/F!16)*X^3 + a*X^4;
    f := ExactQuotient(h^2 + (X - 1)^9, X^4);
    return f, h;
end function;

function Contact9RootPolynomialFinite(F, s, eps)
    PF<X> := PolynomialRing(F);
    r := 1 - s^2;
    if r eq 0 then
        return false, PF!0, PF!0, F!0, F!0;
    end if;
    h0 := 1 - (F!9/F!2)*r + (F!63/F!8)*r^2 - (F!105/F!16)*r^3;
    a := (eps*s^9 - h0)/r^4;
    f, h := Contact9PolynomialFinite(F, a);
    if Evaluate(f, r) ne 0 then
        return false, f, h, a, r;
    end if;
    return true, f, h, a, r;
end function;

procedure VerifySamples()
    samples := [Q!0, Q!1, -Q!1, Q!2, -Q!3/2, Q!5/4];
    print "VERIFY contact-9 samples";
    for a in samples do
        f, h := Contact9Polynomial(a);
        print "a", a, "f", f, "h1", Evaluate(h, Q!1), "disc_nonzero", Discriminant(f) ne 0;
        if GoodContactCurve(f, h) then
            C := HyperellipticCurve(f);
            J := Jacobian(C);
            D9 := J![x - 1, Evaluate(h, Q!1)];
            ord := Order(D9);
            fI, L := IntegralModel(f);
            G, phi := TorsionSubgroup(Jacobian(HyperellipticCurve(fI)));
            simple, pcert, Lp := SimpleCertificate(fI);
            print "  marked_order", ord, "torsion", Invariants(G),
                  "simple", simple, "pcert", pcert;
        end if;
    end for;
end procedure;

procedure FamilySearch()
    vals := RationalParametersOfHeight(height);
    checked := 0;
    smooth := 0;
    marked9 := 0;
    exact_tests := 0;
    hits := 0;
    torsion_counts := AssociativeArray();

    print "CONTACT9 direct family search";
    print "height", height, "params", #vals, "max_exact", max_exact;
    for a in vals do
        checked +:= 1;
        f, h := Contact9Polynomial(a);
        if not GoodContactCurve(f, h) then
            continue;
        end if;
        smooth +:= 1;
        C := HyperellipticCurve(f);
        J := Jacobian(C);
        D9 := J![x - 1, Evaluate(h, Q!1)];
        if Order(D9) ne 9 then
            continue;
        end if;
        marked9 +:= 1;
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
        if TorsionOrder(invs) ge 18 then
            hits +:= 1;
            simple, pcert, Lp := SimpleCertificate(fI);
            print "FAMILY_HIT", "a", a, "torsion", invs,
                  "order", TorsionOrder(invs), "simple", simple, "pcert", pcert;
            print "  f =", fI;
        end if;
    end for;

    print "DONE contact9 family";
    print "checked", checked, "smooth", smooth, "marked9", marked9,
          "exact_tests", exact_tests, "hits_ge18", hits;
    print "TORSION_COUNTS";
    for key in Sort([ k : k in Keys(torsion_counts) ]) do
        print " ", key, torsion_counts[key];
    end for;
end procedure;

procedure RootFinite()
    print "CONTACT9 rational-root finite diagnostic";
    print "target", target;
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43] do
        F := GF(p);
        total := 0;
        good := 0;
        pass_target := 0;
        samples := [];
        for s in F do
            for eps in [F!-1, F!1] do
                total +:= 1;
                ok, f, h, a, r := Contact9RootPolynomialFinite(F, s, eps);
                if not ok or Degree(f) ne 5 or Discriminant(f) eq 0 or Evaluate(h, F!1) eq 0 then
                    continue;
                end if;
                good +:= 1;
                C := HyperellipticCurve(f);
                jp := Integers()!#Jacobian(C);
                required := RequiredPrimeToP(target, p);
                if required eq 1 or jp mod required eq 0 then
                    pass_target +:= 1;
                    if #samples lt 10 then
                        Append(~samples, <Integers()!s, Integers()!eps, Integers()!a,
                                         Integers()!r, jp>);
                    end if;
                end if;
            end for;
        end for;
        print "p", p, "total", total, "good", good,
              "pass_target", pass_target, "samples", samples;
    end for;
end procedure;

procedure RootSearch()
    vals := RationalParametersOfHeight(height);
    checked := 0;
    smooth := 0;
    marked9 := 0;
    target_survivors := 0;
    exact_tests := 0;
    hits := 0;
    first_kill := AssociativeArray();
    torsion_counts := AssociativeArray();

    print "CONTACT9 rational-root search";
    print "height", height, "params", #vals, "target", target,
          "max_exact", max_exact;
    print "filter_primes", odd_primes;

    for s in vals do
        for eps in [Q!-1, Q!1] do
            checked +:= 1;
            if progress_interval gt 0 and checked mod progress_interval eq 0 then
                print "progress", "checked", checked, "smooth", smooth,
                      "marked9", marked9, "target_survivors", target_survivors,
                      "exact_tests", exact_tests, "hits", hits;
            end if;

            ok, f, h, a, r := Contact9RootPolynomial(s, eps);
            if not ok or not GoodContactCurve(f, h) then
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

            pass, pbad, nbad, used := PassesTargetReduction(f, target);
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
                print "ROOT_HIT", "s", s, "eps", eps, "a", a, "r", r,
                      "torsion", invs, "order", TorsionOrder(invs),
                      "simple", simple, "pcert", pcert;
                print "  f =", fI;
                print "  used", used;
            else
                print "ROOT_SURVIVOR_NOT_TARGET", "s", s, "eps", eps,
                      "a", a, "r", r, "torsion", invs, "used", used;
            end if;
        end for;
    end for;

    print "DONE contact9 root search";
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

if mode eq "verify" then
    VerifySamples();
elif mode eq "family" then
    FamilySearch();
elif mode eq "root_finite" then
    RootFinite();
elif mode eq "root" then
    RootSearch();
else
    print "unknown mode", mode;
end if;

quit;
