//////////////////////////////////////////////////////////////////////
//  p=3 boundary diagnostic for the contact-7 rational-root plus 3 route.
//
//  The finite (s,b,eps) root chart has no good open F_3 points.  This
//  script separates the rational height search by 3-adic boundary branch
//  and then applies the necessary prime-to-p target-42 reductions away
//  from p=3.
//
//  Typical runs:
//      magma -b mode:=finite code/contact7_root_plus3_boundary3.m
//      magma -b mode:=search height:=20 code/contact7_root_plus3_boundary3.m
//////////////////////////////////////////////////////////////////////

if not assigned mode then
    mode := "finite";
end if;
if not assigned height then
    height := 20;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned progress_interval then
    progress_interval := 50000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Q := Rationals();
P<x> := PolynomialRing(Q);

target := 42;
filter_primes := [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97];

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

function VpQ(q, p)
    if q eq 0 then
        return 999;
    end if;
    return Valuation(Numerator(q), p) - Valuation(Denominator(q), p);
end function;

function IsPIntegral(q, p)
    return Valuation(Denominator(q), p) eq 0;
end function;

function ResidueString(q, p)
    if not IsPIntegral(q, p) then
        return "inf";
    end if;
    return IntegerToString(Integers()!(GF(p)!q));
end function;

function PolynomialIntegralAtP(f, p)
    for i in [0..Degree(f)] do
        if not IsPIntegral(Coefficient(f, i), p) then
            return false;
        end if;
    end for;
    return true;
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

function GoodHyperellipticPolynomial(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function RequiredPrimeToP(n, p)
    required := n;
    while required mod p eq 0 do
        required div:= p;
    end while;
    return required;
end function;

function ReductionAt3Tag(f, h)
    if not PolynomialIntegralAtP(f, 3) or not PolynomialIntegralAtP(h, 3) then
        return "red3=nonintegral", false, 0;
    end if;

    F := GF(3);
    fp := ChangeRing(f, F);
    hp := ChangeRing(h, F);
    h1 := Evaluate(hp, F!1);

    if Degree(fp) ne 5 then
        return Sprintf("red3=bad:deg%o:h1%o", Degree(fp), Integers()!h1), false, 0;
    end if;

    disc := Discriminant(fp);
    if disc eq 0 or h1 eq 0 then
        return Sprintf("red3=bad:h1%o:disc%o", Integers()!h1, Integers()!disc), false, 0;
    end if;

    C := HyperellipticCurve(fp);
    n := Integers()!#Jacobian(C);
    return Sprintf("red3=good:#J%o:pass14%o", n, n mod 14 eq 0), true, n;
end function;

function BranchAt3Tag(s, b, eps, a, r)
    parts := [];
    vs := VpQ(s, 3);
    vr := VpQ(r, 3);
    va := VpQ(a, 3);
    vb := VpQ(b, 3);
    ves := VpQ(eps*s, 3);

    if vs lt 0 then
        Append(~parts, "s=inf");
    elif vr gt 0 then
        Append(~parts, "r=0");
        if ves ge 0 and ResidueString(eps*s, 3) eq "1" then
            Append(~parts, "cancel=yes");
        else
            Append(~parts, "cancel=no");
        end if;
    elif vs gt 0 then
        Append(~parts, "s=0");
    else
        Append(~parts, "s=unit");
    end if;

    Append(~parts, "sres=" cat ResidueString(s, 3));
    Append(~parts, "epss=" cat ResidueString(eps*s, 3));
    Append(~parts, "ares=" cat ResidueString(a, 3));
    Append(~parts, "bres=" cat ResidueString(b, 3));

    if va lt 0 then
        Append(~parts, "a=inf");
    end if;
    if vb lt 0 then
        Append(~parts, "b=inf");
    end if;

    return Join(parts, ",");
end function;

function PassesTargetAwayFrom3(f)
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

procedure Increment(~A, key)
    if not IsDefined(A, key) then
        A[key] := 0;
    end if;
    A[key] +:= 1;
end procedure;

procedure RunFinite()
    F := GF(3);
    PF<X> := PolynomialRing(F);
    print "FINITE p=3 boundary summary for contact-7 rational-root plus 3";

    print "finite root-chart branch s=0, r=1";
    for b0 in F do
        h := 1 + X + (1 - b0)*X^2 + b0*X^3;
        f := ExactQuotient(h^2 + (X - 1)^7, X^2);
        print "  b", Integers()!b0, "h1", Integers()!Evaluate(h, F!1),
              "disc", Integers()!Discriminant(f),
              "factor", Factorization(f);
    end for;

    print "r=0 cancellation limit: eps*s == 1, a -> 1";
    for b0 in F do
        h := 1 + X + X^2 + b0*X^3;
        f := ExactQuotient(h^2 + (X - 1)^7, X^2);
        line := Sprintf("  b %o h1 %o disc %o factor %o",
                        Integers()!b0, Integers()!Evaluate(h, F!1),
                        Integers()!Discriminant(f), Factorization(f));
        if Degree(f) eq 5 and Discriminant(f) ne 0 then
            C := HyperellipticCurve(f);
            line cat:= Sprintf(" #J %o L %o", #Jacobian(C), LPolynomial(C));
        end if;
        print line;
    end for;
end procedure;

procedure RunSearch()
    params := RationalParametersOfHeight(height);
    checked := 0;
    smooth := 0;
    tag_counts := AssociativeArray();
    tag_good3 := AssociativeArray();
    tag_first_kill := AssociativeArray();
    tag_survivors := AssociativeArray();
    survivors := 0;

    print "SEARCH p=3 boundary tags for contact-7 rational-root plus 3";
    print "height", height, "params", #params, "target", target;
    print "filter_primes_away_from_3", filter_primes;

    for s in params do
        for b in params do
            for eps in [Q!-1, Q!1] do
                checked +:= 1;
                if progress_interval gt 0 and checked mod progress_interval eq 0 then
                    print "progress", "checked", checked, "smooth", smooth,
                          "survivors_away3", survivors;
                end if;

                ok, f, h, a, r := Contact7RootPolynomial(s, b, eps);
                if not ok or not GoodHyperellipticPolynomial(f) then
                    continue;
                end if;
                if Evaluate(h, Q!1) eq 0 then
                    continue;
                end if;
                smooth +:= 1;

                branch := BranchAt3Tag(s, b, eps, a, r);
                redtag, good3, n3 := ReductionAt3Tag(f, h);
                tag := branch cat " | " cat redtag;
                Increment(~tag_counts, tag);
                if good3 then
                    Increment(~tag_good3, tag);
                end if;

                pass, pbad, nbad, used := PassesTargetAwayFrom3(f);
                if not pass then
                    Increment(~tag_first_kill, tag cat " | firstkill=" cat IntegerToString(pbad));
                    continue;
                end if;

                survivors +:= 1;
                Increment(~tag_survivors, tag);
                print "SURVIVOR_AWAY3", "s", s, "b", b, "eps", eps,
                      "a", a, "r", r, "tag", tag, "used", used;
            end for;
        end for;
    end for;

    print "DONE";
    print "checked", checked, "smooth", smooth, "survivors_away3", survivors;
    print "TAG_COUNTS";
    for key in Sort([ k : k in Keys(tag_counts) ]) do
        print " ", key, tag_counts[key];
    end for;
    print "TAG_GOOD3";
    for key in Sort([ k : k in Keys(tag_good3) ]) do
        print " ", key, tag_good3[key];
    end for;
    print "TAG_SURVIVORS_AWAY3";
    for key in Sort([ k : k in Keys(tag_survivors) ]) do
        print " ", key, tag_survivors[key];
    end for;
    print "TAG_FIRST_KILL";
    for key in Sort([ k : k in Keys(tag_first_kill) ]) do
        print " ", key, tag_first_kill[key];
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
