//////////////////////////////////////////////////////////////////////
//  Boundary-at-3 search for contact-7 plus possible rational 5-torsion.
//
//  The open contact-7 surface has no good F_3 point with
//      5 | #J(F_3).
//  Thus a rational 35-torsion example in this family must reduce to
//  the bad boundary at p=3.  This script classifies the p=3 boundary
//  type and tests the remaining good primes for the necessary 5-torsion
//  condition.
//
//  Typical runs:
//      magma -b height:=20 code/contact7_plus5_boundary3_search.m
//      magma -b height:=30 only_tag:="finite:a=1,b=1,Q5=0" \
//          code/contact7_plus5_boundary3_search.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 20;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned max_exact then
    max_exact := 500;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned progress_interval then
    progress_interval := 50000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;
if not assigned only_tag then
    only_tag := "";
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

filter_primes := [7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97];

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

function Contact7Polynomial(a, b)
    h := 1 - (Q!7/2)*x + a*x^2 + b*x^3;
    return ExactQuotient(h^2 + (x - 1)^7, x^2), h;
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
    for p in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149] do
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

function Q5ComponentMod3(aa, bb)
    F := GF(3);
    a := F!aa; b := F!bb;
    q := 432*a^4 - 64*a^3*b^3 + 1008*a^3*b + 3024*a^3
        - 448*a^2*b^3 + 224*a^2*b^2 + 21168*a^2*b - 32536*a^2
        - 2016*a*b^4 + 4480*a*b^3 + 38416*a*b^2 - 109760*a*b
        + 78890*a - 864*b^5 + 5936*b^4 + 7056*b^3
        - 96040*b^2 + 120050*b - 60025;
    return q eq 0;
end function;

function BoundaryTagAt3(a, b, f, h)
    if Valuation(Denominator(a), 3) gt 0 or Valuation(Denominator(b), 3) gt 0 then
        tags := [];
        if Valuation(Denominator(a), 3) gt 0 then
            Append(~tags, "a=inf");
        else
            Append(~tags, Sprintf("a=%o", Integers()!(GF(3)!a)));
        end if;
        if Valuation(Denominator(b), 3) gt 0 then
            Append(~tags, "b=inf");
        else
            Append(~tags, Sprintf("b=%o", Integers()!(GF(3)!b)));
        end if;
        return true, "infinity:" cat Join(tags, ",");
    end if;

    F := GF(3);
    aa := Integers()!(F!a);
    bb := Integers()!(F!b);
    fp := ChangeRing(f, F);
    h1 := Evaluate(ChangeRing(h, F), F!1);
    good := Degree(fp) eq 5 and Discriminant(fp) ne 0 and h1 ne 0;
    if good then
        return false, "good";
    end if;

    tags := [Sprintf("a=%o", aa), Sprintf("b=%o", bb)];
    if h1 eq 0 then
        Append(~tags, "h1=0");
    end if;
    if F!(2*aa + 2*bb - 5) eq 0 then
        Append(~tags, "L=0");
    end if;
    if Q5ComponentMod3(aa, bb) then
        Append(~tags, "Q5=0");
    end if;
    return true, "finite:" cat Join(tags, ",");
end function;

function PassesFiveAwayFrom3(f)
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
            if n mod 5 ne 0 then
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

params := RationalParametersOfHeight(height);

print "contact-7 plus 5: p=3 boundary-focused search";
print "height", height, "params", #params, "only_tag", only_tag;
print "filter_primes", filter_primes;

checked := 0;
smooth := 0;
boundary := 0;
tag_counts := AssociativeArray();
tag_first_kill := AssociativeArray();
survivors := 0;
exact_tests := 0;
hits := 0;

for a in params do
    for b in params do
        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "progress", "checked", checked, "smooth", smooth,
                  "boundary", boundary, "survivors", survivors,
                  "exact_tests", exact_tests, "hits", hits;
        end if;

        f, h := Contact7Polynomial(a, b);
        if not GoodHyperellipticPolynomial(f) or Evaluate(h, Q!1) eq 0 then
            continue;
        end if;
        smooth +:= 1;

        is_boundary, tag := BoundaryTagAt3(a, b, f, h);
        if not is_boundary then
            continue;
        end if;
        if only_tag ne "" and tag ne only_tag then
            continue;
        end if;
        boundary +:= 1;
        Increment(~tag_counts, tag);

        pass, pbad, nbad, used := PassesFiveAwayFrom3(f);
        if not pass then
            Increment(~tag_first_kill, tag cat " | " cat IntegerToString(pbad));
            continue;
        end if;

        survivors +:= 1;
        if exact_tests ge max_exact then
            continue;
        end if;

        fI, L := IntegralModel(f);
        C := HyperellipticCurve(fI);
        J := Jacobian(C);
        G, phi := TorsionSubgroup(J);
        invs := Invariants(G);
        exact_tests +:= 1;
        print "SURVIVOR", "a", a, "b", b, "tag", tag,
              "torsion", invs, "used", used;
        if HasInvariantDivisibleBy(invs, 35) then
            hits +:= 1;
            simple, pcert, Lp := SimpleCertificate(fI);
            print "HIT35", "a", a, "b", b, "tag", tag,
                  "torsion", invs, "order", TorsionOrder(invs),
                  "simple", simple, "pcert", pcert;
            print "  f =", fI;
        end if;
    end for;
end for;

print "DONE";
print "checked", checked;
print "smooth", smooth;
print "boundary", boundary;
print "survivors", survivors;
print "exact_tests", exact_tests;
print "hits", hits;

print "TAG_COUNTS";
for key in Sort([ k : k in Keys(tag_counts) ]) do
    print key, tag_counts[key];
end for;

print "TAG_FIRST_KILL";
for key in Sort([ k : k in Keys(tag_first_kill) ]) do
    print key, tag_first_kill[key];
end for;

quit;
