//////////////////////////////////////////////////////////////////////
//  Rational-root subfamily of the contact-7 family.
//
//  In contact7_family_search.m we use
//
//      h = 1 - (7/2)*x + a*x^2 + b*x^3,
//      f = (h^2 + (x - 1)^7)/x^2,
//
//  and P=(1,h(1)) gives a rational 7-torsion divisor.
//
//  To force a rational Weierstrass root x=r, write r=1-s^2.  Then
//      -(r-1)^7 = s^14,
//  so f(r)=0 is equivalent to h(r)=eps*s^7, eps=+/-1.
//  For r != 0 this is linear in a:
//
//      a = (eps*s^7 - 1 + (7/2)*r - b*r^3)/r^2.
//
//  This script searches that rational-root subfamily for torsion
//  divisible by a target such as 28 or 56.
//
//  Typical runs:
//      magma -b height:=12 target:=56 code/contact7_root_even_search.m
//      magma -b height:=12 target:=28 code/contact7_root_even_search.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 12;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned target then
    target := 56;
elif Type(target) eq MonStgElt then
    target := StringToInteger(target);
end if;
if not assigned max_exact then
    max_exact := 5000;
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

function PassesTargetReduction(f, target)
    used := [];
    for p in filter_primes do
        try
            fp := ChangeRing(f, GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
                continue;
            end if;
            C := HyperellipticCurve(fp);
            n := Integers()!#Jacobian(C);
            Append(~used, <p, n>);
            required := target;
            while required mod p eq 0 do
                required div:= p;
            end while;
            if required gt 1 and n mod required ne 0 then
                return false, p, n, used;
            end if;
        catch e
            continue;
        end try;
    end for;
    return true, 0, 0, used;
end function;

params := RationalParametersOfHeight(height);

checked := 0;
smooth := 0;
verified7 := 0;
target_survivors := 0;
exact_tests := 0;
hits := [];
first_kill := AssociativeArray();
torsion_counts := AssociativeArray();

print "CONTACT7 rational-root even search";
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
            smooth +:= 1;

            yP := Evaluate(h, Q!1);
            if yP eq 0 then
                continue;
            end if;

            if verify_all then
                C := HyperellipticCurve(f);
                J := Jacobian(C);
                D := J![x - 1, yP];
                if Order(D) ne 7 then
                    continue;
                end if;
            end if;
            verified7 +:= 1;

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
            CI := HyperellipticCurve(fI);
            JI := Jacobian(CI);
            G, phi := TorsionSubgroup(JI);
            invs := Invariants(G);
            exact_tests +:= 1;
            invkey := Sprint(invs);
            if IsDefined(torsion_counts, invkey) then
                torsion_counts[invkey] +:= 1;
            else
                torsion_counts[invkey] := 1;
            end if;

            ord := TorsionOrder(invs);
            if ord mod target eq 0 then
                simple, pcert, Lp := SimpleCertificate(fI);
                Append(~hits, <s,b,eps,a,r,invs,simple,pcert,fI>);
                print "HIT", "s", s, "b", b, "eps", eps, "a", a, "r", r,
                      "torsion", invs, "order", ord, "simple", simple,
                      "pcert", pcert;
                print "  f =", fI;
                print "  factorization =", Factorization(f);
            else
                print "SURVIVOR_NOT_TARGET", "s", s, "b", b, "eps", eps,
                      "a", a, "r", r, "torsion", invs, "used", used;
            end if;
        end for;
    end for;
end for;

print "DONE contact7 rational-root even search";
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

quit;
