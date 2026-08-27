//////////////////////////////////////////////////////////////////////
//  Larger-torsion search in the contact-5/order-20 family.
//
//  The family is
//
//      h = 1 + t*x + ((t^2 - 1)/2)*x^2,
//      f = h^2 - ((t + 1)^4/4)*x^5.
//
//  It always has a rational 5-torsion class and a rational order-4
//  class H halving the Weierstrass class at x=1, so generically it has
//  a point of order 20.  This script searches for larger torsion by:
//
//      1. factoring the residual quartic f/(x-1), hence finding extra
//         rational 2-torsion;
//      2. exact-computing torsion for reducible cases;
//      3. reporting only torsion order > 40, plus factor-type counts.
//
//  Typical run:
//      magma -b height:=1000 code/contact5_large_torsion_search.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 1000;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned max_exact then
    max_exact := 100000;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned progress_interval then
    progress_interval := 100000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Q := Rationals();
P<x> := PolynomialRing(Q);

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            t := Q!num/den;
            key := Sprint(t);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, t);
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

function FamilyPolynomial(t)
    b := (t^2 - 1)/2;
    h := 1 + t*x + b*x^2;
    f := h^2 - ((t + 1)^4/4)*x^5;
    return f, h, b;
end function;

function FactorTypeString(fac)
    degs := Sort([ Degree(ff[1]) : ff in fac ]);
    return Join([ IntegerToString(d) : d in degs ], "+");
end function;

function TorsionOrder(invs)
    if #invs eq 0 then
        return 1;
    end if;
    return &*invs;
end function;

function SimpleCertificate(f)
    C := HyperellipticCurve(f);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97] do
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

params := RationalParametersOfHeight(height);
type_counts := AssociativeArray();
torsion_counts := AssociativeArray();

checked := 0;
smooth := 0;
reducible := 0;
exact_tests := 0;
large_hits := [];

print "Contact-5/order-20 larger torsion search";
print "height", height, "parameters", #params, "max_exact", max_exact;

for t in params do
    if t in {Q!-1, Q!-3} then
        continue;
    end if;

    f, h, b := FamilyPolynomial(t);
    checked +:= 1;
    if progress_interval gt 0 and checked mod progress_interval eq 0 then
        print "progress", "checked", checked, "smooth", smooth,
              "reducible", reducible, "exact_tests", exact_tests,
              "large_hits", #large_hits;
    end if;

    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        continue;
    end if;
    smooth +:= 1;

    q := ExactQuotient(f, x - 1);
    facq := Factorization(q);
    ftype := FactorTypeString(facq);
    if IsDefined(type_counts, ftype) then
        type_counts[ftype] +:= 1;
    else
        type_counts[ftype] := 1;
    end if;

    if #facq eq 1 then
        continue;
    end if;
    reducible +:= 1;
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

    ord := TorsionOrder(invs);
    if ord gt 40 then
        simple, pcert, Lp := SimpleCertificate(fI);
        Append(~large_hits, <t,b,invs,ftype,simple,pcert,fI,facq>);
        print "LARGE", "t", t, "b", b, "torsion", invs, "order", ord,
              "factor_type", ftype, "simple", simple, "pcert", pcert;
        print "  f =", fI;
        print "  quartic_factorization =", facq;
    end if;
end for;

print "DONE height", height;
print "checked", checked, "smooth", smooth, "reducible", reducible,
      "exact_tests", exact_tests, "large_hits", #large_hits;
print "Factor type counts";
for key in Sort([ k : k in Keys(type_counts) ]) do
    print " ", key, type_counts[key];
end for;
print "Torsion counts among exact reducible cases";
for key in Sort([ k : k in Keys(torsion_counts) ]) do
    print " ", key, torsion_counts[key];
end for;

quit;
