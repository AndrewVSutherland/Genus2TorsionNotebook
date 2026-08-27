//////////////////////////////////////////////////////////////////////
//  Boundary-focused search for M(2,2,2,2) plus possible 7-torsion.
//
//  In the normalized full-split chart
//
//      y^2 = x*(x-1)*(x-a)*(x-b)*(x-c),
//
//  the good affine finite checks show no 7-possible points at p=5 and
//  p=13.  Therefore a rational example in this chart must be on the
//  boundary at both primes: one of a,b,c is 0, 1, infinity, or two of
//  a,b,c collide modulo p.
//
//  This script filters to those simultaneous boundary triples and then
//  tests the remaining good primes for the necessary condition
//
//      7 | #J(F_p).
//
//  It also records the boundary collision signatures at p=5 and p=13.
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 12;
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
if not assigned combo_min_count then
    combo_min_count := 1000000;
elif Type(combo_min_count) eq MonStgElt then
    combo_min_count := StringToInteger(combo_min_count);
end if;

Q := Rationals();
Z := Integers();
Qx<x> := PolynomialRing(Q);

boundary_primes := [5,13];
filter_primes := [11,17,19,23,29,31,37,41,43,47,53,59,61,67,71];

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            r := Q!num/Q!den;
            key := Sprint(r);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, r);
            end if;
        end for;
    end for;
    return vals;
end function;

function FullSplitPolynomial(a, b, c)
    return x*(x-1)*(x-a)*(x-b)*(x-c);
end function;

function GoodHyperellipticPolynomial(f)
    return Degree(f) in {5,6} and Discriminant(f) ne 0;
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return Qx!(L^2*f), L;
end function;

function HasSevenInvariants(invs)
    return &or [ n mod 7 eq 0 : n in invs ];
end function;

function IsInfinityModP(q, p)
    return Valuation(Denominator(q), p) gt 0;
end function;

function ResidueModP(q, p)
    return GF(p)!q;
end function;

function BoundaryTagsAtPrime(a, b, c, p)
    tags := [];
    vals := [a,b,c];
    names := ["a","b","c"];

    finite := [];
    residues := [];
    for i in [1..3] do
        q := vals[i];
        name := names[i];
        if IsInfinityModP(q, p) then
            Append(~tags, name cat "=inf");
            Append(~finite, false);
            Append(~residues, GF(p)!0);
        else
            r := ResidueModP(q, p);
            Append(~finite, true);
            Append(~residues, r);
            if r eq 0 then
                Append(~tags, name cat "=0");
            end if;
            if r eq 1 then
                Append(~tags, name cat "=1");
            end if;
        end if;
    end for;

    for i in [1..3] do
        for j in [i+1..3] do
            if finite[i] and finite[j] and residues[i] eq residues[j] then
                Append(~tags, names[i] cat "=" cat names[j]);
            end if;
        end for;
    end for;

    if #tags eq 0 then
        return false, "good";
    end if;
    return true, Join(Sort(tags), "&");
end function;

function ReductionStatus(f, p)
    try
        fp := ChangeRing(f, GF(p));
    catch e
        return "bad", false, false;
    end try;
    if not GoodHyperellipticPolynomial(fp) then
        return "bad", false, false;
    end if;
    C := HyperellipticCurve(fp);
    pass7 := (#Jacobian(C) mod 7) eq 0;
    return pass7 select "pass7" else "kill7", true, pass7;
end function;

procedure Increment(~A, key)
    if not IsDefined(A, key) then
        A[key] := 0;
    end if;
    A[key] +:= 1;
end procedure;

params0 := RationalParametersOfHeight(height);
params := [u : u in params0 | u ne 0 and u ne 1];

print "M(2,2,2,2)+7 boundary-focused search";
print "height", height, "params", #params;
print "boundary_primes", boundary_primes, "filter_primes", filter_primes;
print "combo_min_count", combo_min_count;

checked := 0;
smooth := 0;
boundary_5 := 0;
boundary_5_13 := 0;
reduction_survivors := 0;
exact_tests := 0;
hits := [];

tag_counts := AssociativeArray();
combo_counts := AssociativeArray();
survivor_combo_counts := AssociativeArray();
first_kill_counts := AssociativeArray(Integers());
for p in filter_primes do
    first_kill_counts[p] := 0;
end for;

for i in [1..#params] do
    if #hits ge max_hits then break; end if;
    for j in [i+1..#params] do
        if #hits ge max_hits then break; end if;
        for k in [j+1..#params] do
            a := params[i]; b := params[j]; c := params[k];
            if #(Set([a,b,c,Q!0,Q!1])) ne 5 then
                continue;
            end if;
            checked +:= 1;
            f := FullSplitPolynomial(a,b,c);
            if not GoodHyperellipticPolynomial(f) then
                continue;
            end if;
            smooth +:= 1;

            is5, tag5 := BoundaryTagsAtPrime(a,b,c,5);
            is13, tag13 := BoundaryTagsAtPrime(a,b,c,13);
            if is5 then
                boundary_5 +:= 1;
            end if;
            if not (is5 and is13) then
                continue;
            end if;
            boundary_5_13 +:= 1;
            Increment(~tag_counts, "5:" cat tag5);
            Increment(~tag_counts, "13:" cat tag13);
            combo := "5:" cat tag5 cat " | 13:" cat tag13;
            Increment(~combo_counts, combo);

            first_fail := 0;
            for p in filter_primes do
                status, good, pass7 := ReductionStatus(f, p);
                if good and not pass7 then
                    first_fail := p;
                    break;
                end if;
            end for;
            if first_fail ne 0 then
                first_kill_counts[first_fail] +:= 1;
                continue;
            end if;

            reduction_survivors +:= 1;
            Increment(~survivor_combo_counts, combo);
            fI, L := IntegralModelPolynomial(f);
            C := HyperellipticCurve(fI);
            J := Jacobian(C);
            G, phi := TorsionSubgroup(J);
            invs := Invariants(G);
            exact_tests +:= 1;
            print "REDUCTION_SURVIVOR", [a,b,c], "tags", combo,
                  "torsion", invs, "f", fI;
            if HasSevenInvariants(invs) then
                Append(~hits, <a,b,c,invs,fI>);
                print "HIT", [a,b,c], "torsion", invs, "f", fI;
                if #hits ge max_hits then
                    break;
                end if;
            end if;
            if exact_tests ge max_tests then
                break i;
            end if;
        end for;
    end for;
end for;

print "DONE";
print "checked", checked;
print "smooth", smooth;
print "boundary_5", boundary_5;
print "boundary_5_13", boundary_5_13;
print "reduction_survivors", reduction_survivors;
print "exact_tests", exact_tests;
print "hits", #hits;

print "FIRST_KILL_AFTER_BOUNDARY";
for p in filter_primes do
    print p, first_kill_counts[p];
end for;

print "BOUNDARY_TAG_COUNTS";
for key in Sort(SetToSequence(Keys(tag_counts))) do
    print key, tag_counts[key];
end for;

print "BOUNDARY_COMBO_COUNTS";
for key in Sort(SetToSequence(Keys(combo_counts))) do
    if combo_counts[key] ge combo_min_count then
        print key, combo_counts[key];
    end if;
end for;

print "SURVIVOR_COMBO_COUNTS";
for key in Sort(SetToSequence(Keys(survivor_combo_counts))) do
    print key, survivor_combo_counts[key];
end for;
