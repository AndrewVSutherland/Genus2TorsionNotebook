//////////////////////////////////////////////////////////////////////
//  Refined boundary search for M(2,2,2,2) plus possible 7-torsion.
//
//  The open full-split chart is impossible at p=5 and p=13 for
//  rational 7-torsion.  After forcing boundary at both primes, the
//  height-12 search was mostly killed by the next good primes, especially
//  p=11 and p=17.  This script separates those two primes explicitly:
//
//      1. require boundary at p=5 and p=13;
//      2. classify p=11 and p=17 as good-pass, good-kill, or boundary;
//      3. keep only triples that are not killed at 11 or 17;
//      4. test the remaining good primes and run exact torsion if needed.
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

Q := Rationals();
Qx<x> := PolynomialRing(Q);

post_primes := [19,23,29,31,37,41,43,47,53,59,61,67,71];

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

function IsBoundaryAtPrime(a, b, c, p)
    vals := [a,b,c];
    finite := [];
    residues := [];
    for q in vals do
        if IsInfinityModP(q, p) then
            return true;
        end if;
        Append(~finite, true);
        Append(~residues, GF(p)!q);
    end for;
    for r in residues do
        if r eq 0 or r eq 1 then
            return true;
        end if;
    end for;
    return #(Set(residues)) lt 3;
end function;

function ReductionStatus(f, p)
    try
        fp := ChangeRing(f, GF(p));
    catch e
        return "boundary", false, false;
    end try;
    if not GoodHyperellipticPolynomial(fp) then
        return "boundary", false, false;
    end if;
    C := HyperellipticCurve(fp);
    pass7 := (#Jacobian(C) mod 7) eq 0;
    return pass7 select "good_pass7" else "good_kill7", true, pass7;
end function;

procedure Increment(~A, key)
    if not IsDefined(A, key) then
        A[key] := 0;
    end if;
    A[key] +:= 1;
end procedure;

params0 := RationalParametersOfHeight(height);
params := [u : u in params0 | u ne 0 and u ne 1];

print "M(2,2,2,2)+7 refined boundary search at 5,13 then 11,17";
print "height", height, "params", #params;
print "post_primes", post_primes;

checked := 0;
smooth := 0;
boundary_5_13 := 0;
compatible_11 := 0;
compatible_11_17 := 0;
post_survivors := 0;
exact_tests := 0;
hits := [];

status_counts := AssociativeArray();
pair_status_counts := AssociativeArray();
first_kill_counts := AssociativeArray(Integers());
for p in post_primes do
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

            if not (IsBoundaryAtPrime(a,b,c,5) and IsBoundaryAtPrime(a,b,c,13)) then
                continue;
            end if;
            boundary_5_13 +:= 1;

            status11, good11, pass11 := ReductionStatus(f, 11);
            Increment(~status_counts, "11:" cat status11);
            if good11 and not pass11 then
                continue;
            end if;
            compatible_11 +:= 1;

            status17, good17, pass17 := ReductionStatus(f, 17);
            Increment(~status_counts, "17:" cat status17);
            pair_status := "11:" cat status11 cat " | 17:" cat status17;
            Increment(~pair_status_counts, pair_status);
            if good17 and not pass17 then
                continue;
            end if;
            compatible_11_17 +:= 1;

            first_fail := 0;
            for p in post_primes do
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

            post_survivors +:= 1;
            fI, L := IntegralModelPolynomial(f);
            C := HyperellipticCurve(fI);
            J := Jacobian(C);
            G, phi := TorsionSubgroup(J);
            invs := Invariants(G);
            exact_tests +:= 1;
            print "POST_SURVIVOR", [a,b,c], pair_status, "torsion", invs, "f", fI;
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
print "boundary_5_13", boundary_5_13;
print "compatible_11", compatible_11;
print "compatible_11_17", compatible_11_17;
print "post_survivors", post_survivors;
print "exact_tests", exact_tests;
print "hits", #hits;

print "STATUS_COUNTS";
for key in Sort(SetToSequence(Keys(status_counts))) do
    print key, status_counts[key];
end for;

print "PAIR_STATUS_COUNTS_AFTER_11_COMPATIBLE";
for key in Sort(SetToSequence(Keys(pair_status_counts))) do
    print key, pair_status_counts[key];
end for;

print "FIRST_KILL_AFTER_11_17";
for p in post_primes do
    print p, first_kill_counts[p];
end for;
