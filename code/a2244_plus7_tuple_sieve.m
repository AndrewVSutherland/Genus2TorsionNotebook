//////////////////////////////////////////////////////////////////////
//  A(2,2,4,4) tuple file plus possible rational 7-torsion.
//
//  Input tuples [a,b,c,d] represent
//
//      C: y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2).
//
//  The file data/tor2244_bank.txt is an existing source of
//  A(2,2,4,4) examples.  This script applies the necessary condition
//  for rational 7-torsion: for every good prime p != 7, 7 | #J(F_p).
//////////////////////////////////////////////////////////////////////

if not assigned input_file then
    input_file := "data/tor2244_bank.txt";
end if;
if not assigned max_rows then
    max_rows := 0;
elif Type(max_rows) eq MonStgElt then
    max_rows := StringToInteger(max_rows);
end if;
if not assigned progress_interval then
    progress_interval := 5000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;
if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;
if not assigned exact_tests then
    exact_tests := true;
elif Type(exact_tests) eq MonStgElt then
    exact_tests := exact_tests eq "true" or exact_tests eq "1";
end if;

Q := Rationals();
P<x> := PolynomialRing(Q);
prime_list := [3,5,11,13,17,19,23,29,31,37,41,43,47,53,59,61];

function ReadTupleFile(filename)
    S := Read(filename);
    rows := Split(S, "\n");
    tuples := [];
    for rawrow in rows do
        n := #rawrow;
        if n gt 0 and rawrow[n] eq "\r" then
            n -:= 1;
        end if;
        if n ge 2 and rawrow[1] eq "[" and rawrow[n] eq "]" then
            body := rawrow[2..n-1];
            parts := Split(body, ",");
            tup := [ StringToInteger(part) : part in parts ];
            if #tup eq 4 then
                Append(~tuples, tup);
            end if;
        end if;
    end for;
    return tuples;
end function;

function TuplePolynomial(tup)
    f := x;
    for a in tup do
        f *:= x + (Q!a)^2;
    end for;
    return f;
end function;

function GoodHyperellipticPolynomial(f)
    return Degree(f) in {5,6} and Discriminant(f) ne 0;
end function;

function HasSevenInvariants(invs)
    return &or [ n mod 7 eq 0 : n in invs ];
end function;

function PassesSevenReduction(f, primes)
    for p in primes do
        try
            fp := ChangeRing(f, GF(p));
        catch e
            continue;
        end try;
        if not GoodHyperellipticPolynomial(fp) then
            continue;
        end if;
        C := HyperellipticCurve(fp);
        if (#Jacobian(C) mod 7) ne 0 then
            return false, p;
        end if;
    end for;
    return true, 0;
end function;

tuples := ReadTupleFile(input_file);
if max_rows gt 0 and max_rows lt #tuples then
    tuples := tuples[1..max_rows];
end if;

print "A(2,2,4,4) plus 7 tuple sieve";
print "input_file", input_file;
print "rows", #tuples;
print "prime_list", prime_list;
print "exact_tests", exact_tests;

checked := 0;
smooth := 0;
survivors := 0;
exact := 0;
hits := 0;
first_kill_counts := AssociativeArray(Integers());
for p in prime_list do
    first_kill_counts[p] := 0;
end for;

for tup in tuples do
    checked +:= 1;
    if progress_interval gt 0 and checked mod progress_interval eq 0 then
        print "PROGRESS", checked, "smooth", smooth, "survivors", survivors, "hits", hits;
    end if;
    f := TuplePolynomial(tup);
    if not GoodHyperellipticPolynomial(f) then
        continue;
    end if;
    smooth +:= 1;
    ok, badp := PassesSevenReduction(f, prime_list);
    if not ok then
        first_kill_counts[badp] +:= 1;
        continue;
    end if;
    survivors +:= 1;
    print "SEVEN_REDUCTION_SURVIVOR", tup;

    if exact_tests then
        C := HyperellipticCurve(f);
        J := Jacobian(C);
        G, phi := TorsionSubgroup(J);
        invs := Invariants(G);
        exact +:= 1;
        print "EXACT_TORSION", tup, invs;
        if HasSevenInvariants(invs) then
            hits +:= 1;
            print "HIT", tup, invs, "f", f;
            if hits ge max_hits then
                break;
            end if;
        end if;
    end if;
end for;

print "DONE";
print "checked", checked;
print "smooth", smooth;
print "survivors", survivors;
print "exact_tests_done", exact;
print "hits", hits;
for p in prime_list do
    print "FIRST_KILL", p, first_kill_counts[p];
end for;
