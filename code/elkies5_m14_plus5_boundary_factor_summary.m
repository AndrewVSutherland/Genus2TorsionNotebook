//////////////////////////////////////////////////////////////////////
//  Factor-level boundary summary for M_1(8,4) plus 5.
//
//  The cleared discriminant of the n=1 family factors as
//
//    R^8 w^4 (R-1)^14 (R+1)^2 (w-1)^19 (w+1)^19
//    (w-R)(R+w)
//    (Rw-3R+3w-1)(Rw+3R+3w+1)
//    (R^4-2R^3+R^2w^2-R^2+2Rw^2-w^2)
//    (-2R^2+Rw^2-R+2w^2)^8.
//
//  This script labels which factors occur in the boundary reductions
//  at p=7,11,19 and records which other primes kill the candidates
//  that are boundary/non-killing at all three obstructing primes.
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 50;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

Q := Rationals();
Z := Integers();
Qx<x> := PolynomialRing(Q);
branch_primes := [7,11,19];
filter_primes := [3,7,11,13,17,19,23,29,31,37,41,43];

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            r := Q!num/den;
            key := Sprint(r);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, r);
            end if;
        end for;
    end for;
    return vals;
end function;

function IsSquareQ(q)
    q := Q!q;
    if q lt 0 then
        return false;
    end if;
    return IsSquare(Numerator(q)) and IsSquare(Denominator(q));
end function;

function GoodHyperellipticPolynomial(f)
    return Degree(f) in {5,6} and Discriminant(f) ne 0;
end function;

function M14FamilyPolynomial(R, w)
    m := R;
    n := Q!1;
    t := (2*m^2 + (1-w^2)*m*n - 2*w^2*n^2)/(4*(w^2-1));
    A := n^4*x^2
         + (m^3*n + 4*m^2*t + m*n^3 - 8*m*n*t + 4*n^2*t)*x
         + m^4;
    B := (m*n + 2*n^2 + 4*t)*x^2
         + (m^2 + 4*m*n + n^2 + 8*t)*x
         + (2*m^2 + m*n + 4*t);
    return x*A*B, t;
end function;

function PlusDisc(R, w)
    return -4*(w-R)*(R-1)^2*(R+1)*(R+w)*(w+1)
            *(R*w - 3*R + 3*w - 1);
end function;

function MinusDisc(R, w)
    return 4*(w-1)*(R+1)*(R*w + 3*R + 3*w + 1)
           *(R^4 - 2*R^3 + R^2*w^2 - R^2 + 2*R*w^2 - w^2);
end function;

function M14CoverPossible(R, w)
    return IsSquareQ(PlusDisc(R,w)) or IsSquareQ(MinusDisc(R,w));
end function;

procedure Increment(~A, key)
    if not IsDefined(A, key) then
        A[key] := 0;
    end if;
    A[key] +:= 1;
end procedure;

function JoinOrNone(tags)
    if #tags eq 0 then
        return "none";
    end if;
    return Join(tags, "&");
end function;

function FactorTagsAtPrime(R, w, p)
    tags := [];
    if Valuation(Denominator(R), p) gt 0 then Append(~tags, "R=inf"); end if;
    if Valuation(Denominator(w), p) gt 0 then Append(~tags, "w=inf"); end if;
    if #tags gt 0 then
        return tags;
    end if;

    K := GF(p);
    r := K!R;
    u := K!w;

    if r eq 0 then Append(~tags, "R=0"); end if;
    if u eq 0 then Append(~tags, "w=0"); end if;
    if r eq 1 then Append(~tags, "R=1"); end if;
    if r eq -1 then Append(~tags, "R=-1"); end if;
    if u eq 1 then Append(~tags, "w=1"); end if;
    if u eq -1 then Append(~tags, "w=-1"); end if;
    if u-r eq 0 then Append(~tags, "w-R"); end if;
    if r+u eq 0 then Append(~tags, "R+w"); end if;
    if r*u - 3*r + 3*u - 1 eq 0 then Append(~tags, "Lplus"); end if;
    if r*u + 3*r + 3*u + 1 eq 0 then Append(~tags, "Lminus"); end if;
    if r^4 - 2*r^3 + r^2*u^2 - r^2 + 2*r*u^2 - u^2 eq 0 then
        Append(~tags, "Qminus");
    end if;
    if -2*r^2 + r*u^2 - r + 2*u^2 eq 0 then
        Append(~tags, "Eresultant");
    end if;
    return tags;
end function;

function ReductionStatus(f, p)
    try
        fp := ChangeRing(f, GF(p));
    catch e
        return "bad_coeff", false, false;
    end try;
    if not GoodHyperellipticPolynomial(fp) then
        return "bad_curve", false, false;
    end if;
    C := HyperellipticCurve(fp);
    pass5 := (#Jacobian(C) mod 5) eq 0;
    return pass5 select "good_pass5" else "good_kill5", true, pass5;
end function;

params := RationalParametersOfHeight(height);
print "M_1(8,4) plus 5 factor-level boundary summary";
print "height", height, "params", #params;

checked := 0;
cover := 0;
smooth := 0;
branch_ok := 0;
global_ok := 0;

branch_factor_counts := AssociativeArray();
branch_survivor_factor_counts := AssociativeArray();
first_kill_after_branch := AssociativeArray();

for R in params do
    for w in params do
        if R eq 0 or w in {Q!-1, Q!0, Q!1} then
            continue;
        end if;
        checked +:= 1;
        if not M14CoverPossible(R,w) then
            continue;
        end if;
        cover +:= 1;
        f, t := M14FamilyPolynomial(R,w);
        if not GoodHyperellipticPolynomial(f) then
            continue;
        end if;
        smooth +:= 1;

        ok_branch := true;
        for p in branch_primes do
            status, good, pass5 := ReductionStatus(f, p);
            factags := JoinOrNone(FactorTagsAtPrime(R,w,p));
            key := Sprintf("%o:%o:%o", p, status, factags);
            Increment(~branch_factor_counts, key);
            if good and not pass5 then
                ok_branch := false;
            end if;
        end for;

        if not ok_branch then
            continue;
        end if;
        branch_ok +:= 1;

        for p in branch_primes do
            status, good, pass5 := ReductionStatus(f, p);
            factags := JoinOrNone(FactorTagsAtPrime(R,w,p));
            key := Sprintf("%o:%o:%o", p, status, factags);
            Increment(~branch_survivor_factor_counts, key);
        end for;

        first_fail := 0;
        for p in filter_primes do
            status, good, pass5 := ReductionStatus(f, p);
            if good and not pass5 then
                first_fail := p;
                break;
            end if;
        end for;
        if first_fail eq 0 then
            global_ok +:= 1;
        else
            Increment(~first_kill_after_branch, Sprint(first_fail));
        end if;
    end for;
end for;

print "checked", checked;
print "cover", cover;
print "smooth", smooth;
print "branch_ok_7_11_19", branch_ok;
print "global_ok_all_filter_primes", global_ok;

print "BRANCH_FACTOR_COUNTS_ALL";
for key in Sort(SetToSequence(Keys(branch_factor_counts))) do
    print key, branch_factor_counts[key];
end for;

print "BRANCH_FACTOR_COUNTS_FOR_7_11_19_SURVIVORS";
for key in Sort(SetToSequence(Keys(branch_survivor_factor_counts))) do
    print key, branch_survivor_factor_counts[key];
end for;

print "FIRST_KILL_AFTER_7_11_19_BRANCH";
for key in Sort(SetToSequence(Keys(first_kill_after_branch))) do
    print key, first_kill_after_branch[key];
end for;

print "DONE";
