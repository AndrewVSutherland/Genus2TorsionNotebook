//////////////////////////////////////////////////////////////////////
//  Boundary-branch classifier for the M_1(8,4) plus 5 search.
//
//  The finite good affine chart for the [4,8] tangent-cover family has
//  no 5-divisible reductions at p=7 and p=11.  Therefore any rational
//  [4,8]+5 example in this chart must reduce to a boundary/bad stratum
//  at those primes.  This script classifies which boundary strata occur
//  among rational smooth tangent-cover candidates.
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 50;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

if not assigned max_reports then
    max_reports := 40;
elif Type(max_reports) eq MonStgElt then
    max_reports := StringToInteger(max_reports);
end if;

Q := Rationals();
Z := Integers();
Qx<x> := PolynomialRing(Q);
prime_list := [7,11,19];
all_filter_primes := [3,7,11,13,17,19,23,29,31,37,41,43];

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

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return Qx!(L^2*f), L;
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
    return x*A*B, t, A, B;
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

function ValuationQ(q, p)
    return Valuation(Numerator(q), p) - Valuation(Denominator(q), p);
end function;

function ResidueString(q, p)
    if Valuation(Denominator(q), p) gt 0 then
        return "inf";
    end if;
    K := GF(p);
    return Sprint(K!q);
end function;

procedure Increment(~A, key)
    if not IsDefined(A, key) then
        A[key] := 0;
    end if;
    A[key] +:= 1;
end procedure;

function ClassifyAtPrime(R, w, f, p)
    tags := [];
    vRnum := Valuation(Numerator(R), p);
    vRden := Valuation(Denominator(R), p);
    vwnum := Valuation(Numerator(w), p);
    vwden := Valuation(Denominator(w), p);
    vwminus := Valuation(Numerator(w - 1), p) - Valuation(Denominator(w - 1), p);
    vwplus := Valuation(Numerator(w + 1), p) - Valuation(Denominator(w + 1), p);

    if vRden gt 0 then Append(~tags, "R=inf"); end if;
    if vRnum gt 0 then Append(~tags, "R=0"); end if;
    if vwden gt 0 then Append(~tags, "w=inf"); end if;
    if vwnum gt 0 then Append(~tags, "w=0"); end if;
    if vwminus gt 0 then Append(~tags, "w=1"); end if;
    if vwplus gt 0 then Append(~tags, "w=-1"); end if;

    fp_ok := true;
    try
        fp := ChangeRing(f, GF(p));
    catch e
        fp_ok := false;
    end try;

    if not fp_ok then
        Append(~tags, "coeff_bad");
        return tags, false, false;
    end if;

    good := GoodHyperellipticPolynomial(fp);
    pass5 := false;
    if good then
        C := HyperellipticCurve(fp);
        pass5 := (#Jacobian(C) mod 5) eq 0;
        if pass5 then
            Append(~tags, "good_pass5");
        else
            Append(~tags, "good_kill5");
        end if;
    else
        if Degree(fp) notin {5,6} then
            Append(~tags, "degree_drop");
        end if;
        if Degree(fp) in {5,6} and Discriminant(fp) eq 0 then
            Append(~tags, "disc_zero");
        end if;
    end if;

    return tags, good, pass5;
end function;

function JoinTags(tags)
    if #tags eq 0 then
        return "none";
    end if;
    return Join(tags, "&");
end function;

params := RationalParametersOfHeight(height);

print "M_1(8,4) plus 5 boundary branch classifier";
print "height", height, "params", #params, "branch_primes", prime_list;

checked := 0;
cover := 0;
smooth := 0;
pass_all := 0;
bad_7_11 := 0;
bad_7_11_19_or_pass19 := 0;
reports := 0;

tag_counts := AssociativeArray();
pattern_counts := AssociativeArray();
residue_counts := AssociativeArray();

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
        f, t, A, B := M14FamilyPolynomial(R,w);
        if not GoodHyperellipticPolynomial(f) then
            continue;
        end if;
        smooth +:= 1;

        pattern_parts := [];
        all_ok := true;
        bad_or_pass := AssociativeArray(Integers());
        for p in prime_list do
            tags, good, pass5 := ClassifyAtPrime(R, w, f, p);
            tag := JoinTags(tags);
            Increment(~tag_counts, Sprintf("%o:%o", p, tag));
            Append(~pattern_parts, Sprintf("%o:%o", p, tag));
            bad_or_pass[p] := (not good) or pass5;
            if good and not pass5 then
                all_ok := false;
            end if;
            rres := ResidueString(R, p);
            wres := ResidueString(w, p);
            Increment(~residue_counts, Sprintf("%o:R=%o,w=%o,%o", p, rres, wres, tag));
        end for;
        pattern := Join(pattern_parts, " | ");
        Increment(~pattern_counts, pattern);

        if bad_or_pass[7] and bad_or_pass[11] then
            bad_7_11 +:= 1;
        end if;
        if bad_or_pass[7] and bad_or_pass[11] and bad_or_pass[19] then
            bad_7_11_19_or_pass19 +:= 1;
            if reports lt max_reports then
                reports +:= 1;
                print "SURVIVES_7_11_19_NECESSARY",
                      "R", R, "w", w, "t", t, "pattern", pattern;
            end if;
        end if;
        if all_ok then
            pass_all +:= 1;
        end if;
    end for;
end for;

print "checked", checked;
print "cover", cover;
print "smooth", smooth;
print "pass_all_branch_primes", pass_all;
print "bad_or_pass_7_11", bad_7_11;
print "bad_or_pass_7_11_19", bad_7_11_19_or_pass19;

print "TAG_COUNTS";
keys := Sort(SetToSequence(Keys(tag_counts)));
for key in keys do
    print key, tag_counts[key];
end for;

print "PATTERN_COUNTS";
pattern_keys := Sort(SetToSequence(Keys(pattern_counts)));
for key in pattern_keys do
    print key, pattern_counts[key];
end for;

print "TOP_RESIDUE_COUNTS";
res_keys := Sort(SetToSequence(Keys(residue_counts)));
for key in res_keys do
    // Print all nontrivial boundary residues, and small good-pass samples.
    if (Position(key, "good_kill5") eq 0) or (residue_counts[key] ge 10) then
        print key, residue_counts[key];
    end if;
end for;

print "DONE";
