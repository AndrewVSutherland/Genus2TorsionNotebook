//////////////////////////////////////////////////////////////////////
//  Prime-by-prime diagnostic for M_1(8,4) tangent cover plus 5.
//
//  This enumerates the same (R,w) chart as the m14_search mode in
//  elkies5_4_6_hybrid_search.m, but records which good-reduction primes
//  kill the necessary condition 5 | #J(F_p).
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 30;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

if not assigned max_reports then
    max_reports := 20;
elif Type(max_reports) eq MonStgElt then
    max_reports := StringToInteger(max_reports);
end if;

Q := Rationals();
Z := Integers();
Qx<x> := PolynomialRing(Q);
prime_list := [3,7,11,13,17,19,23,29,31,37,41,43];

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

params := RationalParametersOfHeight(height);

print "M_1(8,4) plus 5 prime diagnostic";
print "height", height, "params", #params, "primes", prime_list;

checked := 0;
cover := 0;
smooth := 0;
pass_all := 0;
reports := 0;

good_counts := AssociativeArray(Integers());
bad_counts := AssociativeArray(Integers());
pass_counts := AssociativeArray(Integers());
kill_counts := AssociativeArray(Integers());
first_kill_counts := AssociativeArray(Integers());

for p in prime_list do
    good_counts[p] := 0;
    bad_counts[p] := 0;
    pass_counts[p] := 0;
    kill_counts[p] := 0;
    first_kill_counts[p] := 0;
end for;

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

        first_fail := 0;
        for p in prime_list do
            try
                fp := ChangeRing(f, GF(p));
            catch e
                bad_counts[p] +:= 1;
                continue;
            end try;
            if not GoodHyperellipticPolynomial(fp) then
                bad_counts[p] +:= 1;
                continue;
            end if;
            good_counts[p] +:= 1;
            C := HyperellipticCurve(fp);
            if (#Jacobian(C) mod 5) eq 0 then
                pass_counts[p] +:= 1;
            else
                kill_counts[p] +:= 1;
                if first_fail eq 0 then
                    first_fail := p;
                end if;
            end if;
        end for;

        if first_fail eq 0 then
            pass_all +:= 1;
            if reports lt max_reports then
                reports +:= 1;
                print "PASS_ALL", "R", R, "w", w, "t", t;
            end if;
        else
            first_kill_counts[first_fail] +:= 1;
        end if;
    end for;
end for;

print "checked", checked;
print "cover", cover;
print "smooth", smooth;
print "pass_all", pass_all;

for p in prime_list do
    print "PRIME", p,
          "good", good_counts[p],
          "bad_or_boundary", bad_counts[p],
          "pass5", pass_counts[p],
          "kill5", kill_counts[p],
          "first_kill", first_kill_counts[p];
end for;

print "DONE";
