//////////////////////////////////////////////////////////////////////
//  Torsion-structure scan inside the M_1(8,4) tangent-cover family.
//
//  The family gives rational torsion containing [4,8].  This script
//  searches exact tangent-cover points and exact-computes J(Q)_tors,
//  looking especially for [8,8] or [4,16].
//
//  Typical run:
//      magma -b height:=50 max_tests:=2000 \
//          code/m18_m14_torsion_structure_search.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned height then
    height := 40;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

if not assigned max_tests then
    max_tests := 1000;
elif Type(max_tests) eq MonStgElt then
    max_tests := StringToInteger(max_tests);
end if;

if not assigned progress_interval then
    progress_interval := 1000000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

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

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function FamilyPolynomial(R, w)
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

function TangentCandidates(f, R, w)
    h := ExactQuotient(f, x);
    c1 := Coefficient(h, 1);
    c2 := Coefficient(h, 2);
    c3 := Coefficient(h, 3);
    c4 := Coefficient(h, 4);

    out := [];
    PR<U> := PolynomialRing(Q);

    for V in [R^2*w, -R^2*w] do
        F := 4*(c3 - 2*c4*U)*(c1 - 2*c4*U*V)
             - (c2 - c4*(U^2 + 2*V))^2;
        for rt in Roots(F) do
            U0 := rt[1];
            M2 := c3 - 2*c4*U0;
            N2 := c1 - 2*c4*U0*V;
            if IsSquareQ(M2) and IsSquareQ(N2) then
                Append(~out, <U0, V, M2, N2>);
            end if;
        end for;
    end for;
    return out;
end function;

function TorsionOrder(invs)
    if #invs eq 0 then
        return 1;
    end if;
    return &*invs;
end function;

function Exponent(invs)
    if #invs eq 0 then
        return 1;
    end if;
    return invs[#invs];
end function;

function TwoAdicExponents(invs)
    vals := Sort([Valuation(n, 2) : n in invs]);
    return Reverse(vals);
end function;

function IsTarget(invs)
    vals := TwoAdicExponents(invs);
    has88 := #vals ge 2 and vals[1] ge 3 and vals[2] ge 3;
    has416 := #vals ge 2 and vals[1] ge 4 and vals[2] ge 2;
    return has88 or has416 or Exponent(invs) ge 16;
end function;

params := RationalParametersOfHeight(height);
checked := 0;
cover := 0;
smooth := 0;
tangent_points := 0;
verified := 0;
exact_tests := 0;
targets := [];
torsion_counts := AssociativeArray();
max_order := 0;
max_exponent := 0;

print "M_1(8,4) torsion-structure scan";
print "height", height, "parameters", #params,
      "max_tests", max_tests, "progress_interval", progress_interval;

for R in params do
    for w in params do
        if R eq 0 or w in {Q!-1, Q!0, Q!1} then
            continue;
        end if;
        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "progress", checked, "cover", cover, "smooth", smooth,
                  "tangent_points", tangent_points, "exact", exact_tests,
                  "targets", #targets, "max_order", max_order,
                  "max_exponent", max_exponent;
        end if;

        if not M14CoverPossible(R,w) then
            continue;
        end if;
        cover +:= 1;

        f, t, A, B := FamilyPolynomial(R,w);
        if Degree(f) ne 5 or Discriminant(f) eq 0 then
            continue;
        end if;
        smooth +:= 1;

        candidates := TangentCandidates(f, R, w);
        if #candidates eq 0 then
            continue;
        end if;
        tangent_points +:= 1;

        fI, L := IntegralModelPolynomial(f);
        C := HyperellipticCurve(fI);
        J := Jacobian(C);
        D := J![x, Q!0];
        divisible, half := IsDivisibleBy(D, 2);
        if not divisible then
            print "WARNING tangent equations did not verify", R, w, t, candidates;
            continue;
        end if;
        verified +:= 1;

        exact_tests +:= 1;
        G, phi := TorsionSubgroup(J);
        invs := Invariants(G);
        key := Sprint(invs);
        if IsDefined(torsion_counts, key) then
            torsion_counts[key] +:= 1;
        else
            torsion_counts[key] := 1;
        end if;
        ord := TorsionOrder(invs);
        exp := Exponent(invs);
        if ord gt max_order then
            max_order := ord;
        end if;
        if exp gt max_exponent then
            max_exponent := exp;
        end if;

        if invs ne [4,8] then
            print "NONSTANDARD", "R", R, "w", w, "t", t,
                  "half_order", Order(half), "torsion", invs,
                  "order", ord, "exponent", exp;
            print "  f =", fI;
        end if;

        if IsTarget(invs) then
            Append(~targets, <R,w,t,invs,fI>);
            print "TARGET", "R", R, "w", w, "t", t,
                  "torsion", invs, "order", ord, "exponent", exp;
            print "  f =", fI;
        end if;

        if exact_tests ge max_tests then
            break R;
        end if;
    end for;
end for;

print "DONE height", height;
print "checked", checked, "cover", cover, "smooth", smooth,
      "tangent_points", tangent_points, "verified", verified,
      "exact_tests", exact_tests, "targets", #targets,
      "max_order", max_order, "max_exponent", max_exponent;
print "TORSION_COUNTS";
for key in Sort([k : k in Keys(torsion_counts)]) do
    print " ", key, torsion_counts[key];
end for;

quit;
