//////////////////////////////////////////////////////////////////////
//  Rational search on the symmetric halving surface.
//
//  Parameters:
//      s = u+v, p = uv,
//      a(X) = X^2 + A*X + B.
//
//  We require:
//      1. s^2 - 4p is a rational square, so u,v are rational;
//      2. f_{s,p}(X) - L*(X+1)*a(X)^2 is a square over QQ[X];
//      3. Magma's exact Jacobian IsDivisibleBy verifies Q is halved.
//
//  Typical run from torsion_jac:
//      magma -b height:=5 code/m3222_halving_sp_search.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 5;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;
if not assigned progress_interval then
    progress_interval := 100000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Qq := Rationals();
P<X> := PolynomialRing(Qq);
local_primes := [7,11,13,17,19,23,29,31,37,41,43];

function RationalParametersOfHeight(Bnd)
    vals := [];
    seen := {};
    for den in [1..Bnd] do
        for num in [-Bnd..Bnd] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            r := Qq!num/den;
            key := Sprint(r);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, r);
            end if;
        end for;
    end for;
    return vals;
end function;

function IsSquareQ(a)
    a := Qq!a;
    if a lt 0 then
        return false;
    end if;
    return IsSquare(Numerator(a)) and IsSquare(Denominator(a));
end function;

function SqrtQ(a)
    return Qq!Sqrt(Numerator(a)) / Qq!Sqrt(Denominator(a));
end function;

function IntegralModelPolynomial(f)
    Lden := 1;
    for i in [0..Degree(f)] do
        Lden := LCM(Lden, Denominator(Coefficient(f, i)));
    end for;
    return P!(Lden^2*f), Lden;
end function;

function FamilyPolynomialSP(s, pp)
    qtilde := -X^2 + (pp*s - s^2 + 2*pp - s - 2)*X
              - (s^2 - pp + s + 1);
    return ((pp-s+1)*X^2 + (2-s)*X + 1)*((s+2)*X + 1)*qtilde;
end function;

function BoundaryDegenerateSP(s, pp)
    if pp eq 0 or s+1 eq 0 or s+2 eq 0 or pp-s+1 eq 0 then
        return true;
    end if;
    if s^2 - 4*pp eq 0 then
        return true;
    end if;
    return false;
end function;

function CurveDataSP(s, pp)
    if BoundaryDegenerateSP(s, pp) then
        return false, _, _;
    end if;
    f := FamilyPolynomialSP(s, pp);
    yq := pp*(s+1);
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false, _, _;
    end if;
    if yq eq 0 or yq^2 ne Evaluate(f, Qq!-1) then
        return false, _, _;
    end if;
    return true, f, yq;
end function;

function ResidualQuarticSP(s, pp, A, B)
    f := FamilyPolynomialSP(s, pp);
    L := Coefficient(f, 5);
    a := X^2 + A*X + B;
    return f - L*(X+1)*a^2;
end function;

function LocalResidualSquareFilter(s, pp, A, B)
    for ell in local_primes do
        F := GF(ell);
        PF<XF> := PolynomialRing(F);
        try
            sF := F!s;
            pF := F!pp;
            AF := F!A;
            BF := F!B;
        catch e
            continue;
        end try;

        qtilde := -XF^2 + (pF*sF - sF^2 + 2*pF - sF - 2)*XF
                  - (sF^2 - pF + sF + 1);
        fF := ((pF-sF+1)*XF^2 + (2-sF)*XF + 1)*((sF+2)*XF + 1)*qtilde;
        LF := Coefficient(fF, 5);
        aF := XF^2 + AF*XF + BF;
        resF := fF - LF*(XF+1)*aF^2;
        if not IsSquare(resF) then
            return false, ell;
        end if;
    end for;
    return true, 0;
end function;

function ExactHalvingSP(s, pp)
    ok, f, yq := CurveDataSP(s, pp);
    if not ok then
        return false, _, _, _, _, _;
    end if;
    fI, Lden := IntegralModelPolynomial(f);
    C := HyperellipticCurve(fI);
    J := Jacobian(C);
    DQ := J![X + 1, Qq!(Lden*yq)];
    ordQ := Order(DQ);
    half := J!0;
    divisible := IsDivisibleBy(DQ, 2);
    if divisible then
        _, half := IsDivisibleBy(DQ, 2);
    end if;
    return true, ordQ, divisible, half, J, fI;
end function;

function TorsionInvariants(J)
    G, phi := TorsionSubgroup(J);
    return Invariants(G);
end function;

params := RationalParametersOfHeight(height);
print "Symmetric halving-surface rational search";
print "height", height, "params", #params, "tuples", #params^4;

base_checked := 0;
base_split := 0;
base_curve_ok := 0;
tuples_checked := 0;
local_pass := 0;
quartic_square := 0;
exact_checked := 0;
hits := [];

for s0 in params do
    for p0 in params do
        if #hits ge max_hits then
            break s0;
        end if;
        base_checked +:= 1;
        delta := s0^2 - 4*p0;
        if delta eq 0 or not IsSquareQ(delta) then
            continue;
        end if;
        base_split +:= 1;
        ok, f, yq := CurveDataSP(s0, p0);
        if not ok then
            continue;
        end if;
        base_curve_ok +:= 1;

        for A0 in params do
            for B0 in params do
                if #hits ge max_hits then
                    break s0;
                end if;
                tuples_checked +:= 1;
                if progress_interval gt 0 and tuples_checked mod progress_interval eq 0 then
                    print "progress tuples", tuples_checked,
                          "base_checked", base_checked,
                          "base_curve_ok", base_curve_ok,
                          "local_pass", local_pass,
                          "quartic_square", quartic_square,
                          "exact", exact_checked,
                          "hits", #hits;
                end if;

                local_ok, badp := LocalResidualSquareFilter(s0, p0, A0, B0);
                if not local_ok then
                    continue;
                end if;
                local_pass +:= 1;

                res := ResidualQuarticSP(s0, p0, A0, B0);
                square_ok, ellpoly := IsSquare(res);
                if not square_ok then
                    continue;
                end if;
                quartic_square +:= 1;

                exact_checked +:= 1;
                exact_ok, ordQ, divisible, half, J, fI := ExactHalvingSP(s0, p0);
                if not exact_ok then
                    continue;
                end if;
                d := SqrtQ(delta);
                u := (s0+d)/2;
                v := (s0-d)/2;

                if divisible then
                    invs := TorsionInvariants(J);
                    Append(~hits, <s0, p0, u, v, A0, B0, ellpoly, ordQ, Order(half), invs, fI>);
                    print "HIT", "s", s0, "p", p0, "u", u, "v", v,
                          "A", A0, "B", B0,
                          "order_Q", ordQ, "half_order", Order(half),
                          "torsion", invs;
                    print "  ell", ellpoly;
                    print "  fI", fI;
                else
                    print "SQUARE SURVIVOR NOT DIVISIBLE",
                          "s", s0, "p", p0, "u", u, "v", v,
                          "A", A0, "B", B0, "ordQ", ordQ;
                end if;
            end for;
        end for;
    end for;
end for;

print "DONE";
print "base_checked", base_checked;
print "base_split", base_split;
print "base_curve_ok", base_curve_ok;
print "tuples_checked", tuples_checked;
print "local_pass", local_pass;
print "quartic_square", quartic_square;
print "exact_checked", exact_checked;
print "hits", #hits;
for H in hits do
    print H;
end for;

quit;
