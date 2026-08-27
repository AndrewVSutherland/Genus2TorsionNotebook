//////////////////////////////////////////////////////////////////////
//  CRT-filtered rational search on the p=0 boundary-exit chart for
//  halving Q in M_1(8,2,2).
//
//  Symmetric parameters:
//      s = u+v, p = uv.
//
//  The local boundary analysis says that any rational point coming from
//  this chart should satisfy, at good 7- and 11-adic reduction,
//
//      p == 0 mod 77,
//      s mod 77 in {10,27,38,76},
//      A-B-1 == 0 mod 77.
//
//  This script searches
//      p = 77*q,  A = B + 1 + 77*r,
//  checks that s^2-4p is a rational square, then tests the residual
//  quartic square condition and verifies exact halving in the Jacobian.
//////////////////////////////////////////////////////////////////////

if not assigned height_s then
    height_s := 40;
elif Type(height_s) eq MonStgElt then
    height_s := StringToInteger(height_s);
end if;
if not assigned height_q then
    height_q := 20;
elif Type(height_q) eq MonStgElt then
    height_q := StringToInteger(height_q);
end if;
if not assigned height_b then
    height_b := 20;
elif Type(height_b) eq MonStgElt then
    height_b := StringToInteger(height_b);
end if;
if not assigned height_r then
    height_r := 5;
elif Type(height_r) eq MonStgElt then
    height_r := StringToInteger(height_r);
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
crt_s_residues := {10, 27, 38, 76};

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

function DefinedMod(q, ell)
    return (Denominator(Qq!q) mod ell) ne 0;
end function;

function ResidueInteger(q, ell)
    F := GF(ell);
    return Integers()!(F!q);
end function;

function PassesCRT(s0, p0, A0, B0)
    for ell in [7,11] do
        if not (DefinedMod(s0, ell) and DefinedMod(p0, ell) and
                DefinedMod(A0, ell) and DefinedMod(B0, ell)) then
            return false;
        end if;
        if ResidueInteger(p0, ell) ne 0 then
            return false;
        end if;
        if ell eq 7 and ResidueInteger(s0, ell) notin {3,6} then
            return false;
        end if;
        if ell eq 11 and ResidueInteger(s0, ell) notin {5,10} then
            return false;
        end if;
        if ResidueInteger(A0 - B0 - 1, ell) ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function FamilyPolynomialSP(s0, p0)
    qtilde := -X^2 + (p0*s0 - s0^2 + 2*p0 - s0 - 2)*X
              - (s0^2 - p0 + s0 + 1);
    return ((p0-s0+1)*X^2 + (2-s0)*X + 1)*((s0+2)*X + 1)*qtilde;
end function;

function CurveDataSP(s0, p0)
    if p0 eq 0 or s0+1 eq 0 or s0+2 eq 0 or p0-s0+1 eq 0 then
        return false, _, _;
    end if;
    if s0^2 - 4*p0 eq 0 then
        return false, _, _;
    end if;
    f := FamilyPolynomialSP(s0, p0);
    yq := p0*(s0+1);
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false, _, _;
    end if;
    if yq eq 0 or yq^2 ne Evaluate(f, Qq!-1) then
        return false, _, _;
    end if;
    return true, f, yq;
end function;

function ResidualQuarticSP(s0, p0, A0, B0)
    f := FamilyPolynomialSP(s0, p0);
    L := Coefficient(f, 5);
    a := X^2 + A0*X + B0;
    return f - L*(X+1)*a^2;
end function;

function LocalResidualSquareFilter(s0, p0, A0, B0)
    for ell in local_primes do
        F := GF(ell);
        PF<XF> := PolynomialRing(F);
        try
            sF := F!s0;
            pF := F!p0;
            AF := F!A0;
            BF := F!B0;
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

function ExactHalvingSP(s0, p0)
    ok, f, yq := CurveDataSP(s0, p0);
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

svals_all := RationalParametersOfHeight(height_s);
qvals_all := RationalParametersOfHeight(height_q);
bvals := RationalParametersOfHeight(height_b);
rvals := RationalParametersOfHeight(height_r);

svals := [s0 : s0 in svals_all |
    DefinedMod(s0,7) and DefinedMod(s0,11) and
    ResidueInteger(s0,7) in {3,6} and ResidueInteger(s0,11) in {5,10}
];
qvals := [q0 : q0 in qvals_all |
    q0 ne 0 and DefinedMod(q0,7) and DefinedMod(q0,11)
];

print "p=0 chart CRT-filtered rational search";
print "height_s", height_s, "height_q", height_q,
      "height_b", height_b, "height_r", height_r;
print "svals", #svals, "qvals", #qvals, "bvals", #bvals, "rvals", #rvals;

base_checked := 0;
base_split := 0;
base_curve_ok := 0;
tuples_checked := 0;
local_pass := 0;
reject_counts := AssociativeArray(Integers());
for ell in local_primes do
    reject_counts[ell] := 0;
end for;
quartic_square := 0;
exact_checked := 0;
hits := [];

for s0 in svals do
    for q0 in qvals do
        if #hits ge max_hits then
            break s0;
        end if;
        p0 := 77*q0;
        base_checked +:= 1;

        delta := s0^2 - 4*p0;
        if not IsSquareQ(delta) then
            continue;
        end if;
        base_split +:= 1;

        ok, f, yq := CurveDataSP(s0, p0);
        if not ok then
            continue;
        end if;
        base_curve_ok +:= 1;

        for B0 in bvals do
            for r0 in rvals do
                if #hits ge max_hits then
                    break s0;
                end if;
                A0 := B0 + 1 + 77*r0;
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

                if not PassesCRT(s0, p0, A0, B0) then
                    continue;
                end if;

                local_ok, badp := LocalResidualSquareFilter(s0, p0, A0, B0);
                if not local_ok then
                    reject_counts[badp] +:= 1;
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
                d := SqrtQ(delta);
                u := (s0+d)/2;
                v := (s0-d)/2;
                if exact_ok and divisible then
                    invs := TorsionInvariants(J);
                    Append(~hits, <s0,p0,u,v,A0,B0,ellpoly,ordQ,Order(half),invs,fI>);
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
for ell in local_primes do
    print "rejected_at", ell, reject_counts[ell];
end for;
print "quartic_square", quartic_square;
print "exact_checked", exact_checked;
print "hits", #hits;
for H in hits do
    print H;
end for;

quit;
