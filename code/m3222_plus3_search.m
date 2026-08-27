//////////////////////////////////////////////////////////////////////
//  M_1(8,2,2) / [2,2,8] family plus possible rational 3-torsion.
//
//  Odd model from the existing m3222 scripts:
//      C_{u,v}: Y^2 = f_{u,v}(X),
//      Q = (-1, u*v*(u+v+1))
//  with Q of order 8 on the nonsingular open family.  The family also
//  has the visible rational 2-torsion from the three linear factors.
//
//  A rational 3-torsion point implies 3 | #J(F_p) at every good prime
//  p != 3.  This script precomputes the open residue classes satisfying
//  that necessary condition, scans rational (u,v), and exact-computes
//  torsion only for survivors.
//
//  Typical run:
//      magma -b height:=30 prime_bound:=101 code/m3222_plus3_search.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 20;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned prime_bound then
    prime_bound := 101;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;
if not assigned max_exact then
    max_exact := 200;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned progress_interval then
    progress_interval := 100000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;
if not assigned max_print then
    max_print := 20;
elif Type(max_print) eq MonStgElt then
    max_print := StringToInteger(max_print);
end if;

Qq := Rationals();
Z := Integers();
P<X> := PolynomialRing(Qq);
prime_list := [ p : p in PrimesUpTo(prime_bound) | p notin {2,3} ];

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
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

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function OddPolynomial(u, v)
    qtilde := -X^2
        + (u^2*v - u^2 + u*v^2 - u - v^2 - v - 2)*X
        - (u^2 + u*v + v^2 + u + v + 1);
    return ((1-u)*X + 1)*((1-v)*X + 1)*((u+v+2)*X + 1)*qtilde;
end function;

function QY(u, v)
    return u*v*(u+v+1);
end function;

function HasBasicDegeneracy(u, v)
    if u eq 0 or v eq 0 or u eq v then
        return true;
    end if;
    if u eq 1 or v eq 1 or u+v+1 eq 0 or u+v+2 eq 0 then
        return true;
    end if;
    return false;
end function;

function CurveData(u, v)
    if HasBasicDegeneracy(u, v) then
        return false, _, _;
    end if;
    f := OddPolynomial(u, v);
    yq := QY(u, v);
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false, _, _;
    end if;
    if yq eq 0 or yq^2 ne Evaluate(f, Qq!-1) then
        return false, _, _;
    end if;
    return true, f, yq;
end function;

function FiniteOpenData(u, v, F, PF)
    XF := PF.1;
    if u eq 0 or v eq 0 or u eq v then
        return false, _, _;
    end if;
    if u eq 1 or v eq 1 or u+v+1 eq 0 or u+v+2 eq 0 then
        return false, _, _;
    end if;

    qtilde := -XF^2
        + (u^2*v - u^2 + u*v^2 - u - v^2 - v - 2)*XF
        - (u^2 + u*v + v^2 + u + v + 1);
    f := ((1-u)*XF + 1)*((1-v)*XF + 1)*((u+v+2)*XF + 1)*qtilde;
    yq := u*v*(u+v+1);
    if Degree(f) ne 5 or Discriminant(f) eq 0 or yq eq 0 then
        return false, _, _;
    end if;
    if yq^2 ne Evaluate(f, F!-1) then
        return false, _, _;
    end if;
    key := Z!u + Characteristic(F)*(Z!v);
    return true, f, key;
end function;

function Allowed3Residues(p)
    F := GF(p);
    PF<XF> := PolynomialRing(F);
    allowed := { Z | };
    nonsingular := 0;
    for u in F do
        for v in F do
            open, f, key := FiniteOpenData(u, v, F, PF);
            if not open then
                continue;
            end if;
            nonsingular +:= 1;
            C := HyperellipticCurve(f);
            N := Z!Evaluate(LPolynomial(C), 1);
            if N mod 3 eq 0 then
                Include(~allowed, key);
            end if;
        end for;
    end for;
    return allowed, nonsingular;
end function;

function RationalOpenKeyAtPrime(u0, v0, p)
    F := GF(p);
    PF<XF> := PolynomialRing(F);
    try
        u := F!u0;
        v := F!v0;
    catch e
        return false, 0;
    end try;
    open, f, key := FiniteOpenData(u, v, F, PF);
    if not open then
        return false, 0;
    end if;
    return true, key;
end function;

function ResidueSieve(u0, v0, residue_data)
    for data in residue_data do
        p := data[1];
        allowed := data[2];
        open, key := RationalOpenKeyAtPrime(u0, v0, p);
        if open and key notin allowed then
            return false, p;
        end if;
    end for;
    return true, 0;
end function;

function ExactTorsionData(u, v)
    ok, f, yq := CurveData(u, v);
    if not ok then
        return false, [], 0, P!0, 0;
    end if;
    fI, L := IntegralModelPolynomial(f);
    if Degree(fI) ne 5 or Discriminant(fI) eq 0 then
        return false, [], 0, fI, L;
    end if;
    C := HyperellipticCurve(fI);
    J := Jacobian(C);
    G, phi := TorsionSubgroup(J);
    return true, Invariants(G), &*Invariants(G), fI, L;
end function;

function BoundGCD(f, prime_list)
    g := 0;
    used := 0;
    for p in prime_list do
        F := GF(p);
        PF<XF> := PolynomialRing(F);
        try
            fp := PF!ChangeRing(f, F);
        catch e
            continue;
        end try;
        if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
            continue;
        end if;
        C := HyperellipticCurve(fp);
        N := Z!Evaluate(LPolynomial(C), 1);
        if g eq 0 then
            g := N;
        else
            g := GCD(g, N);
        end if;
        used +:= 1;
    end for;
    return g, used;
end function;

residue_data := [];
print "M_1(8,2,2) / [2,2,8] plus 3-torsion search";
print "height", height, "prime_bound", prime_bound, "max_exact", max_exact;
print "Precomputing finite 3-residue filters";
for p in prime_list do
    allowed, nonsingular := Allowed3Residues(p);
    Append(~residue_data, <p, allowed>);
    print " p", p, "nonsingular", nonsingular, "allowed3", #allowed;
end for;

params := RationalParametersOfHeight(height);
checked := 0;
curve_ok := 0;
survivors := 0;
exact := 0;
hits := 0;
rejected := AssociativeArray();
bound_counts := AssociativeArray();
printed := 0;

for u in params do
    for v in params do
        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "progress", checked, "curve_ok", curve_ok, "survivors", survivors,
                  "exact", exact, "hits", hits;
        end if;

        ok, f, yq := CurveData(u, v);
        if not ok then
            continue;
        end if;
        curve_ok +:= 1;

        pass, badp := ResidueSieve(u, v, residue_data);
        if not pass then
            if IsDefined(rejected, badp) then
                rejected[badp] +:= 1;
            else
                rejected[badp] := 1;
            end if;
            continue;
        end if;
        survivors +:= 1;

        g, used := BoundGCD(f, prime_list);
        key := Sprint(g);
        if IsDefined(bound_counts, key) then
            bound_counts[key] +:= 1;
        else
            bound_counts[key] := 1;
        end if;

        if exact ge max_exact then
            continue;
        end if;
        exact +:= 1;
        exact_ok, invs, tor_order, fI, L := ExactTorsionData(u, v);
        if exact_ok and tor_order mod 3 eq 0 then
            hits +:= 1;
            print "HIT", "u", u, "v", v, "torsion", invs, "order", tor_order,
                  "bound", g, "used_primes", used;
            print "  f =", fI;
        elif printed lt max_print then
            print "SURVIVOR", "u", u, "v", v, "torsion", invs,
                  "order", tor_order, "bound", g, "used_primes", used;
            printed +:= 1;
        end if;
    end for;
end for;

print "DONE height", height;
print "checked", checked, "curve_ok", curve_ok, "survivors", survivors,
      "exact", exact, "hits", hits;
print "Rejected first at";
for p in Sort([ k : k in Keys(rejected) ]) do
    print " ", p, rejected[p];
end for;
print "Surviving gcd bounds";
for key in Sort([ k : k in Keys(bound_counts) ]) do
    print " ", key, bound_counts[key];
end for;

quit;
