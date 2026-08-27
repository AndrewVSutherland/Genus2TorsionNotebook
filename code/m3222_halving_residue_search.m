//////////////////////////////////////////////////////////////////////
//  Fast rational search for rational halves of the distinguished
//  order-8 class in the odd M_1(8,2,2) family.
//
//  This is the same target as m3222_halving_search.m, but it first
//  precomputes, for several primes, the open finite-field residue
//  classes where Q is divisible by 2.  A rational pair (u,v) is rejected
//  immediately if it has good open reduction modulo p outside that set.
//
//  Typical run from torsion_jac:
//      magma -b height:=50 code/m3222_halving_residue_search.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 30;
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
if not assigned prime_bound then
    prime_bound := 29;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;

Qq := Rationals();
P<X> := PolynomialRing(Qq);
prime_list := [p : p in [7,11,13,17,19,23,29,31,37,41,43] | p le prime_bound];

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

function IsDivisibleBy2Finite(J, D)
    G, phi := AbelianGroup(J);
    a := D @@ phi;
    coords := Eltseq(a);
    invs := Invariants(G);
    for i in [1..#coords] do
        if GCD(2, invs[i]) ne 1 and coords[i] mod 2 ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function FiniteOpenData(u, v, F, PF)
    XF := PF.1;
    if u eq 0 or v eq 0 or u eq v then
        return false, _, _, _;
    end if;
    if u eq 1 or v eq 1 or u+v+1 eq 0 or u+v+2 eq 0 then
        return false, _, _, _;
    end if;

    qtilde := -XF^2
        + (u^2*v - u^2 + u*v^2 - u - v^2 - v - 2)*XF
        - (u^2 + u*v + v^2 + u + v + 1);
    f := ((1-u)*XF + 1)*((1-v)*XF + 1)*((u+v+2)*XF + 1)*qtilde;
    yq := u*v*(u+v+1);

    if Degree(f) ne 5 or Discriminant(f) eq 0 or yq eq 0 then
        return false, _, _, _;
    end if;
    if yq^2 ne Evaluate(f, F!-1) then
        return false, _, _, _;
    end if;

    key := Integers()!u + Characteristic(F)*(Integers()!v);
    return true, f, yq, key;
end function;

function AllowedOpenResidues(p)
    F := GF(p);
    PF<XF> := PolynomialRing(F);
    allowed := { Integers() | };
    nonsingular := 0;
    order8 := 0;

    for u in F do
        for v in F do
            open, f, yq, key := FiniteOpenData(u, v, F, PF);
            if not open then
                continue;
            end if;
            nonsingular +:= 1;
            C := HyperellipticCurve(f);
            J := Jacobian(C);
            DQ := J![XF + F!1, yq];
            if Order(DQ) eq 8 then
                order8 +:= 1;
            end if;
            if IsDivisibleBy2Finite(J, DQ) then
                Include(~allowed, key);
            end if;
        end for;
    end for;

    return allowed, nonsingular, order8;
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

    open, f, yq, key := FiniteOpenData(u, v, F, PF);
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

function ExactHalvingData(u, v)
    ok, f, yq := CurveData(u, v);
    if not ok then
        return false, _, _, _, _, _, _, _, _;
    end if;

    fI, L := IntegralModelPolynomial(f);
    if Degree(fI) ne 5 or Discriminant(fI) eq 0 then
        return false, _, _, _, _, _, _, _, _;
    end if;

    C := HyperellipticCurve(fI);
    J := Jacobian(C);
    DQ := J![X + 1, Qq!(L*yq)];
    ordQ := Order(DQ);
    half := J!0;
    divisible := IsDivisibleBy(DQ, 2);
    if divisible then
        _, half := IsDivisibleBy(DQ, 2);
        ordHalf := Order(half);
    else
        ordHalf := 0;
    end if;

    return true, fI, L, ordQ, divisible, half, ordHalf, J, yq;
end function;

function TorsionInvariants(J)
    G, phi := TorsionSubgroup(J);
    return Invariants(G);
end function;

print "M_1(8,2,2) fast Q-halving residue search";
print "height", height, "prime_bound", prime_bound, "max_hits", max_hits;
print "primes", prime_list;

residue_data := [];
for p in prime_list do
    allowed, nonsingular, order8 := AllowedOpenResidues(p);
    Append(~residue_data, <p, allowed>);
    print "local p", p,
          "open", nonsingular,
          "Q_order8", order8,
          "allowed_open_halving_residues", #allowed;
end for;

params := RationalParametersOfHeight(height);
print "params", #params, "pairs", #params^2;

checked := 0;
curve_ok := 0;
residue_pass := 0;
exact_checked := 0;
hits := [];
reject_counts := AssociativeArray(Integers());
for p in prime_list do
    reject_counts[p] := 0;
end for;

for u0 in params do
    for v0 in params do
        if #hits ge max_hits then
            break u0;
        end if;
        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "progress checked", checked,
                  "curve_ok", curve_ok,
                  "residue_pass", residue_pass,
                  "exact", exact_checked,
                  "hits", #hits;
        end if;

        ok, f, yq := CurveData(u0, v0);
        if not ok then
            continue;
        end if;
        curve_ok +:= 1;

        pass, badp := ResidueSieve(u0, v0, residue_data);
        if not pass then
            reject_counts[badp] +:= 1;
            continue;
        end if;
        residue_pass +:= 1;

        exact_checked +:= 1;
        exact_ok, fI, L, ordQ, divisible, half, ordHalf, J, yq := ExactHalvingData(u0, v0);
        if not exact_ok then
            continue;
        end if;

        if divisible then
            invs := TorsionInvariants(J);
            Append(~hits, <u0, v0, ordQ, ordHalf, invs, fI, L>);
            print "HIT", "u", u0, "v", v0,
                  "order_Q", ordQ, "half_order", ordHalf,
                  "torsion", invs, "L", L;
            print "  fI", fI;
        end if;

        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "progress checked", checked,
                  "curve_ok", curve_ok,
                  "residue_pass", residue_pass,
                  "exact", exact_checked,
                  "hits", #hits;
        end if;
    end for;
end for;

print "DONE";
print "checked", checked;
print "curve_ok", curve_ok;
print "residue_pass", residue_pass;
print "exact_checked", exact_checked;
print "hits", #hits;
for p in prime_list do
    print "rejected_at", p, reject_counts[p];
end for;
for H in hits do
    print H;
end for;

quit;
