//////////////////////////////////////////////////////////////////////
//  Classify residue-filter survivors for halving the distinguished
//  order-8 class in the odd M_1(8,2,2) family.
//
//  This uses the same necessary finite-field condition as
//  m3222_halving_residue_search.m: at a good open prime, Q must be
//  divisible by 2 in J(F_p).  Since the open allowed sets at 7 and 11
//  are empty, every rational survivor is forced onto boundary at both
//  primes.  This script records which boundary components occur.
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 30;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned prime_bound then
    prime_bound := 43;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;
if not assigned progress_interval then
    progress_interval := 100000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;
if not assigned max_print then
    max_print := 40;
elif Type(max_print) eq MonStgElt then
    max_print := StringToInteger(max_print);
end if;

Qq := Rationals();
Z := Integers();
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

function OddPolynomial(u, v)
    qtilde := -X^2
        + (u^2*v - u^2 + u*v^2 - u - v^2 - v - 2)*X
        - (u^2 + u*v + v^2 + u + v + 1);
    return ((1-u)*X + 1)*((1-v)*X + 1)*((u+v+2)*X + 1)*qtilde;
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
    yq := u*v*(u+v+1);
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

    key := Z!u + Characteristic(F)*(Z!v);
    return true, f, yq, key;
end function;

function AllowedOpenResidues(p)
    F := GF(p);
    PF<XF> := PolynomialRing(F);
    allowed := { Z | };
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

function BoundaryTag(u0, v0, p)
    F := GF(p);
    try
        u := F!u0;
        v := F!v0;
    catch e
        return "denom";
    end try;
    PF<XF> := PolynomialRing(F);
    tags := [];
    if u eq 0 then Append(~tags, "u0"); end if;
    if v eq 0 then Append(~tags, "v0"); end if;
    if u eq v then Append(~tags, "u=v"); end if;
    if u eq 1 then Append(~tags, "u1"); end if;
    if v eq 1 then Append(~tags, "v1"); end if;
    if u+v+1 eq 0 then Append(~tags, "u+v+1"); end if;
    if u+v+2 eq 0 then Append(~tags, "u+v+2"); end if;

    open, f, yq, key := FiniteOpenData(u, v, F, PF);
    if open then
        Append(~tags, "open");
    elif #tags eq 0 then
        Append(~tags, "disc_or_other");
    end if;
    return Join(tags, "&");
end function;

procedure Bump(~A, key)
    if IsDefined(A, key) then
        A[key] +:= 1;
    else
        A[key] := 1;
    end if;
end procedure;

print "M_1(8,2,2) halving residue-survivor classifier";
print "height", height, "prime_bound", prime_bound, "primes", prime_list;

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
checked := 0;
curve_ok := 0;
survivors := 0;
tag_counts := AssociativeArray();
tag7_counts := AssociativeArray();
tag11_counts := AssociativeArray();
samples := [];

for u0 in params do
    for v0 in params do
        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "progress", checked, "curve_ok", curve_ok, "survivors", survivors;
        end if;

        ok, f, yq := CurveData(u0, v0);
        if not ok then
            continue;
        end if;
        curve_ok +:= 1;

        pass, badp := ResidueSieve(u0, v0, residue_data);
        if not pass then
            continue;
        end if;
        survivors +:= 1;
        tag7 := BoundaryTag(u0, v0, 7);
        tag11 := BoundaryTag(u0, v0, 11);
        key := Sprintf("7:%o | 11:%o", tag7, tag11);
        Bump(~tag_counts, key);
        Bump(~tag7_counts, tag7);
        Bump(~tag11_counts, tag11);
        if #samples lt max_print then
            Append(~samples, <u0, v0, tag7, tag11>);
        end if;
    end for;
end for;

print "DONE";
print "checked", checked;
print "curve_ok", curve_ok;
print "survivors", survivors;
print "tag7_counts";
for k in Sort(Setseq(Keys(tag7_counts))) do
    print k, tag7_counts[k];
end for;
print "tag11_counts";
for k in Sort(Setseq(Keys(tag11_counts))) do
    print k, tag11_counts[k];
end for;
print "joint_tag_counts";
for k in Sort(Setseq(Keys(tag_counts))) do
    print k, tag_counts[k];
end for;
print "samples";
for s in samples do
    print s;
end for;

quit;
