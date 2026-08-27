//////////////////////////////////////////////////////////////////////
//  Search in the odd-degree M_1(8,2,2) model for rational halves of
//  the distinguished order-8 class.
//
//  In the notation of paper/NotesAndTodo.tex, Section (16,2,2),
//
//      C_{u,v}: Y^2 = f_{u,v}(X),
//      Q = (-1, u*v*(u+v+1)).
//
//  The divisor [Q - infinity] should have order 8 on the nonsingular
//  open of the family.  If it is divisible by 2 over QQ, then its half
//  has order 16, giving the desired [2,2,2,16]-type behavior together
//  with the rational Weierstrass points already present in the family.
//
//  Typical runs from torsion_jac:
//      magma -b mode:=check_one u:=2 v:=3 code/m3222_halving_search.m
//      magma -b mode:=halve_q height:=12 code/m3222_halving_search.m
//      magma -b mode:=cantor32 height:=8 code/m3222_halving_search.m
//////////////////////////////////////////////////////////////////////

if not assigned mode then
    mode := "halve_q";
end if;
if not assigned height then
    height := 10;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;
if not assigned progress_interval then
    progress_interval := 500;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;
if not assigned use_local_sieve then
    use_local_sieve := true;
elif Type(use_local_sieve) eq MonStgElt then
    use_local_sieve := use_local_sieve in {"true", "True", "1", "yes"};
end if;

Qq := Rationals();
P<X> := PolynomialRing(Qq);
local_primes := [7,11,13,17,19,23,29,31,37,41,43,3,5];

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

function IsSquareQ(a)
    a := Qq!a;
    if a lt 0 then
        return false;
    end if;
    return IsSquare(Numerator(a)) and IsSquare(Denominator(a));
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
        return false, _, _, _;
    end if;

    f := OddPolynomial(u, v);
    yq := QY(u, v);

    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false, _, _, _;
    end if;
    if yq eq 0 or yq^2 ne Evaluate(f, Qq!-1) then
        return false, _, _, _;
    end if;

    return true, f, yq, Discriminant(f);
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

function LocalHalvingSieve(f, yq)
    good_primes := 0;
    for p in local_primes do
        F := GF(p);
        PF<XF> := PolynomialRing(F);

        try
            fp := PF!ChangeRing(f, F);
            yp := F!yq;
        catch e
            continue;
        end try;

        if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
            continue;
        end if;
        if yp^2 ne Evaluate(fp, F!-1) then
            continue;
        end if;

        C := HyperellipticCurve(fp);
        J := Jacobian(C);
        D := J![XF + F!1, yp];
        good_primes +:= 1;
        if not IsDivisibleBy2Finite(J, D) then
            return false, p, good_primes;
        end if;
    end for;
    return true, 0, good_primes;
end function;

function ExactHalvingData(u, v)
    ok, f, yq, disc := CurveData(u, v);
    if not ok then
        return false, _, _, _, _, _, _, _;
    end if;

    fI, L := IntegralModelPolynomial(f);
    if Degree(fI) ne 5 or Discriminant(fI) eq 0 then
        return false, _, _, _, _, _, _, _;
    end if;
    yI := Qq!(L*yq);

    C := HyperellipticCurve(fI);
    J := Jacobian(C);
    DQ := J![X + 1, yI];

    ordQ := Order(DQ);
    half := J!0;
    divisible := IsDivisibleBy(DQ, 2);
    if divisible then
        _, half := IsDivisibleBy(DQ, 2);
        ordHalf := Order(half);
    else
        ordHalf := 0;
    end if;

    return true, f, fI, L, DQ, ordQ, divisible, half, ordHalf, J;
end function;

function TorsionInvariants(J)
    G, phi := TorsionSubgroup(J);
    return Invariants(G);
end function;

function DerivativeEquations(f, t)
    f0 := Evaluate(f, t);
    f1 := Evaluate(Derivative(f, 1), t);
    f2 := Evaluate(Derivative(f, 2), t);
    f3 := Evaluate(Derivative(f, 3), t);
    f4 := Evaluate(Derivative(f, 4), t);
    a5 := Coefficient(f, 5);

    E2 := 4*f0*f3 - 6*f0*f1*f2 + 3*f1^3;
    E3 := 8*f0^3*f4
        - 3*(2*f0*f2 - f1^2)^2
        - 192*f0^3*a5*(t+1);

    return E2, E3, f0;
end function;

if mode eq "check_one" then
    if not assigned u or not assigned v then
        print "For mode:=check_one, pass u:=... v:=...";
        quit;
    end if;
    if Type(u) eq MonStgElt then
        u := StringToRational(u);
    else
        u := Qq!u;
    end if;
    if Type(v) eq MonStgElt then
        v := StringToRational(v);
    else
        v := Qq!v;
    end if;

    ok, f, yq, disc := CurveData(u, v);
    print "check_one u", u, "v", v, "curve_ok", ok;
    if ok then
        print "f", f;
        print "Q", <Qq!-1, yq>;
        if use_local_sieve then
            local_ok, badp, goodp := LocalHalvingSieve(f, yq);
            print "local_halving_possible", local_ok, "bad_prime", badp, "good_primes", goodp;
        end if;
        exact_ok, f0, fI, L, DQ, ordQ, divisible, half, ordHalf, J := ExactHalvingData(u, v);
        print "exact_ok", exact_ok;
        print "L", L;
        print "order_Q", ordQ;
        print "divisible_by_2", divisible;
        print "half_order", ordHalf;
        if divisible then
            print "torsion", TorsionInvariants(J);
        end if;
    end if;
    quit;
end if;

if mode eq "halve_q" then
    params := RationalParametersOfHeight(height);
    print "M_1(8,2,2) Q-halving search";
    print "height", height, "params", #params, "pairs", #params^2;
    print "use_local_sieve", use_local_sieve, "max_hits", max_hits;

    checked := 0;
    nonsingular := 0;
    local_pass := 0;
    exact_checked := 0;
    order8 := 0;
    divisible_count := 0;
    hits := [];

    for u0 in params do
        for v0 in params do
            if #hits ge max_hits then
                break u0;
            end if;
            checked +:= 1;

            ok, f, yq, disc := CurveData(u0, v0);
            if not ok then
                continue;
            end if;
            nonsingular +:= 1;

            if use_local_sieve then
                local_ok, badp, goodp := LocalHalvingSieve(f, yq);
                if not local_ok then
                    continue;
                end if;
            else
                local_ok := true;
                badp := 0;
                goodp := 0;
            end if;
            local_pass +:= 1;

            exact_checked +:= 1;
            exact_ok, f0, fI, L, DQ, ordQ, divisible, half, ordHalf, J := ExactHalvingData(u0, v0);
            if not exact_ok then
                continue;
            end if;
            if ordQ eq 8 then
                order8 +:= 1;
            else
                print "WARNING Q order not 8", "u", u0, "v", v0, "order", ordQ;
            end if;

            if divisible then
                divisible_count +:= 1;
                invs := TorsionInvariants(J);
                Append(~hits, <u0, v0, ordQ, ordHalf, invs, fI, L>);
                print "HIT", "u", u0, "v", v0,
                      "order_Q", ordQ, "half_order", ordHalf,
                      "torsion", invs, "L", L;
                print "  fI", fI;
            end if;

            if progress_interval gt 0 and nonsingular mod progress_interval eq 0 then
                print "progress checked", checked,
                      "nonsingular", nonsingular,
                      "local_pass", local_pass,
                      "exact", exact_checked,
                      "order8", order8,
                      "hits", #hits;
            end if;
        end for;
    end for;

    print "DONE";
    print "checked", checked;
    print "nonsingular", nonsingular;
    print "local_pass", local_pass;
    print "exact_checked", exact_checked;
    print "order8", order8;
    print "divisible", divisible_count;
    print "hits", #hits;
    for H in hits do
        print H;
    end for;
    quit;
end if;

if mode eq "cantor32" then
    params := RationalParametersOfHeight(height);
    print "M_1(8,2,2) claimed Cantor-32 equation scan";
    print "height", height, "params", #params, "triples", #params^3;

    checked := 0;
    e23 := 0;
    square_points := 0;
    verified32 := 0;
    hits := [];

    for u0 in params do
        for v0 in params do
            ok, f, yq, disc := CurveData(u0, v0);
            if not ok then
                continue;
            end if;

            for t0 in params do
                if #hits ge max_hits then
                    break u0;
                end if;
                checked +:= 1;

                E2, E3, f0 := DerivativeEquations(f, t0);
                if E2 ne 0 or E3 ne 0 then
                    continue;
                end if;
                e23 +:= 1;
                if not IsSquareQ(f0) then
                    continue;
                end if;
                square_points +:= 1;
                s0 := Sqrt(f0);

                fI, L := IntegralModelPolynomial(f);
                if Degree(fI) ne 5 or Discriminant(fI) eq 0 then
                    continue;
                end if;

                C := HyperellipticCurve(fI);
                J := Jacobian(C);
                D := J![X - t0, Qq!(L*s0)];
                DQ := J![X + 1, Qq!(L*yq)];
                ordD := Order(D);
                four_ok := (4*D eq DQ) or (4*D eq -DQ);
                if four_ok then
                    invs := TorsionInvariants(J);
                    Append(~hits, <u0, v0, t0, s0, ordD, invs, fI, L>);
                    if ordD eq 32 then
                        verified32 +:= 1;
                    end if;
                    print "CANTOR HIT", "u", u0, "v", v0,
                          "t", t0, "s", s0,
                          "order_D", ordD, "four_to_Q", four_ok,
                          "torsion", invs, "L", L;
                    print "  fI", fI;
                else
                    print "WARNING derivative equations did not verify",
                          "u", u0, "v", v0, "t", t0, "s", s0,
                          "order_D", ordD;
                end if;
            end for;
        end for;
    end for;

    print "DONE";
    print "checked", checked;
    print "E2E3_zero", e23;
    print "square_points", square_points;
    print "verified32", verified32;
    print "hits", #hits;
    for H in hits do
        print H;
    end for;
    quit;
end if;

print "Unknown mode", mode;
quit;
