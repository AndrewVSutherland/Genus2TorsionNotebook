//////////////////////////////////////////////////////////////////////
//  Targeted q^2=1 diagonal check inside the D-square + F0-square
//  surface for M(2,2,4,8).
//
//  On this diagonal, sigma = rho/tau.  Then
//
//      F0 = rho^2,
//      C  = rho^2/d^2        on the D-square surface,
//      F2 = rho^2*F1.
//
//  Thus the apparent F1/F2/C successes in the pair-fiber diagnostic are
//  dependent; the only remaining independent full-cover tests are F1
//  and F4 (then F3 follows from the product identity).
//////////////////////////////////////////////////////////////////////

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

if not assigned dn_height then
    dn_height := 30;
elif Type(dn_height) eq MonStgElt then
    dn_height := StringToInteger(dn_height);
end if;

if not assigned max_reports then
    max_reports := 30;
elif Type(max_reports) eq MonStgElt then
    max_reports := StringToInteger(max_reports);
end if;

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

function IsNonzeroSquareQ(q)
    if q eq 0 then
        return false, Q!0;
    end if;
    num := Z!Numerator(q);
    den := Z!Denominator(q);
    if num lt 0 then
        return false, Q!0;
    end if;
    okN, rtN := IsSquare(num);
    if not okN then return false, Q!0; end if;
    okD, rtD := IsSquare(den);
    if not okD then return false, Q!0; end if;
    return true, Q!rtN/Q!rtD;
end function;

function DSurfaceParam(d, n)
    den := n^2 + d^2 - 1;
    if den eq 0 then
        return false, Q!0, Q!0;
    end if;
    tau := ((1 - n)^2 - d^2)/den;
    rho := 1 + n*(tau - 1);
    if rho eq 0 or tau eq 0 or tau^2 eq 1 then
        return false, Q!0, Q!0;
    end if;
    return true, rho, tau;
end function;

function CoverValues(rho, sigma, tau)
    F0 := rho*sigma*tau;
    F1 := (1 + rho)*(1 + sigma)*(1 + tau);
    F2 := rho*(1 + rho)*(rho + sigma)*(rho + tau);
    F3 := sigma*(1 + sigma)*(rho + sigma)*(sigma + tau);
    F4 := tau*(1 + tau)*(rho + tau)*(sigma + tau);
    return F0, F1, F2, F3, F4;
end function;

function RatHeight(q)
    return Max(Abs(Z!Numerator(q)), Z!Denominator(q));
end function;

function TupleHeight(tup)
    return Max([ RatHeight(Q!u) : u in tup ]);
end function;

function TupleKey(tup)
    roots := Sort([ Q!tup[i]^2 : i in [1..4] ]);
    return Sprint(roots);
end function;

function IrreducibleFrobeniusCertificateFromTuple(tup)
    a := Q!tup[1]; b := Q!tup[2]; c := Q!tup[3]; d := Q!tup[4];
    f := x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2);
    C2 := HyperellipticCurve(f);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97] do
        try
            fp := ChangeRing(f, GF(p));
            if Discriminant(fp) eq 0 then
                continue;
            end if;
            Lp := LPolynomial(ChangeRing(C2, GF(p)));
            fac := Factorization(Lp);
            if #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4 then
                return true, p, Lp;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false, 0, P!0;
end function;

params := RationalParametersOfHeight(dn_height);

print "M2248 q^2=1 diagonal search";
print "dn_height", dn_height, "params", #params;
print "identity", "sigma=rho/tau, F0=rho^2, C=rho^2/d^2, F2=rho^2*F1";

dn_checked := 0;
valid_surface := 0;
f1_square := 0;
f4_square := 0;
both_square := 0;
simple_certified := 0;
seen := {};
reports := 0;

for d in params do
    if d eq 0 then continue; end if;
    for n in params do
        if n eq 0 then continue; end if;
        dn_checked +:= 1;
        ok, rho, tau := DSurfaceParam(d, n);
        if not ok then
            continue;
        end if;
        sigma := rho/tau;
        if rho eq -1 or tau eq -1 or sigma eq -1 or rho eq tau or sigma eq rho or sigma eq tau then
            continue;
        end if;
        valid_surface +:= 1;

        F0, F1, F2, F3, F4 := CoverValues(rho, sigma, tau);
        assert F0 eq rho^2;
        assert F2 eq rho^2*F1;

        Cval := (sigma^2 - rho^2)/(sigma^2 - 1);
        assert Cval eq rho^2/d^2;

        ok1, r1 := IsNonzeroSquareQ(F1);
        ok4, r4 := IsNonzeroSquareQ(F4);
        if ok1 then f1_square +:= 1; end if;
        if ok4 then f4_square +:= 1; end if;
        if ok1 and ok4 then
            both_square +:= 1;
            ok3, r3 := IsNonzeroSquareQ(F3);
            assert ok3;
            c := rho/d;
            tup := [rho, Q!1, c, d];
            if #(Set([ u^2 : u in tup ])) ne 4 then
                continue;
            end if;
            key := TupleKey(tup);
            if key in seen then
                continue;
            end if;
            Include(~seen, key);
            ok_simple, pcert, Lp := IrreducibleFrobeniusCertificateFromTuple(tup);
            if ok_simple then
                simple_certified +:= 1;
            end if;
            if reports lt max_reports then
                reports +:= 1;
                print "Q1_FULL",
                      "d", d,
                      "n", n,
                      "rho", rho,
                      "sigma", sigma,
                      "tau", tau,
                      "tuple_height", TupleHeight(tup),
                      "simple_cert", ok_simple,
                      "cert_prime", pcert,
                      "tuple", tup;
            end if;
        end if;
    end for;
end for;

print "dn_checked", dn_checked;
print "valid_surface", valid_surface;
print "f1_square", f1_square;
print "f4_square", f4_square;
print "both_square", both_square;
print "unique_full", #seen;
print "simple_certified", simple_certified;
print "DONE";
