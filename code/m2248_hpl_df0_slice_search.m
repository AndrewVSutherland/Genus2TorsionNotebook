//////////////////////////////////////////////////////////////////////
//  HPL-centered line slices on the D-square + F0-square surface.
//
//  The rational surface coordinates are:
//      d  = branch square root for D,
//      n  = slope in the conic parametrization of (rho,tau),
//      q  = square root parameter with sigma=q^2*rho/tau.
//
//  This keeps D=d^2 and F0=(q*rho)^2 automatic.  Instead of searching
//  small absolute d,n,q, this searches multiplicative line slices
//  through the high-height HPL point:
//
//      d = d0*(1 + ed*u),
//      n = n0*(1 + en*u),
//      q = q0*(1 + eq*u).
//
//  We then test C=c^2 and F1..F4 exactly.
//////////////////////////////////////////////////////////////////////

load "code/m2248_sieve.m";

if not assigned dir_bound then
    dir_bound := 3;
elif Type(dir_bound) eq MonStgElt then
    dir_bound := StringToInteger(dir_bound);
end if;

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

if not assigned max_c_reports then
    max_c_reports := 40;
elif Type(max_c_reports) eq MonStgElt then
    max_c_reports := StringToInteger(max_c_reports);
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

function TupleKey(tup)
    roots := Sort([ Q!tup[i]^2 : i in [1..4] ]);
    return Sprint(roots);
end function;

function RatHeight(q)
    return Max(Abs(Z!Numerator(q)), Z!Denominator(q));
end function;

function TupleHeight(tup)
    return Max([ RatHeight(Q!u) : u in tup ]);
end function;

rho0 := Q!58466134224 / Q!53109477625;
sigma0 := Q!719363573659505664 / Q!749082246897952705;
tau0 := Q!307598400 / Q!352612321;
cB := Q!58466134224 / Q!30294861575;
dB := Q!72946054224 / Q!53109477625;
qB := tau0;
nB := (rho0 - 1)/(tau0 - 1);

ok_hpl, rho_check, tau_check := DSurfaceParam(dB, nB);
assert ok_hpl and rho_check eq rho0 and tau_check eq tau0;
assert sigma0 eq qB^2*rho0/tau0;
assert #M2248WitnessesForTupleAllPermutations([rho0, Q!1, cB, dB], true) gt 0;

params := RationalParametersOfHeight(height);
print "HPL D-square + F0-square line-slice search";
print "dir_bound", dir_bound, "height", height, "parameters", #params;
print "base_heights", "d", RatHeight(dB), "n", RatHeight(nB), "q", RatHeight(qB);

checked := 0;
valid_slice := 0;
c_square := 0;
cover_square := 0;
exact_sources := 0;
simple_hits := 0;
c_reports := 0;
hits := [];
seen := {};
stop := false;

for ed in [-dir_bound..dir_bound] do
    if stop then break; end if;
    for en in [-dir_bound..dir_bound] do
        if stop then break; end if;
        for eqd in [-dir_bound..dir_bound] do
            if stop then break; end if;
            if ed eq 0 and en eq 0 and eqd eq 0 then
                continue;
            end if;

            for u in params do
                if u eq 0 then
                    continue;
                end if;
                checked +:= 1;

                d := dB*(1 + ed*u);
                n := nB*(1 + en*u);
                q := qB*(1 + eqd*u);
                if d eq 0 or n eq 0 or q eq 0 then
                    continue;
                end if;

                ok_surface, rho, tau := DSurfaceParam(d, n);
                if not ok_surface then
                    continue;
                end if;
                if rho eq -1 or tau eq -1 or rho eq tau then
                    continue;
                end if;
                valid_slice +:= 1;

                sigma := q^2*rho/tau;
                if sigma eq 0 or sigma^2 eq 1 or sigma eq -1 or sigma eq -rho or sigma eq -tau then
                    continue;
                end if;

                Cden := sigma^2 - 1;
                if Cden eq 0 then
                    continue;
                end if;
                Cval := (sigma^2 - rho^2)/Cden;
                okC, c := IsNonzeroSquareQ(Cval);
                if not okC then
                    continue;
                end if;
                tup := [rho, Q!1, c, d];
                if #(Set([ w^2 : w in tup ])) ne 4 then
                    continue;
                end if;
                c_square +:= 1;

                if c_reports lt max_c_reports then
                    c_reports +:= 1;
                    print "C_SQUARE",
                          "dir", [ed,en,eqd],
                          "u", u,
                          "rho", rho,
                          "sigma", sigma,
                          "tau", tau,
                          "tuple_height", TupleHeight(tup);
                end if;

                F0, F1, F2, F3, F4 := CoverValues(rho, sigma, tau);
                ok0, q0 := IsNonzeroSquareQ(F0);
                if not ok0 then continue; end if;
                ok1, q1 := IsNonzeroSquareQ(F1);
                if not ok1 then continue; end if;
                ok2, q2 := IsNonzeroSquareQ(F2);
                if not ok2 then continue; end if;
                ok3, q3 := IsNonzeroSquareQ(F3);
                if not ok3 then continue; end if;
                ok4, q4 := IsNonzeroSquareQ(F4);
                if not ok4 then continue; end if;
                cover_square +:= 1;

                key := TupleKey(tup);
                if key in seen then
                    continue;
                end if;
                Include(~seen, key);

                witnesses := M2248WitnessesForTupleAllPermutations(tup, true);
                if #witnesses eq 0 then
                    continue;
                end if;
                exact_sources +:= 1;
                simple, pcert, Lp := IrreducibleFrobeniusCertificateFromTuple(tup);
                if simple then
                    simple_hits +:= 1;
                end if;
                Append(~hits, <ed,en,eqd,u,rho,sigma,tau,tup,simple,pcert,Lp,#witnesses>);
                print "HIT",
                      "dir", [ed,en,eqd],
                      "u", u,
                      "rho", rho,
                      "sigma", sigma,
                      "tau", tau,
                      "tuple", tup,
                      "tuple_height", TupleHeight(tup),
                      "witnesses", #witnesses,
                      "simple", simple,
                      "pcert", pcert,
                      "L", Lp;
                if #hits ge max_hits then
                    stop := true;
                    break;
                end if;
            end for;
        end for;
    end for;
end for;

print "DONE";
print "checked", checked;
print "valid_slice", valid_slice;
print "c_square", c_square;
print "cover_square", cover_square;
print "exact_sources", exact_sources;
print "simple_hits", simple_hits;
print "hits", #hits;
