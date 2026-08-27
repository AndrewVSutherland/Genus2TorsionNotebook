//////////////////////////////////////////////////////////////////////
//  M(2,2,4,8) search on a rational D-square + F0-square surface.
//
//  The condition
//
//      D = d^2 = (tau^2-rho^2)/(tau^2-1)
//
//  is awkward as a fixed-rho genus-one curve, but it is rational as a
//  surface.  For fixed d, the equation
//
//      rho^2 = (1-d^2)*tau^2 + d^2
//
//  is a conic in (rho,tau), with boundary point (rho,tau)=(1,1).
//  A line rho = 1 + n*(tau-1) gives a two-parameter parametrization.
//
//  Then set
//
//      sigma = q^2*rho/tau,
//
//  so F0=rho*sigma*tau=(q*rho)^2 is automatic.  We test C=c^2 and
//  the four remaining full-cover square conditions F1..F4.
//
//  Typical run:
//      magma -b height:=10 max_hits:=20 code/m2248_df0_surface_search.m \
//          > data/m2248_df0_surface_h10.txt
//////////////////////////////////////////////////////////////////////

load "code/m2248_sieve.m";

if not assigned height then
    height := 8;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;

if not assigned progress_interval then
    progress_interval := 500000;
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

// Correct normalized HPL witness, used only as a parametrization sanity check.
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
print "M2248 D-square + F0-square rational surface search";
print "height", height, "parameters", #params, "max_hits", max_hits;
print "HPL_param_heights", "d", RatHeight(dB), "n", RatHeight(nB), "q", RatHeight(qB);

checked := 0;
valid_surface := 0;
c_square := 0;
cover_square := 0;
exact_sources := 0;
simple_hits := 0;
hits := [];
seen := {};
stop := false;

for d in params do
    if stop then break; end if;
    if d eq 0 then
        continue;
    end if;

    for n in params do
        if stop then break; end if;
        if n eq 0 then
            continue;
        end if;

        ok_surface, rho, tau := DSurfaceParam(d, n);
        if not ok_surface then
            continue;
        end if;
        if rho eq -1 or tau eq -1 or rho eq tau then
            continue;
        end if;
        valid_surface +:= 1;

        for q in params do
            if q eq 0 then
                continue;
            end if;
            checked +:= 1;

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
            if #(Set([ u^2 : u in tup ])) ne 4 then
                continue;
            end if;
            c_square +:= 1;

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
            Append(~hits, <d,n,q,rho,sigma,tau,tup,simple,pcert,Lp,#witnesses>);
            print "HIT",
                  "d", d,
                  "n", n,
                  "q", q,
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

            if progress_interval gt 0 and checked mod progress_interval eq 0 then
                print "checked", checked,
                      "valid_surface", valid_surface,
                      "c_square", c_square,
                      "cover_square", cover_square,
                      "hits", #hits;
            end if;
        end for;
    end for;
end for;

print "DONE";
print "checked", checked;
print "valid_surface", valid_surface;
print "c_square", c_square;
print "cover_square", cover_square;
print "exact_sources", exact_sources;
print "simple_hits", simple_hits;
print "hits", #hits;
