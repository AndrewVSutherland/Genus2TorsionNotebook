//////////////////////////////////////////////////////////////////////
//  M(2,2,4,8) search forcing D, F0, and F1.
//
//  Coordinates:
//      D=d^2 is forced by the rational surface parametrization
//      F0 is forced by sigma=q^2*rho/tau
//
//  For fixed (d,n), hence fixed (rho,tau), the condition
//
//      F1 = (1+rho)(1+sigma)(1+tau)
//
//  with sigma=q^2*rho/tau is the conic
//
//      y^2 = A*(1+B*q^2),  A=(1+rho)(1+tau), B=rho/tau.
//
//  This script first searches for a small seed point (q0,y0) on that
//  conic, then parametrizes the conic by lines of slope m through the
//  seed.  Thus D, F0, and F1 are automatic, and only C,F2,F3,F4 remain.
//////////////////////////////////////////////////////////////////////

load "code/m2248_sieve.m";

if not assigned dn_height then
    dn_height := 8;
elif Type(dn_height) eq MonStgElt then
    dn_height := StringToInteger(dn_height);
end if;

if not assigned seed_height then
    seed_height := 6;
elif Type(seed_height) eq MonStgElt then
    seed_height := StringToInteger(seed_height);
end if;

if not assigned slope_height then
    slope_height := 8;
elif Type(slope_height) eq MonStgElt then
    slope_height := StringToInteger(slope_height);
end if;

if not assigned max_seeds_per_dn then
    max_seeds_per_dn := 2;
elif Type(max_seeds_per_dn) eq MonStgElt then
    max_seeds_per_dn := StringToInteger(max_seeds_per_dn);
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

function OtherConicPoint(A, B, q0, y0, m)
    den := m^2 - A*B;
    if den eq 0 then
        return false, Q!0, Q!0;
    end if;
    q1 := -2*m*(y0 - m*q0)/den - q0;
    y1 := y0 + m*(q1 - q0);
    if q1 eq 0 then
        return false, Q!0, Q!0;
    end if;
    return true, q1, y1;
end function;

// HPL sanity check.
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
F00, F10, F20, F30, F40 := CoverValues(rho0, sigma0, tau0);
okF1HPL, yHPL := IsNonzeroSquareQ(F10);
assert okF1HPL;
assert #M2248WitnessesForTupleAllPermutations([rho0, Q!1, cB, dB], true) gt 0;

dn_params := RationalParametersOfHeight(dn_height);
seed_params := RationalParametersOfHeight(seed_height);
slope_params := RationalParametersOfHeight(slope_height);

print "M2248 D+F0+F1 conic search";
print "dn_height", dn_height,
      "seed_height", seed_height,
      "slope_height", slope_height,
      "dn_params", #dn_params,
      "seed_params", #seed_params,
      "slope_params", #slope_params,
      "max_hits", max_hits;
print "HPL_param_heights", "d", RatHeight(dB), "n", RatHeight(nB), "q", RatHeight(qB);

dn_checked := 0;
valid_surface := 0;
f1_seed_dns := 0;
f1_seeds := 0;
f1_param_points := 0;
c_square := 0;
f234_square := 0;
exact_sources := 0;
simple_hits := 0;
c_reports := 0;
hits := [];
seen := {};
stop := false;

for d in dn_params do
    if stop then break; end if;
    if d eq 0 then continue; end if;
    for n in dn_params do
        if stop then break; end if;
        if n eq 0 then continue; end if;
        dn_checked +:= 1;

        ok_surface, rho, tau := DSurfaceParam(d, n);
        if not ok_surface then
            continue;
        end if;
        if rho eq -1 or tau eq -1 or rho eq tau then
            continue;
        end if;
        valid_surface +:= 1;

        A := (1 + rho)*(1 + tau);
        B := rho/tau;
        seeds := [];
        for q0 in seed_params do
            if q0 eq 0 then continue; end if;
            sigma0_local := q0^2*rho/tau;
            if sigma0_local eq 0 or sigma0_local^2 eq 1 then
                continue;
            end if;
            F1val := A*(1 + B*q0^2);
            okSeed, y0 := IsNonzeroSquareQ(F1val);
            if okSeed then
                Append(~seeds, <q0, y0>);
                f1_seeds +:= 1;
                if #seeds ge max_seeds_per_dn then
                    break;
                end if;
            end if;
        end for;
        if #seeds eq 0 then
            continue;
        end if;
        f1_seed_dns +:= 1;

        for seed in seeds do
            q0 := seed[1];
            y0 := seed[2];
            for m in slope_params do
                ok_q, q, y := OtherConicPoint(A, B, q0, y0, m);
                if not ok_q then
                    continue;
                end if;
                f1_param_points +:= 1;

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

                if c_reports lt max_c_reports then
                    c_reports +:= 1;
                    print "C_SQUARE",
                          "d", d,
                          "n", n,
                          "seed_q", q0,
                          "slope", m,
                          "q", q,
                          "rho", rho,
                          "sigma", sigma,
                          "tau", tau,
                          "tuple_height", TupleHeight(tup);
                end if;

                F0, F1, F2, F3, F4 := CoverValues(rho, sigma, tau);
                ok0, r0 := IsNonzeroSquareQ(F0);
                if not ok0 then continue; end if;
                ok1, r1 := IsNonzeroSquareQ(F1);
                if not ok1 then continue; end if;
                ok2, r2 := IsNonzeroSquareQ(F2);
                if not ok2 then continue; end if;
                ok3, r3 := IsNonzeroSquareQ(F3);
                if not ok3 then continue; end if;
                ok4, r4 := IsNonzeroSquareQ(F4);
                if not ok4 then continue; end if;
                f234_square +:= 1;

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
                Append(~hits, <d,n,q0,m,q,rho,sigma,tau,tup,simple,pcert,Lp,#witnesses>);
                print "HIT",
                      "d", d,
                      "n", n,
                      "seed_q", q0,
                      "slope", m,
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
            end for;
        end for;
    end for;
end for;

print "DONE";
print "dn_checked", dn_checked;
print "valid_surface", valid_surface;
print "f1_seed_dns", f1_seed_dns;
print "f1_seeds", f1_seeds;
print "f1_param_points", f1_param_points;
print "c_square", c_square;
print "f234_square", f234_square;
print "exact_sources", exact_sources;
print "simple_hits", simple_hits;
print "hits", #hits;
