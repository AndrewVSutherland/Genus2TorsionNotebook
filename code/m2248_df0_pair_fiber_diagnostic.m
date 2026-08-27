//////////////////////////////////////////////////////////////////////
//  Symbolic/experimental diagnostic for simultaneous Fi/Fj fibers.
//
//  On the D-square + F0-square surface, write
//
//      sigma = q^2*rho/tau.
//
//  For fixed (d,n), hence fixed (rho,tau), each of F1,F2,F4 is a
//  conic in q:
//
//      y^2 = A*(1+B*q^2).
//
//  If one conic has a rational point (q0,y0), parametrizing it by
//  slope m and imposing a second conic gives a quartic
//
//      W^2 = A2*((m^2-A1*B1)^2
//                + B2*(q0*m^2 - 2*y0*m + q0*A1*B1)^2).
//
//  This script records the squarefree degree/genus of these pair
//  fibers and searches small slopes for simultaneous Fi/Fj points.
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
    max_seeds_per_dn := 1;
elif Type(max_seeds_per_dn) eq MonStgElt then
    max_seeds_per_dn := StringToInteger(max_seeds_per_dn);
end if;

if not assigned max_pair_reports then
    max_pair_reports := 30;
elif Type(max_pair_reports) eq MonStgElt then
    max_pair_reports := StringToInteger(max_pair_reports);
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);
R<m> := PolynomialRing(Q);

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

function TargetConicAB(rho, tau, tgt)
    if tgt eq "F1" then
        return true, (1 + rho)*(1 + tau), rho/tau;
    elif tgt eq "F2" then
        return true, (1 + rho)*(rho + tau), Q!1/tau;
    elif tgt eq "F4" then
        return true, (1 + tau)*(rho + tau), rho/(tau^2);
    end if;
    return false, Q!0, Q!0;
end function;

function OtherConicPoint(A, B, q0, y0, s)
    den := s^2 - A*B;
    if den eq 0 then
        return false, Q!0, Q!0;
    end if;
    q1 := -2*s*(y0 - s*q0)/den - q0;
    y1 := y0 + s*(q1 - q0);
    if q1 eq 0 then
        return false, Q!0, Q!0;
    end if;
    return true, q1, y1;
end function;

function PairQuartic(A1, B1, A2, B2, q0, y0)
    den := m^2 - A1*B1;
    num := q0*m^2 - 2*y0*m + q0*A1*B1;
    return A2*(den^2 + B2*num^2);
end function;

function SquarefreePart(f)
    if f eq 0 then
        return R!0;
    end if;
    fac := Factorization(f);
    sf := R!1;
    for ff in fac do
        if ff[2] mod 2 eq 1 then
            sf *:= ff[1];
        end if;
    end for;
    return sf;
end function;

function GenusFromSquarefreeDegree(deg)
    if deg le 0 then
        return -1;
    end if;
    return Floor((deg - 1)/2);
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

pairs := [
    <"F1","F2">,
    <"F1","F4">,
    <"F2","F4">
];

dn_params := RationalParametersOfHeight(dn_height);
seed_params := RationalParametersOfHeight(seed_height);
slope_params := RationalParametersOfHeight(slope_height);

print "M2248 D+F0 simultaneous pair-fiber diagnostic";
print "dn_height", dn_height,
      "seed_height", seed_height,
      "slope_height", slope_height,
      "dn_params", #dn_params,
      "seed_params", #seed_params,
      "slope_params", #slope_params;

// HPL fiber degree/genus sanity check for each pair.
rhoH := Q!58466134224 / Q!53109477625;
sigmaH := Q!719363573659505664 / Q!749082246897952705;
tauH := Q!307598400 / Q!352612321;
dH := Q!72946054224 / Q!53109477625;
qH := tauH;
nH := (rhoH - 1)/(tauH - 1);
okH, rho_check, tau_check := DSurfaceParam(dH, nH);
assert okH and rho_check eq rhoH and tau_check eq tauH;
assert sigmaH eq qH^2*rhoH/tauH;

for pair in pairs do
    first := pair[1];
    second := pair[2];
    ok1H, A1H, B1H := TargetConicAB(rhoH, tauH, first);
    ok2H, A2H, B2H := TargetConicAB(rhoH, tauH, second);
    assert ok1H and ok2H;
    okYH, yH := IsNonzeroSquareQ(A1H*(1 + B1H*qH^2));
    assert okYH;
    fH := PairQuartic(A1H, B1H, A2H, B2H, qH, yH);
    sfH := SquarefreePart(fH);
    degH := Degree(sfH);
    print "HPL_PAIR_FIBER", first, second,
          "degree", Degree(fH),
          "squarefree_degree", degH,
          "genus", GenusFromSquarefreeDegree(degH);
end for;

for pair in pairs do
    first := pair[1];
    second := pair[2];
    print "PAIR_START", first, second;

    dn_checked := 0;
    valid_surface := 0;
    seed_dns := 0;
    seeds := 0;
    genus_counts := AssociativeArray(Integers());
    degree_counts := AssociativeArray(Integers());
    pair_slope_points := 0;
    c_square := 0;
    all_cover_square := 0;
    reports := 0;
    seen_c := {};

    for d in dn_params do
        if d eq 0 then continue; end if;
        for n in dn_params do
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

            ok1, A1, B1 := TargetConicAB(rho, tau, first);
            ok2, A2, B2 := TargetConicAB(rho, tau, second);
            if not (ok1 and ok2) then
                continue;
            end if;

            local_seeds := [];
            for q0 in seed_params do
                if q0 eq 0 then continue; end if;
                sigma0 := q0^2*rho/tau;
                if sigma0 eq 0 or sigma0^2 eq 1 then
                    continue;
                end if;
                okSeed, y0 := IsNonzeroSquareQ(A1*(1 + B1*q0^2));
                if okSeed then
                    Append(~local_seeds, <q0, y0>);
                    seeds +:= 1;
                    if #local_seeds ge max_seeds_per_dn then
                        break;
                    end if;
                end if;
            end for;
            if #local_seeds eq 0 then
                continue;
            end if;
            seed_dns +:= 1;

            for seed in local_seeds do
                q0 := seed[1];
                y0 := seed[2];
                f := PairQuartic(A1, B1, A2, B2, q0, y0);
                sf := SquarefreePart(f);
                deg := Degree(sf);
                genus := GenusFromSquarefreeDegree(deg);
                if not IsDefined(degree_counts, deg) then
                    degree_counts[deg] := 0;
                end if;
                degree_counts[deg] +:= 1;
                if not IsDefined(genus_counts, genus) then
                    genus_counts[genus] := 0;
                end if;
                genus_counts[genus] +:= 1;

                for s in slope_params do
                    ok_q, q, y := OtherConicPoint(A1, B1, q0, y0, s);
                    if not ok_q then
                        continue;
                    end if;
                    okSecond, y2 := IsNonzeroSquareQ(A2*(1 + B2*q^2));
                    if not okSecond then
                        continue;
                    end if;
                    pair_slope_points +:= 1;

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
                    if okC then
                        tup := [rho, Q!1, c, d];
                        if #(Set([ u^2 : u in tup ])) eq 4 then
                            key := TupleKey(tup);
                            if key notin seen_c then
                                Include(~seen_c, key);
                                c_square +:= 1;
                                if reports lt max_pair_reports then
                                    reports +:= 1;
                                    print "PAIR_C_SQUARE",
                                          first, second,
                                          "d", d,
                                          "n", n,
                                          "q0", q0,
                                          "slope", s,
                                          "q", q,
                                          "rho", rho,
                                          "sigma", sigma,
                                          "tau", tau,
                                          "tuple_height", TupleHeight(tup),
                                          "genus", genus;
                                end if;
                            end if;
                            F0, F1, F2, F3, F4 := CoverValues(rho, sigma, tau);
                            ok0, r0 := IsNonzeroSquareQ(F0);
                            okF1, rF1 := IsNonzeroSquareQ(F1);
                            okF2, rF2 := IsNonzeroSquareQ(F2);
                            okF3, rF3 := IsNonzeroSquareQ(F3);
                            okF4, rF4 := IsNonzeroSquareQ(F4);
                            if ok0 and okF1 and okF2 and okF3 and okF4 then
                                all_cover_square +:= 1;
                                print "PAIR_FULL_COVER",
                                      first, second,
                                      "d", d,
                                      "n", n,
                                      "q", q,
                                      "tuple", tup;
                            end if;
                        end if;
                    end if;
                end for;
            end for;
        end for;
    end for;

    print "PAIR_DONE", first, second;
    print "dn_checked", dn_checked;
    print "valid_surface", valid_surface;
    print "seed_dns", seed_dns;
    print "seeds", seeds;
    deg_keys := Sort(SetToSequence(Keys(degree_counts)));
    genus_keys := Sort(SetToSequence(Keys(genus_counts)));
    for deg in deg_keys do
        print "DEGREE_COUNT", first, second, deg, degree_counts[deg];
    end for;
    for genus in genus_keys do
        print "GENUS_COUNT", first, second, genus, genus_counts[genus];
    end for;
    print "pair_slope_points", pair_slope_points;
    print "c_square_unique", c_square;
    print "all_cover_square", all_cover_square;
end for;

print "DONE";
