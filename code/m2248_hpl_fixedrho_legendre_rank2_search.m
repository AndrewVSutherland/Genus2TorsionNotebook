//////////////////////////////////////////////////////////////////////
//  Fixed-rho HPL Legendre rank-2 lattice search for M(2,2,4,8).
//
//  A 2-descent on the Legendre model
//
//      E0: V^2 = U(U-1)(U-rho0^2)
//
//  with torsion and the HPL point removed finds a rational point on a
//  remaining cover.  Its image is an independent rank direction P1.
//  This script searches small combinations
//
//      a*P0 + b*P1 + T,  T in E0[2],
//
//  keeps only points with U in Q^2, lifts them back to the fixed-rho
//  quartic, and then applies the residual M(2,2,4,8) square tests.
//////////////////////////////////////////////////////////////////////

load "code/m2248_sieve.m";

if not assigned coeff_bound then
    coeff_bound := 8;
elif Type(coeff_bound) eq MonStgElt then
    coeff_bound := StringToInteger(coeff_bound);
end if;

if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;

if not assigned max_lift_reports then
    max_lift_reports := 80;
elif Type(max_lift_reports) eq MonStgElt then
    max_lift_reports := StringToInteger(max_lift_reports);
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

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
Y0 := dB*(tau0^2-1);
V0 := tau0*Y0;

okqbar, qbar0 := IsNonzeroSquareQ(sigma0*tau0/rho0);
assert okqbar;
assert #M2248WitnessesForTupleAllPermutations([rho0, Q!1, cB, dB], true) gt 0;

E0 := EllipticCurve([Q!0, -(Q!1+rho0^2), Q!0, rho0^2, Q!0]);
P0 := E0![tau0^2, V0, 1];

// Independent point found by 2-descent after removing torsion and P0.
P1 := E0![
    Q!185261445034489541 / Q!75208331264762500,
    Q!2019033491444667568388406032021 / Q!950637494391639082835996875000,
    Q!1
];
assert P0 in E0;
assert P1 in E0;

rank_lo, rank_hi := RankBounds(E0);
Tors, tors_map := TorsionSubgroup(E0);

torsion_pts := [
    <"O", E0!0>,
    <"T0", E0![Q!0, Q!0, 1]>,
    <"T1", E0![Q!1, Q!0, 1]>,
    <"Trho", E0![rho0^2, Q!0, 1]>
];

print "HPL fixed-rho Legendre rank-2 search";
print "coeff_bound", coeff_bound, "max_hits", max_hits;
print "rho0", rho0;
print "sigma0", sigma0;
print "tau0", tau0;
print "cB", cB;
print "dB", dB;
print "qbar0", qbar0;
print "legendre_curve", E0;
print "rank_bounds", rank_lo, rank_hi;
print "torsion_invariants", Invariants(Tors);
print "P0", P0;
print "P1", P1;
print "P1_source", "TwoDescent RemoveTorsion=true RemoveGens={P0}, cover 2 point (-67 : 181175348718000 : 75)";

checked := 0;
distinct_E0 := 0;
u_square := 0;
valid_lifts := 0;
c_square := 0;
cover_square := 0;
a2244_sources := 0;
exact_sources := 0;
simple_hits := 0;
lift_reports := 0;
hits := [];
seenE0 := {};
seenTup := {};

coeffs := [];
for radius in [0..coeff_bound] do
    for a in [-radius..radius] do
        for b in [-radius..radius] do
            if Max(Abs(a), Abs(b)) ne radius then
                continue;
            end if;
            if a eq 0 and b eq 0 then
                continue;
            end if;
            Append(~coeffs, <a,b>);
        end for;
    end for;
end for;

for coeff in coeffs do
    a := coeff[1];
    b := coeff[2];
    for tors in torsion_pts do
        checked +:= 1;
        Qe := a*P0 + b*P1 + tors[2];
        if Qe eq E0!0 then
            continue;
        end if;
        ekey := Sprint(Qe);
        if ekey in seenE0 then
            continue;
        end if;
        Include(~seenE0, ekey);
        distinct_E0 +:= 1;

        if Qe[3] eq 0 then
            continue;
        end if;
        U := Q!(Qe[1]/Qe[3]);
        V := Q!(Qe[2]/Qe[3]);
        okU, tau_abs := IsNonzeroSquareQ(U);
        if not okU then
            continue;
        end if;
        u_square +:= 1;

        for tau in [tau_abs, -tau_abs] do
            if tau eq 0 or tau^2 eq 1 then
                continue;
            end if;
            Y := V/tau;
            dbranch := Y/(tau^2-1);
            if dbranch eq 0 then
                continue;
            end if;
            valid_lifts +:= 1;

            sigma := qbar0^2 * rho0 / tau;
            if sigma eq 0 or sigma^2 eq 1 then
                continue;
            end if;
            Cden := sigma^2 - 1;
            if Cden eq 0 then
                continue;
            end if;
            Cval := (sigma^2 - rho0^2)/Cden;
            okC, cbranch := IsNonzeroSquareQ(Cval);
            if not okC then
                continue;
            end if;
            c_square +:= 1;

            tup := [rho0, Q!1, cbranch, dbranch];
            if #(Set([ u^2 : u in tup ])) ne 4 then
                continue;
            end if;
            tkey := TupleKey(tup);
            if tkey in seenTup then
                continue;
            end if;
            Include(~seenTup, tkey);

            if lift_reports lt max_lift_reports then
                lift_reports +:= 1;
                print "C_SQUARE",
                      "coeffs", a, b,
                      "tors", tors[1],
                      "tau", tau,
                      "sigma", sigma,
                      "tuple_height", TupleHeight(tup);
            end if;

            F0, F1, F2, F3, F4 := CoverValues(rho0, sigma, tau);
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

            if M2248OrderedTupleHasA2244Squares(tup) then
                a2244_sources +:= 1;
            end if;

            witnesses := M2248WitnessesForTupleAllPermutations(tup, true);
            if #witnesses eq 0 then
                continue;
            end if;
            exact_sources +:= 1;

            simple, pcert, Lp := IrreducibleFrobeniusCertificateFromTuple(tup);
            if simple then
                simple_hits +:= 1;
            end if;
            Append(~hits, <a, b, tors[1], tau, sigma, tup, simple, pcert, Lp, #witnesses>);

            print "HIT",
                  "coeffs", a, b,
                  "tors", tors[1],
                  "tau", tau,
                  "sigma", sigma,
                  "tuple", tup,
                  "tuple_height", TupleHeight(tup),
                  "witnesses", #witnesses,
                  "simple", simple,
                  "pcert", pcert,
                  "L", Lp;

            if #hits ge max_hits then
                break coeff;
            end if;
        end for;
    end for;
end for;

print "DONE";
print "checked", checked;
print "distinct_E0", distinct_E0;
print "u_square", u_square;
print "valid_lifts", valid_lifts;
print "c_square", c_square;
print "cover_square", cover_square;
print "a2244_sources", a2244_sources;
print "exact_sources", exact_sources;
print "simple_hits", simple_hits;
print "hits", #hits;
