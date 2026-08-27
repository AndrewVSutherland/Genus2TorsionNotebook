//////////////////////////////////////////////////////////////////////
//  Fixed-rho HPL elliptic-fibration search for M(2,2,4,8).
//
//  The corrected fixed-rho square condition
//
//      D = d^2 = (tau^2-rho^2)/(tau^2-1)
//
//  is not a conic in tau.  It is the genus-one quartic
//
//      Y^2 = (tau^2-rho^2)(tau^2-1),     Y = d*(tau^2-1).
//
//  For the normalized HPL witness this elliptic curve has rank 2.
//  This script searches bounded integer combinations of natural rational
//  points on that elliptic model.  For each new tau it keeps
//
//      sigma*tau/rho = sigma0*tau0/rho0
//
//  fixed, so the intermediate condition F0=rho*sigma*tau is automatic.
//  It then tests C=c^2 and the remaining full-cover square conditions.
//
//  Typical run:
//
//      magma -b coeff_bound:=8 max_hits:=20 \
//          code/m2248_hpl_fixedrho_elliptic_search.m \
//          > data/m2248_hpl_fixedrho_elliptic_summary.txt
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

if not assigned max_square_reports then
    max_square_reports := 20;
elif Type(max_square_reports) eq MonStgElt then
    max_square_reports := StringToInteger(max_square_reports);
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

// HPL normalized full witness.  The branch tuple is obtained from
// [1,b,c,d] by dividing through by b, so the second branch value is 1.
rho0 := Q!58466134224 / Q!53109477625;
sigma0 := Q!719363573659505664 / Q!749082246897952705;
tau0 := Q!307598400 / Q!352612321;
cB := Q!58466134224 / Q!30294861575;
dB := Q!72946054224 / Q!53109477625;

F00, F10, F20, F30, F40 := CoverValues(rho0, sigma0, tau0);
okF0, qF0 := IsNonzeroSquareQ(F00);
assert okF0;
okqbar, qbar0 := IsNonzeroSquareQ(sigma0*tau0/rho0);
assert okqbar;
assert #M2248WitnessesForTupleAllPermutations([rho0, Q!1, cB, dB], true) gt 0;

R<t> := PolynomialRing(Q);
quartic := (t^2-rho0^2)*(t^2-1);
C := HyperellipticCurve(quartic);
Y0 := dB*(tau0^2-1);
base_pt := C![tau0, Y0, 1];
assert base_pt in C;

E, phi := EllipticCurve(C, base_pt);
psi := Inverse(phi);

rank_lo, rank_hi := RankBounds(E);
Tors, tors_map := TorsionSubgroup(E);

print "HPL fixed-rho elliptic search";
print "coeff_bound", coeff_bound, "max_hits", max_hits;
print "rho0", rho0;
print "sigma0", sigma0;
print "tau0", tau0;
print "cB", cB;
print "dB", dB;
print "qbar0", qbar0;
print "elliptic_curve", E;
print "rank_bounds", rank_lo, rank_hi;
print "torsion_invariants", Invariants(Tors);

pointsC := [
    <"rho", C![rho0, Q!0, 1]>,
    <"negrho", C![-rho0, Q!0, 1]>,
    <"one", C![Q!1, Q!0, 1]>,
    <"negone", C![-Q!1, Q!0, 1]>,
    <"hpl_conj", C![tau0, -Y0, 1]>
];

candidates := [];
for item in pointsC do
    Append(~candidates, <item[1], phi(item[2])>);
end for;

for item in candidates do
    print "candidate", item[1], item[2];
end for;

function PullbackTauD(Qe)
    try
        Pc := psi(Qe);
    catch e
        return false, Q!0, Q!0, Q!0;
    end try;

    if Pc[3] eq 0 then
        return false, Q!0, Q!0, Q!0;
    end if;
    tau := Q!(Pc[1]/Pc[3]);
    Y := Q!(Pc[2]/(Pc[3]^2));
    if tau eq 0 or tau^2 eq 1 then
        return false, Q!0, Q!0, Q!0;
    end if;
    dbranch := Y/(tau^2-1);
    if dbranch eq 0 then
        return false, Q!0, Q!0, Q!0;
    end if;
    return true, tau, Y, dbranch;
end function;

checked := 0;
distinct_E := 0;
affine_pullbacks := 0;
valid_cover_coords := 0;
c_square := 0;
cover_square := 0;
a2244_sources := 0;
exact_sources := 0;
simple_hits := 0;
square_reports := 0;
hits := [];
seenE := {};
seenTup := {};

// Include single-generator multiples and pairwise rank-lattice probes.
jobs := [];
for i in [1..#candidates] do
    Append(~jobs, <i, i>);
end for;
for i in [1..#candidates] do
    for j in [i+1..#candidates] do
        Append(~jobs, <i, j>);
    end for;
end for;

for job in jobs do
    i := job[1];
    j := job[2];
    name_i := candidates[i][1];
    name_j := candidates[j][1];
    Pi := candidates[i][2];
    Pj := candidates[j][2];

    for ai in [-coeff_bound..coeff_bound] do
        for bj in [-coeff_bound..coeff_bound] do
            if ai eq 0 and bj eq 0 then
                continue;
            end if;
            if i eq j and bj ne 0 then
                continue;
            end if;

            checked +:= 1;
            Qe := ai*Pi + bj*Pj;
            if Qe eq E!0 then
                continue;
            end if;
            ekey := Sprint(Qe);
            if ekey in seenE then
                continue;
            end if;
            Include(~seenE, ekey);
            distinct_E +:= 1;

            ok_pb, tau, Y, dbranch := PullbackTauD(Qe);
            if not ok_pb then
                continue;
            end if;
            affine_pullbacks +:= 1;

            // Keep sigma*tau/rho equal to its HPL value.  Then F0 is
            // automatically rho0^2*qbar0^2, hence square.
            sigma := qbar0^2 * rho0 / tau;
            if sigma eq 0 or sigma^2 eq 1 then
                continue;
            end if;
            Cden := sigma^2 - 1;
            if Cden eq 0 then
                continue;
            end if;
            Cval := (sigma^2 - rho0^2)/Cden;
            valid_cover_coords +:= 1;

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

            if square_reports lt max_square_reports then
                square_reports +:= 1;
                print "C_SQUARE",
                      "pair", name_i, name_j,
                      "coeffs", ai, bj,
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
            Append(~hits, <name_i, name_j, ai, bj, tau, sigma, tup, simple, pcert, Lp, #witnesses>);

            print "HIT",
                  "pair", name_i, name_j,
                  "coeffs", ai, bj,
                  "tau", tau,
                  "sigma", sigma,
                  "tuple", tup,
                  "tuple_height", TupleHeight(tup),
                  "witnesses", #witnesses,
                  "simple", simple,
                  "pcert", pcert,
                  "L", Lp;

            if #hits ge max_hits then
                break job;
            end if;
        end for;
    end for;
end for;

print "DONE";
print "checked", checked;
print "distinct_E", distinct_E;
print "affine_pullbacks", affine_pullbacks;
print "valid_cover_coords", valid_cover_coords;
print "c_square", c_square;
print "cover_square", cover_square;
print "a2244_sources", a2244_sources;
print "exact_sources", exact_sources;
print "simple_hits", simple_hits;
print "hits", #hits;
