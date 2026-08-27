//////////////////////////////////////////////////////////////////////
//  Agent FRONT 3: HPL simultaneous pair-fiber sharpened computation.
//
//  Fixed at the normalized HPL (rho,tau,d) on the D-square + F0-square
//  surface.  For each pair among F1,F2,F4 this script:
//
//    * builds the reciprocal genus-one quartic obtained by forcing the
//      first conic and imposing the second;
//    * converts it to an elliptic curve with the HPL point as base point;
//    * records torsion/minimal-model data, with optional RankBounds;
//    * searches rational slopes on the first conic and tests the residual
//      C-square and remaining full-cover conditions.
//////////////////////////////////////////////////////////////////////

Q := Rationals();
Z := Integers();
R<m> := PolynomialRing(Q);
Px<x> := PolynomialRing(Q);

if not assigned slope_height then
    slope_height := 120;
elif Type(slope_height) eq MonStgElt then
    slope_height := StringToInteger(slope_height);
end if;

if not assigned do_rank then
    do_rank := false;
elif Type(do_rank) eq MonStgElt then
    do_rank := do_rank eq "true" or do_rank eq "1";
end if;

if not assigned max_reports then
    max_reports := 20;
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

function RatHeight(q)
    return Max(Abs(Z!Numerator(q)), Z!Denominator(q));
end function;

function TupleHeight(tup)
    return Max([ RatHeight(Q!u) : u in tup ]);
end function;

function TupleKey(tup)
    roots := Sort([ Q!tup[i]^2 : i in [1..#tup] ]);
    return Sprint(roots);
end function;

function CoeffOrZero(f, i)
    if i gt Degree(f) then
        return Q!0;
    end if;
    return Coefficient(f, i);
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

function PairQuartic(A1, B1, A2, B2, q0, y0)
    K := A1*B1;
    den := m^2 - K;
    num := q0*m^2 - 2*y0*m + q0*K;
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

function OtherConicPoint(A, B, q0, y0, s)
    K := A*B;
    den := s^2 - K;
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
    return false, 0, Px!0;
end function;

procedure TestQ(~stats, pair_label, q, slope_label, rho, tau, dH, forced)
    sigma := q^2*rho/tau;
    sigma_key := Sprint(sigma);
    if IsDefined(stats`seen_sigma, sigma_key) then
        return;
    end if;
    stats`seen_sigma[sigma_key] := true;
    stats`pair_points +:= 1;

    if sigma eq 0 or sigma^2 eq 1 or sigma eq -1 or
       sigma eq -rho or sigma eq -tau then
        return;
    end if;

    Cden := sigma^2 - 1;
    if Cden eq 0 then
        return;
    end if;
    Cval := (sigma^2 - rho^2)/Cden;
    okC, c := IsNonzeroSquareQ(Cval);
    if okC then
        tup := [rho, Q!1, c, dH];
        if #(Set([ u^2 : u in tup ])) eq 4 then
            key := TupleKey(tup);
            if key notin stats`seen_c then
                Include(~stats`seen_c, key);
                stats`c_square +:= 1;
                if stats`reports lt stats`max_reports then
                    stats`reports +:= 1;
                    print "C_SQUARE", pair_label,
                          "slope", slope_label,
                          "forced", forced,
                          "q_height", RatHeight(q),
                          "tuple_height", TupleHeight(tup),
                          "q", q,
                          "sigma", sigma;
                end if;
            end if;
            F0, F1, F2, F3, F4 := CoverValues(rho, sigma, tau);
            ok0, r0 := IsNonzeroSquareQ(F0);
            ok1, r1 := IsNonzeroSquareQ(F1);
            ok2, r2 := IsNonzeroSquareQ(F2);
            ok3, r3 := IsNonzeroSquareQ(F3);
            ok4, r4 := IsNonzeroSquareQ(F4);
            if ok0 and ok1 and ok2 and ok3 and ok4 then
                if key notin stats`seen_full then
                    Include(~stats`seen_full, key);
                    stats`full_cover +:= 1;
                    ok_simple, pcert, Lp := IrreducibleFrobeniusCertificateFromTuple(tup);
                    if ok_simple then
                        stats`simple_certified +:= 1;
                    end if;
                    print "FULL_COVER", pair_label,
                          "slope", slope_label,
                          "forced", forced,
                          "q_height", RatHeight(q),
                          "tuple_height", TupleHeight(tup),
                          "simple_cert", ok_simple,
                          "cert_prime", pcert,
                          "tuple", tup;
                end if;
            end if;
        end if;
    end if;
end procedure;

pairs := [
    <"F1","F2">,
    <"F1","F4">,
    <"F2","F4">
];

rhoH := Q!58466134224 / Q!53109477625;
tauH := Q!307598400 / Q!352612321;
dH := Q!72946054224 / Q!53109477625;
qH := tauH;

slopes := RationalParametersOfHeight(slope_height);

print "AGENT M2248 HPL pair-fiber search";
print "slope_height", slope_height, "num_slopes", #slopes, "do_rank", do_rank;
print "rhoH", rhoH;
print "tauH", tauH;
print "dH", dH;
print "qH", qH;

StatsFmt := recformat<
    seen_sigma,
    seen_c,
    seen_full,
    pair_points,
    c_square,
    full_cover,
    simple_certified,
    reports,
    max_reports
>;

for pair in pairs do
    first := pair[1];
    second := pair[2];
    pair_label := first cat "/" cat second;
    print "PAIR_START", pair_label;

    ok1, A1, B1 := TargetConicAB(rhoH, tauH, first);
    ok2, A2, B2 := TargetConicAB(rhoH, tauH, second);
    assert ok1 and ok2;

    okY1, y1 := IsNonzeroSquareQ(A1*(1 + B1*qH^2));
    okY2, y2 := IsNonzeroSquareQ(A2*(1 + B2*qH^2));
    assert okY1 and okY2;

    K := A1*B1;
    f := PairQuartic(A1, B1, A2, B2, qH, y1);
    sf := SquarefreePart(f);
    tangent_m := K*qH/y1;
    den0 := tangent_m^2 - K;
    P0W := y2*den0;
    C := HyperellipticCurve(f);
    P0 := C![tangent_m, P0W, 1];
    assert P0 in C;
    E, phi := EllipticCurve(C, P0);

    print "quartic_degree", Degree(f);
    print "squarefree_degree", Degree(sf);
    print "genus", Floor((Degree(sf)-1)/2);
    print "reciprocal_check_a0_eq_a4K2", CoeffOrZero(f,0) eq CoeffOrZero(f,4)*K^2;
    print "reciprocal_check_a1_eq_a3K", CoeffOrZero(f,1) eq CoeffOrZero(f,3)*K;
    print "hpl_tangent_slope", tangent_m;
    print "elliptic_model", E;
    try
        Emin := MinimalModel(E);
        print "minimal_model", Emin;
        print "minimal_a_invariants", aInvariants(Emin);
        print "minimal_discriminant", Discriminant(Emin);
        T, tors_map := TorsionSubgroup(Emin);
        print "torsion_invariants", Invariants(T);
        if do_rank then
            lo, hi := RankBounds(Emin);
            print "rank_bounds", lo, hi;
        end if;
    catch e
        print "minimal_or_rank_error", e`Object;
        T, tors_map := TorsionSubgroup(E);
        print "torsion_invariants_unminimal", Invariants(T);
        if do_rank then
            lo, hi := RankBounds(E);
            print "rank_bounds_unminimal", lo, hi;
        end if;
    end try;

    stats := rec<StatsFmt |
        seen_sigma := AssociativeArray(),
        seen_c := {},
        seen_full := {},
        pair_points := 0,
        c_square := 0,
        full_cover := 0,
        simple_certified := 0,
        reports := 0,
        max_reports := max_reports
    >;

    TestQ(~stats, pair_label, qH, "HPL", rhoH, tauH, dH, true);

    for s in slopes do
        ok_q, q, y := OtherConicPoint(A1, B1, qH, y1, s);
        if not ok_q then
            continue;
        end if;
        okSecond, y_second := IsNonzeroSquareQ(A2*(1 + B2*q^2));
        if okSecond then
            TestQ(~stats, pair_label, q, Sprint(s), rhoH, tauH, dH, false);
        end if;
    end for;

    print "PAIR_DONE", pair_label;
    print "pair_points", stats`pair_points;
    print "c_square_unique", stats`c_square;
    print "full_cover", stats`full_cover;
    print "simple_certified", stats`simple_certified;
end for;

print "DONE";
