//////////////////////////////////////////////////////////////////////
//  Agent FRONT 3B: HPL pair-fiber rank/descent/point probe.
//
//  This is intentionally separate from the shared m2248_pair_* scripts.
//  It fixes the normalized HPL (rho,tau,d) and studies the three
//  simultaneous pair fibers F1/F2, F1/F4, F2/F4.
//
//  Optional -b parameters:
//      pair_index       0 for all, or 1..3
//      slope_height     rational slope height for direct pair-fiber search
//      multiple_bound   multiples of known rational pair-fiber points
//      cover_bound      Points(... : Bound := cover_bound) on 2-coverings
//      do_rank          true/false
//      do_twodescent    true/false
//      max_reports      cap detailed candidate prints per pair
//////////////////////////////////////////////////////////////////////

Q := Rationals();
Z := Integers();
R<m> := PolynomialRing(Q);
Px<x> := PolynomialRing(Q);

if not assigned pair_index then
    pair_index := 0;
elif Type(pair_index) eq MonStgElt then
    pair_index := StringToInteger(pair_index);
end if;

if not assigned slope_height then
    slope_height := 120;
elif Type(slope_height) eq MonStgElt then
    slope_height := StringToInteger(slope_height);
end if;

if not assigned multiple_bound then
    multiple_bound := 40;
elif Type(multiple_bound) eq MonStgElt then
    multiple_bound := StringToInteger(multiple_bound);
end if;

if not assigned cover_bound then
    cover_bound := 200;
elif Type(cover_bound) eq MonStgElt then
    cover_bound := StringToInteger(cover_bound);
end if;

if not assigned do_rank then
    do_rank := false;
elif Type(do_rank) eq MonStgElt then
    do_rank := do_rank eq "true" or do_rank eq "1";
end if;

if not assigned do_twodescent then
    do_twodescent := false;
elif Type(do_twodescent) eq MonStgElt then
    do_twodescent := do_twodescent eq "true" or do_twodescent eq "1";
end if;

if not assigned max_reports then
    max_reports := 24;
elif Type(max_reports) eq MonStgElt then
    max_reports := StringToInteger(max_reports);
end if;

function BoolString(b)
    return b select "true" else "false";
end function;

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            r := Q!num / Q!den;
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
    return true, Q!rtN / Q!rtD;
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

function QFromSlope(A, B, q0, y0, s)
    K := A*B;
    den := s^2 - K;
    if den eq 0 then
        return false, Q!0;
    end if;
    q1 := (q0*s^2 - 2*y0*s + q0*K)/den;
    if q1 eq 0 then
        return false, Q!0;
    end if;
    return true, q1;
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

StatsFmt := recformat<
    seen_sigma,
    seen_tuple,
    pair_points,
    c_square_unique,
    full_cover_unique,
    simple_certified,
    reports
>;

procedure TestQ(~stats, pair_label, q, source, rho, tau, dH)
    sigma := q^2*rho/tau;
    sigma_key := Sprint(sigma);
    if sigma_key in stats`seen_sigma then
        return;
    end if;
    Include(~stats`seen_sigma, sigma_key);
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
    if not okC then
        return;
    end if;

    tup := [rho, Q!1, c, dH];
    if #(Set([ u^2 : u in tup ])) ne 4 then
        return;
    end if;

    key := TupleKey(tup);
    if key in stats`seen_tuple then
        return;
    end if;
    Include(~stats`seen_tuple, key);
    stats`c_square_unique +:= 1;

    F0, F1, F2, F3, F4 := CoverValues(rho, sigma, tau);
    ok0, r0 := IsNonzeroSquareQ(F0);
    ok1, r1 := IsNonzeroSquareQ(F1);
    ok2, r2 := IsNonzeroSquareQ(F2);
    ok3, r3 := IsNonzeroSquareQ(F3);
    ok4, r4 := IsNonzeroSquareQ(F4);
    all_ok := ok0 and ok1 and ok2 and ok3 and ok4;

    if stats`reports lt max_reports then
        stats`reports +:= 1;
        print "C_SQUARE",
              pair_label,
              "source", source,
              "q_height", RatHeight(q),
              "tuple_height", TupleHeight(tup),
              "all_cover", BoolString(all_ok),
              "q", q,
              "sigma", sigma;
    end if;

    if all_ok then
        stats`full_cover_unique +:= 1;
        ok_simple, pcert, Lp := IrreducibleFrobeniusCertificateFromTuple(tup);
        if ok_simple then
            stats`simple_certified +:= 1;
        end if;
        print "FULL_COVER",
              pair_label,
              "source", source,
              "q_height", RatHeight(q),
              "tuple_height", TupleHeight(tup),
              "simple_cert", BoolString(ok_simple),
              "cert_prime", pcert,
              "tuple", tup;
    end if;
end procedure;

function SmallTorsionOrder(P, bound)
    if P eq Parent(P)!0 then
        return 1;
    end if;
    for n in [2..bound] do
        if n*P eq Parent(P)!0 then
            return n;
        end if;
    end for;
    return 0;
end function;

function TorsionPointsFromSubgroup(E)
    pts := [E!0];
    try
        G, mp := TorsionSubgroup(E);
        pts := [];
        for g in G do
            Append(~pts, mp(g));
        end for;
        return pts;
    catch e
        return [E!0];
    end try;
end function;

pairs := [
    <"F1","F2">,
    <"F1","F4">,
    <"F2","F4">
];

rhoH := Q!58466134224 / Q!53109477625;
sigmaH := Q!719363573659505664 / Q!749082246897952705;
tauH := Q!307598400 / Q!352612321;
cH := Q!58466134224 / Q!30294861575;
dH := Q!72946054224 / Q!53109477625;
qH := tauH;

slopes := RationalParametersOfHeight(slope_height);

print "AGENT M2248 HPL pair rank probe";
print "pair_index", pair_index,
      "slope_height", slope_height,
      "num_slopes", #slopes,
      "multiple_bound", multiple_bound,
      "cover_bound", cover_bound,
      "do_rank", BoolString(do_rank),
      "do_twodescent", BoolString(do_twodescent);
print "rhoH", rhoH;
print "sigmaH", sigmaH;
print "tauH", tauH;
print "cH", cH;
print "dH", dH;

for idx in [1..#pairs] do
    if pair_index ne 0 and pair_index ne idx then
        continue;
    end if;

    pair := pairs[idx];
    first := pair[1];
    second := pair[2];
    pair_label := first cat "/" cat second;
    print "PAIR_START", idx, pair_label;

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
    Pbase := C![tangent_m, P0W, 1];
    assert Pbase in C;
    E, phi := EllipticCurve(C, Pbase);

    print "quartic_degree", Degree(f);
    print "squarefree_degree", Degree(sf);
    print "genus", Floor((Degree(sf)-1)/2);
    print "reciprocal_check_a0_eq_a4K2", BoolString(CoeffOrZero(f,0) eq CoeffOrZero(f,4)*K^2);
    print "reciprocal_check_a1_eq_a3K", BoolString(CoeffOrZero(f,1) eq CoeffOrZero(f,3)*K);
    print "hpl_tangent_slope", tangent_m;
    print "elliptic_model", E;

    Emin := E;
    try
        Emin := MinimalModel(E);
        print "minimal_model", Emin;
        print "minimal_a_invariants", aInvariants(Emin);
        print "minimal_discriminant", Discriminant(Emin);
    catch e
        print "minimal_model_error", e`Object;
    end try;

    try
        T, tors_map := TorsionSubgroup(Emin);
        print "torsion_invariants", Invariants(T);
    catch e
        print "torsion_error", e`Object;
    end try;

    try
        print "root_number", RootNumber(Emin);
    catch e
        print "root_number_error", e`Object;
    end try;

    try
        print "analytic_rank", AnalyticRank(Emin);
    catch e
        print "analytic_rank_error", e`Object;
    end try;

    knownC := [
        <"base", Pbase>,
        <"base_neg", C![tangent_m, -P0W, 1]>
    ];
    if K ne 0 then
        Append(~knownC, <"slope_0_plus", C![Q!0, -K*y2, 1]>);
        Append(~knownC, <"slope_0_minus", C![Q!0, K*y2, 1]>);
    end if;

    knownE := [];
    for rec in knownC do
        name := rec[1];
        Pc := rec[2];
        if Pc in C then
            try
                Pe := phi(Pc);
                Append(~knownE, <name, Pe>);
                print "known_point_image",
                      name,
                      "small_order_le_32", SmallTorsionOrder(Pe, 32),
                      "point", Pe;
            catch e
                print "known_point_image_error", name, e`Object;
            end try;
        else
            print "known_point_not_on_curve", name;
        end if;
    end for;

    if do_rank then
        try
            lo, hi := RankBounds(Emin);
            print "rank_bounds", lo, hi;
        catch e
            print "rank_bounds_error", e`Object;
        end try;
    end if;

    descent_points := [];
    if do_twodescent then
        try
            SetVerbose("TwoDescent", 1);
            covers, maps := TwoDescent(E : RemoveTorsion := true);
            print "twodescent_num_covers", #covers;
            for i in [1..#covers] do
                print "TWODESCENT_COVER", i, covers[i];
                try
                    pts := Points(covers[i] : Bound := cover_bound);
                    print "TWODESCENT_POINTS", i, #pts;
                    for j in [1..Min(#pts, 4)] do
                        try
                            Pe := maps[i](pts[j]);
                            Append(~descent_points, <Sprintf("td_%o_%o", i, j), Pe>);
                            print "TWODESCENT_IMAGE", i, j, Pe;
                        catch e
                            print "TWODESCENT_IMAGE_ERROR", i, j, e`Object;
                        end try;
                    end for;
                catch e
                    print "TWODESCENT_POINTS_ERROR", i, e`Object;
                end try;
            end for;
        catch e
            print "twodescent_error", e`Object;
        end try;
    end if;

    stats := rec<StatsFmt |
        seen_sigma := {},
        seen_tuple := {},
        pair_points := 0,
        c_square_unique := 0,
        full_cover_unique := 0,
        simple_certified := 0,
        reports := 0
    >;

    TestQ(~stats, pair_label, qH, "HPL", rhoH, tauH, dH);

    slope_pair_points := 0;
    for s in slopes do
        ok_q, q, y := OtherConicPoint(A1, B1, qH, y1, s);
        if not ok_q then
            continue;
        end if;
        okSecond, y_second := IsNonzeroSquareQ(A2*(1 + B2*q^2));
        if okSecond then
            slope_pair_points +:= 1;
            TestQ(~stats, pair_label, q, "slope_" cat Sprint(s), rhoH, tauH, dH);
        end if;
    end for;

    print "slope_search_pair_points_raw", slope_pair_points;

    multiple_raw := 0;
    inv_ok := false;
    try
        psi := Inverse(phi);
        inv_ok := true;
        tors_pts := TorsionPointsFromSubgroup(E);
        gen_candidates := knownE cat descent_points;
        seenE := {};
        for gen in gen_candidates do
            gen_name := gen[1];
            Pe0 := gen[2];
            ord := SmallTorsionOrder(Pe0, 32);
            if ord ne 0 then
                continue;
            end if;
            for n in [-multiple_bound..multiple_bound] do
                if n eq 0 then
                    continue;
                end if;
                for Tpt in tors_pts do
                    Pe := n*Pe0 + Tpt;
                    ekey := Sprint(Pe);
                    if ekey in seenE then
                        continue;
                    end if;
                    Include(~seenE, ekey);
                    try
                        Pc := psi(Pe);
                        if Pc[3] eq 0 then
                            continue;
                        end if;
                        s := Q!(Pc[1]/Pc[3]);
                        ok_mq, q_m := QFromSlope(A1, B1, qH, y1, s);
                        if not ok_mq then
                            continue;
                        end if;
                        okSecondM, y_second_m := IsNonzeroSquareQ(A2*(1 + B2*q_m^2));
                        if okSecondM then
                            multiple_raw +:= 1;
                            TestQ(~stats, pair_label, q_m,
                                  gen_name cat "_mult_" cat Sprint(n),
                                  rhoH, tauH, dH);
                        end if;
                    catch e
                        continue;
                    end try;
                end for;
            end for;
        end for;
    catch e
        print "inverse_map_error", e`Object;
    end try;

    print "inverse_map_available", BoolString(inv_ok);
    print "multiple_search_pair_points_raw", multiple_raw;
    print "PAIR_DONE", idx, pair_label;
    print "pair_points_unique", stats`pair_points;
    print "c_square_unique", stats`c_square_unique;
    print "full_cover_unique", stats`full_cover_unique;
    print "simple_certified", stats`simple_certified;
end for;

print "DONE";
