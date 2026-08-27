//////////////////////////////////////////////////////////////////////
//  Z/35 lane: b=0 pole blow-up for the simultaneous
//  contact-7 / point-contact-5 equations.
//
//  This returns to the original coefficient equations in
//      (a,b,c0,c1,c2,r)
//  for
//      q^2 - f = -(x-r)^5,     q=c0+c1*x+c2*x^2,
//      f=(h^2+(x-1)^7)/x^2,   h=1-(7/2)x+a*x^2+b*x^3.
//
//  It deliberately does not use the old d,e formulas that divide by b
//  or c2.  The main local charts are the two b=0,r=1 centers
//
//      (a,b,c0,c1,c2,r) = (1,0,t,t,t,1),  t=1,2  over F_3,
//
//  with
//
//      a=1+3*A, b=3*B, c0=t+3*C0, c1=t+3*C1,
//      c2=t+3*C2, r=1+3*R.
//
//  Optional: also inspect the finite c2=0 chart.  The original
//  equations already show that the residual d=-e center has no finite
//  F_3 lift in (a,b,c0,c1,c2,r).
//
//  Typical run:
//      magma code/agent_Z35_b0_pole_blowup.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned do_c2_chart then
    do_c2_chart := 1;
elif Type(do_c2_chart) eq MonStgElt then
    do_c2_chart := StringToInteger(do_c2_chart);
end if;

if not assigned do_good_prime_table then
    do_good_prime_table := 1;
elif Type(do_good_prime_table) eq MonStgElt then
    do_good_prime_table := StringToInteger(do_good_prime_table);
end if;

Z := Integers();
Q := Rationals();

R<a,b,c0,c1,c2,r> := PolynomialRing(Z, 6);

// Integral coefficient equations for q^2 - f + (x-r)^5 = 0.
G4 := c2^2 - b^2 - 5*r + 7;
G3 := 2*c1*c2 - 2*a*b - 21 + 10*r^2;
G2 := c1^2 + 2*c0*c2 - a^2 + 7*b + 35 - 10*r^3;
G1 := 2*c0*c1 + 7*a - 2*b - 35 + 5*r^4;
G0 := 4*c0^2 - 8*a + 35 - 4*r^5;
Gs := [G4, G3, G2, G1, G0];
Gnames := ["G4_top", "G3", "G2", "G1", "4*G0"];

function V3Content(poly)
    coeffs := [ Z!c : c in Coefficients(poly) | c ne 0 ];
    if #coeffs eq 0 then
        return 999;
    end if;
    return Minimum([ Valuation(c, 3) : c in coeffs ]);
end function;

function DivideByPower(poly, e)
    if e eq 0 then
        return poly;
    end if;
    P := Parent(poly);
    mons := Monomials(poly);
    coeffs := Coefficients(poly);
    if #mons eq 0 then
        return P!0;
    end if;
    return &+[ P!(ExactQuotient(Z!coeffs[i], 3^e))*mons[i] : i in [1..#mons] ];
end function;

function TupleString(vals)
    return Sprint(< vals[i] : i in [1..#vals] >);
end function;

procedure Increment(~A, key)
    if not IsDefined(A, key) then
        A[key] := 0;
    end if;
    A[key] +:= 1;
end procedure;

function Contact7F(F, aa, bb)
    P<x> := PolynomialRing(F);
    h := 1 - (F!7/F!2)*x + aa*x^2 + bb*x^3;
    f := ExactQuotient(h^2 + (x - 1)^7, x^2);
    return f, h;
end function;

function EvalIntMod(poly, vals, m)
    return (Z!Evaluate(poly, vals)) mod m;
end function;

function EvalAllZeroMod(polys, vals, m)
    return &and [ EvalIntMod(poly, vals, m) eq 0 : poly in polys ];
end function;

function LiftList(polys, base_vals, base_mod, next_mod)
    lifts := [];
    for u1 in [0..2] do
    for u2 in [0..2] do
    for u3 in [0..2] do
    for u4 in [0..2] do
    for u5 in [0..2] do
    for u6 in [0..2] do
        corr := [u1,u2,u3,u4,u5,u6];
        vals := [ base_vals[i] + base_mod*corr[i] : i in [1..6] ];
        if EvalAllZeroMod(polys, vals, next_mod) then
            Append(~lifts, vals);
        end if;
    end for;
    end for;
    end for;
    end for;
    end for;
    end for;
    return lifts;
end function;

function RankAt(polysF, varsF, vals)
    F := Parent(vals[1]);
    rows := [];
    for poly in polysF do
        Append(~rows, [ Evaluate(Derivative(poly, i), vals) : i in [1..#varsF] ]);
    end for;
    return Rank(Matrix(F, rows));
end function;

function BChartScaledEquations(t)
    S<A,B,C0,C1,C2,RR> := PolynomialRing(Z, 6);
    subs := [1 + 3*A, 3*B, t + 3*C0, t + 3*C1, t + 3*C2, 1 + 3*RR];
    raw := [ Evaluate(g, subs) : g in Gs ];
    scales := [ V3Content(g) : g in raw ];
    scaled := [ DivideByPower(raw[i], scales[i]) : i in [1..#raw] ];
    return S, [A,B,C0,C1,C2,RR], raw, scaled, scales;
end function;

procedure AnalyzeBChart(t)
    S, vars, raw, Hs, scales := BChartScaledEquations(t);
    F := GF(3);
    SF := ChangeRing(S, F);
    HF := [ SF!h : h in Hs ];
    varsF := [ SF.i : i in [1..6] ];

    print "B0_CHART";
    print "center_t", t, "center", <1,0,t,t,t,1>;
    print "variables", "A B C0 C1 C2 R";
    print "scales_v3", [ <Gnames[i], scales[i], Degree(Hs[i]), #Terms(Hs[i])> : i in [1..#Hs] ];
    print "scaled_mod3_factorizations";
    for i in [1..#Hs] do
        print " ", Gnames[i], Factorization(HF[i]);
    end for;

    sols := [];
    for A0 in F do
    for B0 in F do
    for D0 in F do
    for D1 in F do
    for D2 in F do
    for R0 in F do
        valsF := [A0,B0,D0,D1,D2,R0];
        if &and [ Evaluate(h, valsF) eq 0 : h in HF ] then
            Append(~sols, valsF);
        end if;
    end for;
    end for;
    end for;
    end for;
    end for;
    end for;

    rank_counts := AssociativeArray();
    lift9_counts := AssociativeArray();
    lift27_counts := AssociativeArray();
    rows := [];
    total9 := 0;
    total27 := 0;
    liftable9 := 0;
    liftable27 := 0;

    for valsF in sols do
        vals := [ Z!v : v in valsF ];
        rk := RankAt(HF, varsF, valsF);
        Increment(~rank_counts, Sprint(rk));

        lifts9 := LiftList(Hs, vals, 3, 9);
        n9 := #lifts9;
        total9 +:= n9;
        if n9 gt 0 then
            liftable9 +:= 1;
        end if;
        Increment(~lift9_counts, Sprint(n9));

        n27 := 0;
        for lift in lifts9 do
            n27 +:= #LiftList(Hs, lift, 9, 27);
        end for;
        total27 +:= n27;
        if n27 gt 0 then
            liftable27 +:= 1;
        end if;
        Increment(~lift27_counts, Sprint(n27));

        Append(~rows, <vals, rk, 6-rk, n9, n27>);
    end for;

    print "F3_DIRECTION_SUMMARY";
    print "directions", #sols, "rank_counts", Sort([ <StringToInteger(k), rank_counts[k]> : k in Keys(rank_counts) ]);
    print "mod9_total_lifts", total9, "directions_liftable_mod9", liftable9;
    print "mod27_total_lifts", total27, "directions_liftable_mod27", liftable27;
    print "mod9_lift_count_histogram", Sort([ <StringToInteger(k), lift9_counts[k]> : k in Keys(lift9_counts) ]);
    print "mod27_lift_count_histogram", Sort([ <StringToInteger(k), lift27_counts[k]> : k in Keys(lift27_counts) ]);
    print "F3_DIRECTION_ROWS";
    print "columns: <A,B,C0,C1,C2,R> jac_rank tangent_dim mod9_lifts mod27_lifts";
    for row in rows do
        print row[1], row[2], row[3], row[4], row[5];
    end for;
end procedure;

procedure AnalyzeFiniteC2Center()
    F := GF(3);
    RF := ChangeRing(R, F);
    Gf := [ RF!g : g in Gs ];
    // d=-e residual center has b=1, c1=c2=r=0, while c0 was a pole
    // in the old chart.  Test all finite a,c0 in the original equations.
    sols := [];
    for aa in F do
        for cc0 in F do
            vals := [aa, F!1, cc0, F!0, F!0, F!0];
            if &and [ Evaluate(g, vals) eq 0 : g in Gf ] then
                Append(~sols, vals);
            end if;
        end for;
    end for;
    print "C2_0_FINITE_CENTER_CHECK";
    print "residual_center_d_e_c1", <2,1,0>;
    print "forced_original_residues", "b=1,c1=0,c2=0,r=0";
    print "finite_original_F3_lifts", #sols, sols;
    if #sols eq 0 then
        print "verdict", "no finite original-variable F3 center";
        print "polar_note", "with a,r integral at this center, 4*G0 contains 4*c0^2 plus integral terms, so a simple v3(c0)<0 escape cannot cancel";
    end if;
end procedure;

function OriginalEquationsHold(F, aa, bb, cc0, cc1, cc2, rr)
    vals := [aa, bb, cc0, cc1, cc2, rr];
    RF := ChangeRing(R, F);
    Gf := [ RF!g : g in Gs ];
    return &and [ Evaluate(g, vals) eq 0 : g in Gf ];
end function;

procedure GoodPrimePointContactTable(primes)
    print "GOOD_PRIME_POINT_CONTACT_TABLE";
    print "columns: p total smooth pass5 b0_points c20_points";
    for p in primes do
        F := GF(p);
        if p in {2,5} then
            continue;
        end if;
        total := 0;
        smooth := 0;
        pass5 := 0;
        b0pts := 0;
        c20pts := 0;
        RF := ChangeRing(R, F);
        Gf := [ RF!g : g in Gs ];

        // Use E4 to solve r.  Use the linear equations in a and c0
        // only when the relevant coefficient is a unit; otherwise
        // enumerate that variable.  The original equations are still
        // checked at the end, including the b=0 and c2=0 poles.
        for bb in F do
        for cc2 in F do
        for cc1 in F do
            rr := (cc2^2 - bb^2 + 7)/5;
            avecs := [];
            if bb ne 0 then
                aa := (2*cc1*cc2 - 21 + 10*rr^2)/(2*bb);
                Append(~avecs, aa);
            else
                for aa in F do
                    if Evaluate(Gf[2], [aa,bb,F!0,cc1,cc2,rr]) eq 0 then
                        Append(~avecs, aa);
                    end if;
                end for;
            end if;

            for aa in avecs do
                c0vecs := [];
                if cc2 ne 0 then
                    cc0 := (aa^2 - 7*bb - 35 + 10*rr^3 - cc1^2)/(2*cc2);
                    Append(~c0vecs, cc0);
                else
                    for cc0 in F do
                        if Evaluate(Gf[3], [aa,bb,cc0,cc1,cc2,rr]) eq 0 then
                            Append(~c0vecs, cc0);
                        end if;
                    end for;
                end if;

                for cc0 in c0vecs do
                    vals := [aa,bb,cc0,cc1,cc2,rr];
                    if not &and [ Evaluate(g, vals) eq 0 : g in Gf ] then
                        continue;
                    end if;
                    total +:= 1;
                    if bb eq 0 then b0pts +:= 1; end if;
                    if cc2 eq 0 then c20pts +:= 1; end if;

                    f, h := Contact7F(F, aa, bb);
                    if Degree(f) eq 5 and Discriminant(f) ne 0 and Evaluate(h, F!1) ne 0 then
                        smooth +:= 1;
                        n := Z!#Jacobian(HyperellipticCurve(f));
                        if n mod 5 eq 0 then
                            pass5 +:= 1;
                        end if;
                    end if;
                end for;
            end for;
        end for;
        end for;
        end for;

        print p, total, smooth, pass5, b0pts, c20pts;
    end for;
end procedure;

print "Z35 b=0 pole blow-up: original simultaneous contact equations";
print "original_equations", Gnames;
print "b=0 centers over F3", <1,0,1,1,1,1>, <1,0,2,2,2,1>;

AnalyzeBChart(1);
AnalyzeBChart(2);

if do_c2_chart ne 0 then
    AnalyzeFiniteC2Center();
end if;

if do_good_prime_table ne 0 then
    GoodPrimePointContactTable([7,11,13,17,19,23,29,31]);
end if;

quit;
