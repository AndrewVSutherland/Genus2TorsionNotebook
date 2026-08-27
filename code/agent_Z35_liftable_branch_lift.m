//////////////////////////////////////////////////////////////////////
//  Z/35 lane: liftable b=0,r=1 pole branches beyond the third-pass
//  mod 27 table.
//
//  This continues agent_Z35_b0_pole_blowup.m without editing it.  For
//  the two centers
//
//      (a,b,c0,c1,c2,r) = (1,0,t,t,t,1),  t=1,2  over F_3,
//
//  use
//
//      a=1+3*A, b=3*B, c0=t+3*C0, c1=t+3*C1,
//      c2=t+3*C2, r=1+3*R.
//
//  The previous worker found 9 liftable first directions for each t.
//  Here each such direction is blown up once more.  The branch-level
//  equations are smooth of rank 3 over F_3, so Hensel gives exact lift
//  counts through arbitrary bounded depth.  The script still constructs
//  one canonical lift on every branch and verifies the congruences.
//
//  Typical bounded run:
//      magma -b max_k:=7 recon_height:=80 \
//          code/agent_Z35_liftable_branch_lift.m \
//          > results/Z35_lift_branch_k7.log
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned max_k then
    max_k := 7;
elif Type(max_k) eq MonStgElt then
    max_k := StringToInteger(max_k);
end if;

if not assigned recon_height then
    recon_height := 80;
elif Type(recon_height) eq MonStgElt then
    recon_height := StringToInteger(recon_height);
end if;

if not assigned print_branch_details then
    print_branch_details := 1;
elif Type(print_branch_details) eq MonStgElt then
    print_branch_details := StringToInteger(print_branch_details);
end if;

if not assigned max_store then
    max_store := 2000000;
elif Type(max_store) eq MonStgElt then
    max_store := StringToInteger(max_store);
end if;

Z := Integers();
Q := Rationals();
F3 := GF(3);

R<a,b,c0,c1,c2,r> := PolynomialRing(Z, 6);

// Integral coefficient equations for q^2 - f + (x-r)^5 = 0.
G4 := c2^2 - b^2 - 5*r + 7;
G3 := 2*c1*c2 - 2*a*b - 21 + 10*r^2;
G2 := c1^2 + 2*c0*c2 - a^2 + 7*b + 35 - 10*r^3;
G1 := 2*c0*c1 + 7*a - 2*b - 35 + 5*r^4;
G0 := 4*c0^2 - 8*a + 35 - 4*r^5;
Gs := [G4, G3, G2, G1, G0];
Gnames := ["G4_top", "G3", "G2", "G1", "4*G0"];
VarNames := ["A", "B", "C0", "C1", "C2", "R"];
DigitNames := ["A1", "B1", "C01", "C11", "C21", "R1"];

LiftDirs := AssociativeArray();
LiftDirs[1] := [
    [0,2,0,1,1,0],
    [0,2,1,2,2,1],
    [0,2,2,0,0,2],
    [1,1,0,2,0,2],
    [1,1,1,0,1,0],
    [1,1,2,1,2,1],
    [2,0,0,0,2,1],
    [2,0,1,1,0,2],
    [2,0,2,2,1,0]
];
LiftDirs[2] := [
    [0,2,0,2,2,2],
    [0,2,1,0,0,1],
    [0,2,2,1,1,0],
    [1,1,0,1,0,1],
    [1,1,1,2,1,0],
    [1,1,2,0,2,2],
    [2,0,0,0,1,0],
    [2,0,1,1,2,2],
    [2,0,2,2,0,1]
];

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

function EvalIntMod(poly, vals, m)
    return (Z!Evaluate(poly, vals)) mod m;
end function;

function EvalAllZeroMod(polys, vals, m)
    return &and [ EvalIntMod(poly, vals, m) eq 0 : poly in polys ];
end function;

procedure Increment(~A, key)
    if not IsDefined(A, key) then
        A[key] := 0;
    end if;
    A[key] +:= 1;
end procedure;

function BChartScaledEquations(t)
    S<A,B,C0,C1,C2,RR> := PolynomialRing(Z, 6);
    subs := [1 + 3*A, 3*B, t + 3*C0, t + 3*C1, t + 3*C2, 1 + 3*RR];
    raw := [ Evaluate(g, subs) : g in Gs ];
    scales := [ V3Content(g) : g in raw ];
    scaled := [ DivideByPower(raw[i], scales[i]) : i in [1..#raw] ];
    return S, [A,B,C0,C1,C2,RR], raw, scaled, scales;
end function;

function RankAt(polysF, vals)
    rows := [];
    for poly in polysF do
        Append(~rows, [ Evaluate(Derivative(poly, i), vals) : i in [1..6] ]);
    end for;
    return Rank(Matrix(Parent(vals[1]), rows));
end function;

function F3Solutions(polys)
    P := Parent(polys[1]);
    PF := ChangeRing(P, F3);
    HF := [ PF!h : h in polys ];
    sols := [];
    for A0 in F3 do
    for B0 in F3 do
    for C00 in F3 do
    for C10 in F3 do
    for C20 in F3 do
    for R0 in F3 do
        vals := [A0,B0,C00,C10,C20,R0];
        if &and [ Evaluate(h, vals) eq 0 : h in HF ] then
            Append(~sols, vals);
        end if;
    end for;
    end for;
    end for;
    end for;
    end for;
    end for;
    return sols, HF;
end function;

function NormalizeEquation(coeffs, c)
    idx := 0;
    for i in [1..#coeffs] do
        if coeffs[i] ne 0 then
            idx := i;
            break;
        end if;
    end for;
    if idx eq 0 then
        return false, coeffs, c;
    end if;
    s := coeffs[idx]^-1;
    return true, [ s*coeffs[i] : i in [1..#coeffs] ], s*c;
end function;

function EquationKey(coeffs, c)
    return Sprint(< [ Z!coeffs[i] : i in [1..#coeffs] ], Z!c >);
end function;

function AffineLinearBasis(points)
    if #points eq 0 then
        return [];
    end if;
    F := Parent(points[1][1]);
    all_eqs := [];
    seen := AssociativeArray();
    for u1 in F do
    for u2 in F do
    for u3 in F do
    for u4 in F do
    for u5 in F do
    for u6 in F do
        coeffs := [u1,u2,u3,u4,u5,u6];
        if &and [ coeffs[i] eq 0 : i in [1..6] ] then
            continue;
        end if;
        for cc in F do
            ok := true;
            for pt in points do
                if &+[ coeffs[i]*pt[i] : i in [1..6] ] + cc ne 0 then
                    ok := false;
                    break;
                end if;
            end for;
            if ok then
                good, nc, nd := NormalizeEquation(coeffs, cc);
                key := EquationKey(nc, nd);
                if good and not IsDefined(seen, key) then
                    seen[key] := true;
                    Append(~all_eqs, <nc, nd, key>);
                end if;
            end if;
        end for;
    end for;
    end for;
    end for;
    end for;
    end for;
    end for;
    Sort(~all_eqs, func<x,y | x[3] lt y[3] select -1 else x[3] gt y[3] select 1 else 0>);

    basis := [];
    for rel in all_eqs do
        old_rank := (#basis eq 0) select 0 else Rank(Matrix(F, [ b[1] : b in basis ]));
        new_rows := [ b[1] : b in basis ];
        Append(~new_rows, rel[1]);
        new_rank := Rank(Matrix(F, new_rows));
        if new_rank gt old_rank then
            Append(~basis, rel);
        end if;
    end for;
    return basis;
end function;

function LinearEquationString(rel, names)
    coeffs := rel[1];
    c := rel[2];
    terms := [];
    for i in [1..#coeffs] do
        ci := Z!coeffs[i];
        if ci eq 1 then
            Append(~terms, names[i]);
        elif ci eq 2 then
            Append(~terms, "2*" cat names[i]);
        end if;
    end for;
    if Z!c eq 1 then
        Append(~terms, "1");
    elif Z!c eq 2 then
        Append(~terms, "2");
    end if;
    if #terms eq 0 then
        return "0 = 0";
    end if;
    return Join(terms, " + ") cat " = 0";
end function;

function BranchEquations(t, dir)
    S, vars, raw, Hs, scales := BChartScaledEquations(t);
    T<A1,B1,C01,C11,C21,R1> := PolynomialRing(Z, 6);
    subs := [ dir[1] + 3*A1, dir[2] + 3*B1, dir[3] + 3*C01,
              dir[4] + 3*C11, dir[5] + 3*C21, dir[6] + 3*R1 ];
    phi := hom<S -> T | subs>;
    raw2 := [ phi(h) : h in Hs ];
    scales2 := [ V3Content(h) : h in raw2 ];
    Ks := [ DivideByPower(raw2[i], 1) : i in [1..#raw2] ];
    return T, Ks, scales2;
end function;

function VecKey(vals)
    return Sprint([ Z!vals[i] : i in [1..#vals] ]);
end function;

function CorrectionTable(polys, residue)
    // Precompute J(residue)*d for all d in F_3^6.  Along a fixed
    // first-direction branch, residue mod 3 is constant, so this avoids
    // re-enumerating 729 corrections for every lift.
    P := Parent(polys[1]);
    PF := ChangeRing(P, F3);
    HF := [ PF!h : h in polys ];
    valsF := [ F3!(residue[i] mod 3) : i in [1..6] ];
    Jrows := [];
    for h in HF do
        Append(~Jrows, [ Evaluate(Derivative(h, j), valsF) : j in [1..6] ]);
    end for;
    table := AssociativeArray();
    for u1 in F3 do
    for u2 in F3 do
    for u3 in F3 do
    for u4 in F3 do
    for u5 in F3 do
    for u6 in F3 do
        d := [u1,u2,u3,u4,u5,u6];
        lhs := [ &+[ Jrows[i][j]*d[j] : j in [1..6] ] : i in [1..#Jrows] ];
        key := VecKey(lhs);
        if not IsDefined(table, key) then
            table[key] := [];
        end if;
        Append(~table[key], [ Z!d[i] : i in [1..6] ]);
    end for;
    end for;
    end for;
    end for;
    end for;
    end for;
    return table;
end function;

function LiftCorrectionsFromTable(polys, vals, base_mod, table)
    rhs := [ F3!(-ExactQuotient(Z!Evaluate(polys[i], vals), base_mod)) : i in [1..#polys] ];
    key := VecKey(rhs);
    if IsDefined(table, key) then
        return table[key];
    end if;
    return [];
end function;

function LiftBranchExact(polys, start_vals, target_k, store_limit)
    residues := [ [ Z!x : x in start_vals ] ];
    base_mod := 3;
    table := CorrectionTable(polys, start_vals);
    rows := [];
    sample_depth := 1;
    sample := residues[1];
    stopped_reason := "completed";

    for kk in [2..target_k] do
        next_residues := [];
        hist := AssociativeArray();
        total := 0;
        can_store := true;
        do_lookahead := kk + 1 eq target_k;
        look_total := 0;
        look_hist := AssociativeArray();

        for vals in residues do
            corrs := LiftCorrectionsFromTable(polys, vals, base_mod, table);
            Increment(~hist, Sprint(#corrs));
            total +:= #corrs;
            if #corrs ne 0 and kk lt target_k then
                for d in corrs do
                    new_val := [ vals[i] + base_mod*d[i] : i in [1..6] ];
                    if do_lookahead then
                        next_corrs := LiftCorrectionsFromTable(polys, new_val, 3*base_mod, table);
                        Increment(~look_hist, Sprint(#next_corrs));
                        look_total +:= #next_corrs;
                    end if;
                    if can_store then
                        if #next_residues lt store_limit then
                            Append(~next_residues, new_val);
                        else
                            can_store := false;
                            next_residues := [];
                        end if;
                    end if;
                end for;
                if can_store and #next_residues gt store_limit then
                    can_store := false;
                    next_residues := [];
                end if;
            end if;
        end for;

        Append(~rows, <kk, 3^kk, total,
                       Sort([ <StringToInteger(k), hist[k]> : k in Keys(hist) ]),
                       can_store>);

        if kk eq target_k then
            if total eq 0 then
                stopped_reason := "no_lifts";
            end if;
            break;
        end if;

        if total eq 0 then
            stopped_reason := "no_lifts";
            break;
        end if;
        if not can_store then
            if do_lookahead then
                Append(~rows, <kk + 1, 3^(kk + 1), look_total,
                               Sort([ <StringToInteger(k), look_hist[k]> : k in Keys(look_hist) ]),
                               false>);
                stopped_reason := "store_limit_after_lookahead";
            else
                stopped_reason := "store_limit";
            end if;
            break;
        end if;

        residues := next_residues;
        base_mod *:= 3;
        sample_depth := kk;
        sample := residues[1];
    end for;

    return rows, stopped_reason, sample_depth, sample;
end function;

function OriginalVals(t, local_vals)
    return [
        1 + 3*local_vals[1],
        3*local_vals[2],
        t + 3*local_vals[3],
        t + 3*local_vals[4],
        t + 3*local_vals[5],
        1 + 3*local_vals[6]
    ];
end function;

function Contact7Data(aa, bb)
    P<x> := PolynomialRing(Q);
    h := 1 - (Q!7/Q!2)*x + (Q!aa)*x^2 + (Q!bb)*x^3;
    f := ExactQuotient(h^2 + (x - 1)^7, x^2);
    return f, h;
end function;

function SmallRationalResidue(n, m, H)
    for den in [1..H] do
        if GCD(den, 3) ne 1 then
            continue;
        end if;
        for num in [-H..H] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            if ((num - n*den) mod m) eq 0 then
                return true, Q!num/Q!den;
            end if;
        end for;
    end for;
    return false, Q!0;
end function;

function ReconstructTuple(vals, m, H)
    rats := [];
    for v in vals do
        ok, rr := SmallRationalResidue(v mod m, m, H);
        if not ok then
            return false, rats;
        end if;
        Append(~rats, rr);
    end for;
    return true, rats;
end function;

function EvalOriginalQ(vals)
    return [ Evaluate(ChangeRing(g, Q), vals) : g in Gs ];
end function;

function DirectionPlaneFormulaOk(t, v)
    A0 := v[1]; B0 := v[2]; C00 := v[3]; C10 := v[4]; C20 := v[5]; R0 := v[6];
    if t eq 1 then
        return (B0 - (2 - A0)) mod 3 eq 0
           and (R0 - (C00 - A0)) mod 3 eq 0
           and (C10 - (C00 + A0 + 1)) mod 3 eq 0
           and (C20 - (C00 + 1 - A0)) mod 3 eq 0;
    else
        return (B0 - (2 - A0)) mod 3 eq 0
           and (R0 - (2 - C00 - A0)) mod 3 eq 0
           and (C10 - (C00 + 2 - A0)) mod 3 eq 0
           and (C20 - (C00 + A0 + 2)) mod 3 eq 0;
    end if;
end function;

print "Z35 liftable b=0,r=1 branch lift";
print "max_k", max_k, "recon_height", recon_height, "max_store", max_store;
print "local equation convention: H(A,B,C0,C1,C2,R)=0 mod 3^k";
print "original coefficient equations then vanish mod 3^(k+1)";
print "variables", VarNames;

print "FIRST_DIRECTION_PLANES";
for t in [1,2] do
    pts := [ [ F3!v[i] : i in [1..6] ] : v in LiftDirs[t] ];
    basis := AffineLinearBasis(pts);
    print "center_t", t, "directions", #LiftDirs[t], "plane_basis";
    for rel in basis do
        print " ", LinearEquationString(rel, VarNames);
    end for;
    print "formula_check", &and [ DirectionPlaneFormulaOk(t, v) : v in LiftDirs[t] ];
    if t eq 1 then
        print "formula", "B=2-A, R=C0-A, C1=C0+A+1, C2=C0+1-A over F3";
    else
        print "formula", "B=2-A, R=2-C0-A, C1=C0+2-A, C2=C0+A+2 over F3";
    end if;
end for;

print "NAIVE_SMOOTH_MAX_TABLE";
print "columns: k modulus max_per_direction_if_unobstructed";
print "This is only a reference bound after the third-pass table, not a verdict.";
for k in [1..max_k] do
    per_dir := (k eq 1) select 1 else 3^(3*(k - 1));
    print k, 3^k, per_dir;
end for;

print "BRANCH_ANALYSIS";
print "columns include: branch_scales_for_H(dir+3*Y), mod3_solutions, rank_histogram, affine_basis";

TotalByK := AssociativeArray();
BranchesByK := AssociativeArray();
ExactBranchesByK := AssociativeArray();
for k in [1..max_k] do
    TotalByK[Sprint(k)] := 0;
    BranchesByK[Sprint(k)] := 0;
    ExactBranchesByK[Sprint(k)] := 0;
end for;
StopCounts := AssociativeArray();

for t in [1,2] do
    S, vars, raw, Hs, scales := BChartScaledEquations(t);
    print "CENTER", t, "chart_scales", scales;
    all_canonical := [];
    for idx in [1..#LiftDirs[t]] do
        dir := LiftDirs[t][idx];
        T, Ks, scales2 := BranchEquations(t, dir);
        sols, KF := F3Solutions(Ks);
        rank_counts := AssociativeArray();
        for pt in sols do
            rk := RankAt(KF, pt);
            Increment(~rank_counts, Sprint(rk));
        end for;
        basis := AffineLinearBasis(sols);

        lift_rows, stopped_reason, sample_depth, sample_vals :=
            LiftBranchExact(Hs, dir, max_k, max_store);
        Increment(~StopCounts, stopped_reason);

        branch_counts := AssociativeArray();
        branch_counts["1"] := 1;
        branch_exact := AssociativeArray();
        branch_exact["1"] := true;
        for row in lift_rows do
            branch_counts[Sprint(row[1])] := row[3];
            branch_exact[Sprint(row[1])] := true;
        end for;
        if #lift_rows gt 0 and stopped_reason eq "no_lifts" then
            last_k := lift_rows[#lift_rows][1];
            if last_k lt max_k then
                for kk in [last_k + 1..max_k] do
                    branch_counts[Sprint(kk)] := 0;
                    branch_exact[Sprint(kk)] := true;
                end for;
            end if;
        end if;

        for kk in [1..max_k] do
            key := Sprint(kk);
            if IsDefined(branch_counts, key) then
                TotalByK[key] +:= branch_counts[key];
                if branch_counts[key] gt 0 then
                    BranchesByK[key] +:= 1;
                end if;
                if IsDefined(branch_exact, key) then
                    ExactBranchesByK[key] +:= 1;
                end if;
            end if;
        end for;
        Append(~all_canonical, sample_vals);

        if print_branch_details ne 0 then
            print "BRANCH", "t", t, "idx", idx, "dir", dir;
            print " branch_scales", scales2;
            print " mod3_solution_count", #sols, "rank_histogram",
                  Sort([ <StringToInteger(k), rank_counts[k]> : k in Keys(rank_counts) ]);
            print " next_digit_affine_basis";
            for rel in basis do
                print "  ", LinearEquationString(rel, DigitNames);
            end for;
            print " exact_lift_rows", lift_rows;
            print " stopped_reason", stopped_reason, "sample_depth", sample_depth;
            print " sample_local_mod", 3^sample_depth, sample_vals;
            orig := OriginalVals(t, sample_vals);
            print " sample_original_mod", 3^(sample_depth + 1), orig;
            print " original_G_residues_mod_3^(sample_depth+1)",
                  [ EvalIntMod(Gs[i], orig, 3^(sample_depth + 1)) : i in [1..#Gs] ];
            f, h := Contact7Data(orig[1], orig[2]);
            disc := Discriminant(f);
            print " sample_h1", Evaluate(h, Q!1), "sample_disc_v3", Valuation(Numerator(disc), 3) - Valuation(Denominator(disc), 3),
                  "sample_disc_nonzero", disc ne 0;
            recon_ok, recon := ReconstructTuple(sample_vals, 3^sample_depth, recon_height);
            print " small_rational_reconstruction_local_height", recon_height, "ok", recon_ok;
            if recon_ok then
                origQ := [ Q!1 + 3*recon[1], 3*recon[2], Q!t + 3*recon[3],
                           Q!t + 3*recon[4], Q!t + 3*recon[5], Q!1 + 3*recon[6] ];
                evalsQ := EvalOriginalQ(origQ);
                print " reconstructed_local", recon;
                print " reconstructed_original", origQ;
                print " reconstructed_exact_zero", &and [ e eq 0 : e in evalsQ ];
            end if;
        end if;
    end for;
end for;

print "SUMMARY";
print "stop_counts", Sort([ <k, StopCounts[k]> : k in Keys(StopCounts) ]);
print "EXACT_TOTAL_LIFT_TABLE";
print "columns: k modulus total_lifts_all_18 nonzero_branches exact_branches";
for k in [1..max_k] do
    key := Sprint(k);
    print k, 3^k, TotalByK[key], BranchesByK[key], ExactBranchesByK[key];
end for;
if max_k ge 4 and BranchesByK["4"] eq 0 then
    print "lift_verdict", "all 18 branches die at the mod 81 step; there are no lifts beyond the third-pass mod 27 table";
else
    print "lift_verdict", "some branches survive past mod 27; inspect exact branch rows and store-limit status";
end if;
print "contact5_verdict", "the five point-contact coefficient equations are the lifted equations; stored samples verify G_i=0 at their recorded precision";
print "rational_reconstruction_verdict", "small reconstruction probe is reported branch-by-branch and is only a heuristic, not a proof of nonexistence";

quit;
