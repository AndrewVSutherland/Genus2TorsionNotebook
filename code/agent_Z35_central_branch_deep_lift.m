//////////////////////////////////////////////////////////////////////
//  Z/35 lane: deep lift of the two surviving central b=0,r=1
//  pole directions.
//
//  This is a disjoint continuation of agent_Z35_liftable_branch_lift.m.
//  It only studies
//
//      t=1, dir=<1,1,1,0,1,0>,
//      t=2, dir=<1,1,1,2,1,0>.
//
//  Convention: after
//
//      a=1+3*A, b=3*B, c0=t+3*C0, c1=t+3*C1,
//      c2=t+3*C2, r=1+3*R,
//
//  the five transformed equations are divided by one power of 3.
//  Thus H=0 mod 3^k implies the original coefficient equations vanish
//  mod 3^(k+1).
//
//  Typical run:
//
//      magma -b max_k:=10 direct_depth:=6 \
//          recon_height:=1000 max_store:=2500000 \
//          code/agent_Z35_central_branch_deep_lift.m \
//          > results/Z35_central_deep_k10.log
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned max_k then
    max_k := 10;
elif Type(max_k) eq MonStgElt then
    max_k := StringToInteger(max_k);
end if;

if not assigned direct_depth then
    direct_depth := 6;
elif Type(direct_depth) eq MonStgElt then
    direct_depth := StringToInteger(direct_depth);
end if;

if not assigned max_store then
    max_store := 2500000;
elif Type(max_store) eq MonStgElt then
    max_store := StringToInteger(max_store);
end if;

if not assigned recon_height then
    recon_height := 1000;
elif Type(recon_height) eq MonStgElt then
    recon_height := StringToInteger(recon_height);
end if;

if not assigned use_compressed then
    use_compressed := 0;
elif Type(use_compressed) eq MonStgElt then
    use_compressed := StringToInteger(use_compressed);
end if;

if not assigned sample_parent_limit then
    sample_parent_limit := 2000;
elif Type(sample_parent_limit) eq MonStgElt then
    sample_parent_limit := StringToInteger(sample_parent_limit);
end if;

Z := Integers();
Q := Rationals();
F3 := GF(3);

R0<a,b,c0,c1,c2,r> := PolynomialRing(Z, 6);

G4 := c2^2 - b^2 - 5*r + 7;
G3 := 2*c1*c2 - 2*a*b - 21 + 10*r^2;
G2 := c1^2 + 2*c0*c2 - a^2 + 7*b + 35 - 10*r^3;
G1 := 2*c0*c1 + 7*a - 2*b - 35 + 5*r^4;
G0 := 4*c0^2 - 8*a + 35 - 4*r^5;
Gs := [G4, G3, G2, G1, G0];

VarNames := ["A", "B", "C0", "C1", "C2", "R"];
DigitNames := ["dA", "dB", "dC0", "dC1", "dC2", "dR"];

CentralDirs := AssociativeArray();
CentralDirs[1] := [1, 1, 1, 0, 1, 0];
CentralDirs[2] := [1, 1, 1, 2, 1, 0];

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
    return &+[ P!(ExactQuotient(Z!coeffs[i], 3^e))*mons[i] :
               i in [1..#mons] ];
end function;

procedure Increment(~A, key)
    if not IsDefined(A, key) then
        A[key] := 0;
    end if;
    A[key] +:= 1;
end procedure;

procedure IncrementBy(~A, key, wt)
    if not IsDefined(A, key) then
        A[key] := 0;
    end if;
    A[key] +:= wt;
end procedure;

function SortedHist(A)
    return Sort([ <k, A[k]> : k in Keys(A) ]);
end function;

function VecKey(vals)
    return Sprint([ Z!vals[i] : i in [1..#vals] ]);
end function;

function BChartScaledEquations(t)
    S<A,B,C0,C1,C2,RR> := PolynomialRing(Z, 6);
    subs := [1 + 3*A, 3*B, t + 3*C0, t + 3*C1, t + 3*C2,
             1 + 3*RR];
    raw := [ Evaluate(g, subs) : g in Gs ];
    scales := [ V3Content(g) : g in raw ];
    scaled := [ DivideByPower(raw[i], scales[i]) : i in [1..#raw] ];
    derivs := [ [ Derivative(scaled[i], j) : j in [1..6] ] :
                i in [1..#scaled] ];
    return S, scaled, derivs, scales;
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
        vals := [A0, B0, C00, C10, C20, R0];
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
        coeffs := [u1, u2, u3, u4, u5, u6];
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

function BranchEquations(Hs, dir)
    S := Parent(Hs[1]);
    T<dA,dB,dC0,dC1,dC2,dR> := PolynomialRing(Z, 6);
    subs := [ dir[1] + 3*dA, dir[2] + 3*dB, dir[3] + 3*dC0,
              dir[4] + 3*dC1, dir[5] + 3*dC2, dir[6] + 3*dR ];
    phi := hom<S -> T | subs>;
    raw2 := [ phi(h) : h in Hs ];
    Ks := [ DivideByPower(raw2[i], 1) : i in [1..#raw2] ];
    return T, Ks;
end function;

function JacobianMod3(Hs, dir)
    P := Parent(Hs[1]);
    PF := ChangeRing(P, F3);
    HF := [ PF!h : h in Hs ];
    vals := [ F3!(dir[i] mod 3) : i in [1..6] ];
    return Matrix(F3, [ [ Evaluate(Derivative(h, j), vals) : j in [1..6] ] :
                        h in HF ]);
end function;

function IntegerRowsFromBasis(B)
    return [ [ Z!(v[i]) : i in [1..Degree(Parent(v))] ] : v in B ];
end function;

function CorrectionTable(polys, residue)
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
        d := [u1, u2, u3, u4, u5, u6];
        lhs := [ &+[ Jrows[i][j]*d[j] : j in [1..6] ] :
                 i in [1..#Jrows] ];
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
    rhs := [ F3!(-ExactQuotient(Z!Evaluate(polys[i], vals), base_mod)) :
             i in [1..#polys] ];
    key := VecKey(rhs);
    if IsDefined(table, key) then
        return table[key];
    end if;
    return [];
end function;

function LiftStoreToDepth(polys, start_vals, target_k, store_limit)
    residues := [ [ Z!x : x in start_vals ] ];
    base_mod := 3;
    table := CorrectionTable(polys, start_vals);
    rows := [ <1, 3, 1, [ <27, 0> ], true> ];
    stopped_reason := "completed";

    for kk in [2..target_k] do
        next_residues := [];
        hist := AssociativeArray();
        total := 0;
        can_store := true;
        for vals in residues do
            corrs := LiftCorrectionsFromTable(polys, vals, base_mod, table);
            Increment(~hist, Sprint(#corrs));
            total +:= #corrs;
            for d in corrs do
                if can_store then
                    if #next_residues lt store_limit then
                        Append(~next_residues, [ vals[i] + base_mod*d[i] : i in [1..6] ]);
                    else
                        can_store := false;
                        next_residues := [];
                    end if;
                end if;
            end for;
        end for;
        Append(~rows, <kk, 3^kk, total,
                       Sort([ <StringToInteger(k), hist[k]> : k in Keys(hist) ]),
                       can_store>);
        if total eq 0 then
            stopped_reason := "no_lifts";
            return rows, stopped_reason, next_residues;
        end if;
        if not can_store then
            stopped_reason := "store_limit";
            return rows, stopped_reason, next_residues;
        end if;
        residues := next_residues;
        base_mod *:= 3;
    end for;
    return rows, stopped_reason, residues;
end function;

function LookaheadOneRow(polys, residues, current_k, start_vals)
    table := CorrectionTable(polys, start_vals);
    base_mod := 3^current_k;
    hist := AssociativeArray();
    total := 0;
    for vals in residues do
        corrs := LiftCorrectionsFromTable(polys, vals, base_mod, table);
        Increment(~hist, Sprint(#corrs));
        total +:= #corrs;
    end for;
    return <current_k + 1, 3^(current_k + 1), total,
            Sort([ <StringToInteger(k), hist[k]> : k in Keys(hist) ]),
            false>;
end function;

function ProjectRowsByObservedParity(known_k, known_count, target_k)
    rows := [];
    count := known_count;
    for kk in [known_k + 1..target_k] do
        multiplier := (kk mod 2 eq 0) select 9 else 27;
        count *:= multiplier;
        Append(~rows, <kk, 3^kk, count, multiplier,
                       (kk mod 2 eq 0) select "obstructed_fraction_1_over_3" else "full_27">);
    end for;
    return rows;
end function;

function AllDigitCorrections()
    digs := [];
    for u1 in [0..2] do
    for u2 in [0..2] do
    for u3 in [0..2] do
    for u4 in [0..2] do
    for u5 in [0..2] do
    for u6 in [0..2] do
        Append(~digs, [u1, u2, u3, u4, u5, u6]);
    end for;
    end for;
    end for;
    end for;
    end for;
    end for;
    return digs;
end function;

DigitCorrections := AllDigitCorrections();

function StateKey(vmod, scaled)
    return Sprint(vmod) cat "|" cat Sprint(scaled);
end function;

procedure AddState(~states, vmod, scaled, wt)
    key := StateKey(vmod, scaled);
    if IsDefined(states, key) then
        old := states[key];
        states[key] := <old[1], old[2], old[3] + wt>;
    else
        states[key] := <vmod, scaled, wt>;
    end if;
end procedure;

function MakeCompressedStates(polys, residues, start_k, target_k)
    rem := target_k - start_k;
    M := 3^rem;
    states := AssociativeArray();
    for vals in residues do
        scaled := [ (ExactQuotient(Z!Evaluate(polys[i], vals), 3^start_k)) mod M :
                    i in [1..#polys] ];
        vmod := [ (Z!vals[i]) mod M : i in [1..6] ];
        AddState(~states, vmod, scaled, 1);
    end for;
    return states;
end function;

function EvalJacobianMod(derivs, vmod, M)
    return [ [ (Z!Evaluate(derivs[i][j], vmod)) mod M : j in [1..6] ] :
             i in [1..#derivs] ];
end function;

function ObsKey(left_rows, scaled)
    obs := [];
    for ell in left_rows do
        Append(~obs, (&+[ ell[i]*(scaled[i] mod 3) : i in [1..#scaled] ]) mod 3);
    end for;
    return Sprint(obs);
end function;

function ObstructionVectorAtDepth(polys, left_rows, vals, depth_k)
    scaled := [ (ExactQuotient(Z!Evaluate(polys[i], vals), 3^depth_k)) mod 3 :
                i in [1..#polys] ];
    return [ (&+[ left_rows[j][i]*scaled[i] : i in [1..#scaled] ]) mod 3 :
             j in [1..#left_rows] ];
end function;

function FirstLiftOneStep(polys, vals, current_k, table)
    corrs := LiftCorrectionsFromTable(polys, vals, 3^current_k, table);
    if #corrs eq 0 then
        return false, vals;
    end if;
    d := corrs[1];
    return true, [ vals[i] + 3^current_k*d[i] : i in [1..6] ];
end function;

function FirstLiftToDepth(polys, vals, current_k, target_k, table)
    cur := vals;
    if current_k gt target_k then
        return false, vals;
    end if;
    if current_k eq target_k then
        return true, cur;
    end if;
    for kk in [current_k..target_k - 1] do
        ok, cur_next := FirstLiftOneStep(polys, cur, kk, table);
        if not ok then
            return false, vals;
        end if;
        cur := cur_next;
    end for;
    return true, cur;
end function;

function RepresentativeAncestorScan(polys, left_rows, parents, start_k, target_k,
                                    start_vals)
    table := CorrectionTable(polys, start_vals);
    good_parents := [];
    obs_hist := AssociativeArray();
    no_rep := 0;
    for vals in parents do
        ok, rep := FirstLiftToDepth(polys, vals, start_k, target_k - 1, table);
        if ok then
            obs := ObstructionVectorAtDepth(polys, left_rows, rep, target_k - 1);
            Increment(~obs_hist, Sprint(obs));
            if &and [ obs[i] eq 0 : i in [1..#obs] ] then
                Append(~good_parents, vals);
            end if;
        else
            no_rep +:= 1;
        end if;
    end for;
    inferred_lifts := #good_parents * 27^(target_k - start_k);
    row := <target_k, 3^target_k, inferred_lifts, #parents, #good_parents,
            no_rep, SortedHist(obs_hist), "representative_parent_recurrence">;
    return good_parents, row;
end function;

function RepresentativeRows(polys, left_rows, residues, start_k, target_k,
                            start_vals)
    rows := [];
    parents := residues;
    for kk in [start_k + 2..target_k] do
        parents, row := RepresentativeAncestorScan(polys, left_rows, parents,
                                                   start_k, kk, start_vals);
        Append(~rows, row);
        if #parents eq 0 then
            break;
        end if;
    end for;
    return rows;
end function;

procedure SampleNextObstruction(polys, left_rows, residues, current_k, start_vals,
                                parent_limit)
    if #residues eq 0 then
        print " sampled_next_obstruction", "no_residues";
        return;
    end if;
    n := Minimum(parent_limit, #residues);
    table := CorrectionTable(polys, start_vals);
    base_mod := 3^current_k;
    parent_hist := AssociativeArray();
    child_obs_hist := AssociativeArray();
    printed_examples := 0;
    print " sampled_next_obstruction_from_depth", current_k,
          "sampled_parent_count", n, "of", #residues;
    for s in [1..n] do
        idx := (n eq 1) select 1 else 1 + ((s - 1)*(#residues - 1)) div (n - 1);
        vals := residues[idx];
        corrs := LiftCorrectionsFromTable(polys, vals, base_mod, table);
        liftable_child_digits := [];
        child_liftable := 0;
        for d in corrs do
            child := [ vals[i] + base_mod*d[i] : i in [1..6] ];
            obs := ObstructionVectorAtDepth(polys, left_rows, child, current_k + 1);
            Increment(~child_obs_hist, Sprint(obs));
            if &and [ obs[i] eq 0 : i in [1..#obs] ] then
                child_liftable +:= 1;
                Append(~liftable_child_digits, [ F3!d[i] : i in [1..6] ]);
            end if;
        end for;
        Increment(~parent_hist, Sprint(child_liftable));
        if printed_examples lt 3 and #liftable_child_digits gt 0 then
            printed_examples +:= 1;
            basis := AffineLinearBasis(liftable_child_digits);
            print "  sampled_parent_index", idx,
                  "next_depth_liftable_child_digits", child_liftable,
                  "affine_basis_for_liftable_next_digit_subset";
            for rel in basis do
                print "   ", LinearEquationString(rel, DigitNames);
            end for;
        end if;
    end for;
    print " sampled_parent_liftable_child_count_hist", SortedHist(parent_hist);
    print " sampled_child_obstruction_hist", SortedHist(child_obs_hist);
end procedure;

function PropagateCompressedStates(derivs, left_rows, states, current_k, target_k)
    rem := target_k - current_k;
    M := 3^rem;
    Mnext := (rem eq 1) select 1 else 3^(rem - 1);
    next_states := AssociativeArray();
    corr_hist := AssociativeArray();
    obs_hist := AssociativeArray();
    total_in := 0;
    total_out := 0;
    liftable_states := 0;
    dead_states := 0;

    for key in Keys(states) do
        st := states[key];
        vmod := st[1];
        scaled := st[2];
        wt := st[3];
        total_in +:= wt;
        IncrementBy(~obs_hist, ObsKey(left_rows, scaled), wt);
        J := EvalJacobianMod(derivs, vmod, M);
        corr_count := 0;
        for d in DigitCorrections do
            tmp := [ (scaled[i] + &+[ J[i][j]*d[j] : j in [1..6] ]) mod M :
                     i in [1..#scaled] ];
            if &and [ (tmp[i] mod 3) eq 0 : i in [1..#tmp] ] then
                corr_count +:= 1;
                if rem gt 1 then
                    scaled_next := [ ExactQuotient(tmp[i], 3) mod Mnext :
                                     i in [1..#tmp] ];
                    vnext := [ vmod[i] mod Mnext : i in [1..6] ];
                    AddState(~next_states, vnext, scaled_next, wt);
                end if;
            end if;
        end for;
        if corr_count eq 0 then
            dead_states +:= wt;
        else
            liftable_states +:= wt;
        end if;
        IncrementBy(~corr_hist, Sprint(corr_count), wt);
        total_out +:= wt*corr_count;
    end for;

    return next_states,
           <current_k + 1, 3^(current_k + 1), total_out,
            #Keys(states), total_in, liftable_states, dead_states,
            SortedHist(corr_hist), SortedHist(obs_hist), #Keys(next_states)>;
end function;

function CompressedCountRows(polys, derivs, residues, start_k, target_k, left_rows)
    rows := [];
    if target_k le start_k then
        return rows;
    end if;
    if start_k lt target_k - start_k then
        return [ <0, 0, 0, 0, 0, 0, 0, [], [], 0> ];
    end if;
    states := MakeCompressedStates(polys, residues, start_k, target_k);
    current_k := start_k;
    while current_k lt target_k do
        states, row := PropagateCompressedStates(derivs, left_rows, states,
                                                 current_k, target_k);
        Append(~rows, row);
        current_k +:= 1;
    end while;
    return rows;
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

function EvalIntMod(poly, vals, m)
    return (Z!Evaluate(poly, vals)) mod m;
end function;

function Contact7Data(aa, bb)
    P<x> := PolynomialRing(Q);
    h := 1 - (Q!7/Q!2)*x + (Q!aa)*x^2 + (Q!bb)*x^3;
    f := ExactQuotient(h^2 + (x - 1)^7, x^2);
    return f, h;
end function;

function SmallRationalResidue(n, m, H)
    nn := (Z!n) mod m;
    for den in [1..H] do
        if GCD(den, 3) ne 1 then
            continue;
        end if;
        for num in [-H..H] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            if ((num - nn*den) mod m) eq 0 then
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

function HeightList(H)
    base := [80, 200, 500, 1000, 2000];
    out := [];
    for h in base do
        if h le H then
            Append(~out, h);
        end if;
    end for;
    if #out eq 0 or out[#out] ne H then
        Append(~out, H);
    end if;
    return out;
end function;

function CorrectionOrder(corrs, mode)
    if mode eq 1 then
        return corrs;
    elif mode eq 2 then
        return Reverse(corrs);
    end if;
    if #corrs le 1 then
        return corrs;
    end if;
    mid := Ceiling(#corrs/2);
    return corrs[[mid..#corrs] cat [1..mid - 1]];
end function;

function FindOneLift(polys, vals, current_k, target_k, table, mode)
    if current_k eq target_k then
        return true, vals;
    end if;
    base_mod := 3^current_k;
    corrs := CorrectionOrder(LiftCorrectionsFromTable(polys, vals, base_mod, table),
                             mode);
    for d in corrs do
        new_vals := [ vals[i] + base_mod*d[i] : i in [1..6] ];
        ok, ans := FindOneLift(polys, new_vals, current_k + 1, target_k,
                               table, mode);
        if ok then
            return true, ans;
        end if;
    end for;
    return false, vals;
end function;

procedure PrintRows(rows)
    for row in rows do
        print row;
    end for;
end procedure;

print "Z35 central b=0,r=1 deep lift";
print "max_k", max_k, "direct_depth", direct_depth, "max_store", max_store,
      "recon_height", recon_height, "use_compressed", use_compressed,
      "sample_parent_limit", sample_parent_limit;
print "local convention: H=0 mod 3^k; original G_i=0 mod 3^(k+1)";
print "central directions", CentralDirs[1], CentralDirs[2];

AllTotals := AssociativeArray();
AllBranches := AssociativeArray();
for kk in [1..max_k] do
    AllTotals[Sprint(kk)] := 0;
    AllBranches[Sprint(kk)] := 0;
end for;

for t in [1,2] do
    dir := CentralDirs[t];
    S, Hs, derivs, scales := BChartScaledEquations(t);
    print "CENTER", t, "dir", dir, "chart_scales", scales;

    J3 := JacobianMod3(Hs, dir);
    left_rows := IntegerRowsFromBasis(Basis(Nullspace(J3)));
    kernel_rows := IntegerRowsFromBasis(Basis(Nullspace(Transpose(J3))));
    print " jacobian_mod3_rank", Rank(J3);
    print " left_obstruction_rows_on_scaled_residuals", left_rows;
    print " homogeneous_correction_kernel_rows", kernel_rows;
    print " correction_coset_formula",
          "every nonempty next-digit set is d0 + span(kernel_rows) over F3";

    T, Ks := BranchEquations(Hs, dir);
    sols, KF := F3Solutions(Ks);
    rank_counts := AssociativeArray();
    for pt in sols do
        Increment(~rank_counts, Sprint(RankAt(KF, pt)));
    end for;
    basis := AffineLinearBasis(sols);
    print " initial_branch_next_digit_solution_count", #sols,
          "rank_histogram", Sort([ <StringToInteger(k), rank_counts[k]> :
                                   k in Keys(rank_counts) ]);
    print " initial_next_digit_affine_basis";
    for rel in basis do
        print "  ", LinearEquationString(rel, DigitNames);
    end for;

    ddepth := Minimum(direct_depth, max_k);
    rows, stop, residues := LiftStoreToDepth(Hs, dir, ddepth, max_store);
    print " direct_lift_rows_to_depth", ddepth;
    PrintRows(rows);
    print " direct_stop", stop, "stored_residues", #residues;

    branch_counts := AssociativeArray();
    for row in rows do
        kk := row[1];
        if kk le max_k then
            branch_counts[Sprint(kk)] := row[3];
        end if;
    end for;

    if stop eq "completed" and max_k gt ddepth then
        look_row := LookaheadOneRow(Hs, residues, ddepth, dir);
        print " exact_one_step_lookahead_row";
        print look_row;
        if look_row[1] le max_k then
            branch_counts[Sprint(look_row[1])] := look_row[3];
        end if;

        SampleNextObstruction(Hs, left_rows, residues, ddepth, dir,
                              sample_parent_limit);

        representative_rows := RepresentativeRows(Hs, left_rows, residues,
                                                  ddepth, max_k, dir);
        print " representative_recurrence_rows";
        print " columns: k modulus inferred_lifts input_parent_count good_parent_count no_representative obstruction_hist status";
        PrintRows(representative_rows);
        for row in representative_rows do
            branch_counts[Sprint(row[1])] := row[3];
        end for;

        if use_compressed ne 0 then
            if ddepth lt max_k - ddepth then
                print " compressed_transition_status",
                      "not attempted: need direct_depth >= max_k-direct_depth";
            else
                print " compressed_transition_status",
                      "attempting optional exact grouped automaton";
                print " warning",
                      "this can be slow if the state partition is too fine";
                comp_rows := CompressedCountRows(Hs, derivs, residues, ddepth,
                                                max_k, left_rows);
                print " compressed_transition_rows";
                print " columns: next_k modulus total_lifts input_state_classes input_weight liftable_weight dead_weight correction_count_hist obstruction_hist next_state_classes";
                PrintRows(comp_rows);
            end if;
        else
            print " compressed_transition_status",
                  "skipped by default; use -b use_compressed:=1 for the experimental exact grouped automaton";
        end if;
    end if;

    for kk in [1..max_k] do
        key := Sprint(kk);
        if IsDefined(branch_counts, key) then
            AllTotals[key] +:= branch_counts[key];
            if branch_counts[key] gt 0 then
                AllBranches[key] +:= 1;
            end if;
        end if;
    end for;

    table := CorrectionTable(Hs, dir);
    print " sample_lifts_and_reconstruction";
    for mode in [1..3] do
        ok, sample := FindOneLift(Hs, dir, 1, max_k, table, mode);
        print " sample_mode", mode, "found", ok;
        if ok then
            print "  sample_local_mod", 3^max_k, sample;
            orig := OriginalVals(t, sample);
            print "  sample_original_mod", 3^(max_k + 1), orig;
            print "  original_G_residues_mod_3^(max_k+1)",
                  [ EvalIntMod(Gs[i], orig, 3^(max_k + 1)) : i in [1..#Gs] ];
            f, h := Contact7Data(orig[1], orig[2]);
            disc := Discriminant(f);
            print "  sample_h1", Evaluate(h, Q!1),
                  "sample_disc_v3",
                  Valuation(Numerator(disc), 3) - Valuation(Denominator(disc), 3),
                  "sample_disc_nonzero", disc ne 0;
            for H in HeightList(recon_height) do
                recon_ok, recon := ReconstructTuple(sample, 3^max_k, H);
                print "  reconstruction_height", H, "ok", recon_ok;
                if recon_ok then
                    origQ := [ Q!1 + 3*recon[1], 3*recon[2],
                               Q!t + 3*recon[3], Q!t + 3*recon[4],
                               Q!t + 3*recon[5], Q!1 + 3*recon[6] ];
                    evalsQ := EvalOriginalQ(origQ);
                    print "   reconstructed_local", recon;
                    print "   reconstructed_original", origQ;
                    print "   reconstructed_exact_zero",
                          &and [ e eq 0 : e in evalsQ ];
                end if;
            end for;
        end if;
    end for;
end for;

print "CENTRAL_TOTAL_LIFT_TABLE";
print "columns: k modulus total_lifts_two_central_branches nonzero_central_branches";
for kk in [1..max_k] do
    print kk, 3^kk, AllTotals[Sprint(kk)], AllBranches[Sprint(kk)];
end for;

print "verdict_marker",
      "central branches survive through requested precision if final table has two nonzero branches";

quit;
