//////////////////////////////////////////////////////////////////////
//  Z/35 lane: compressed obstruction automaton for the two central
//  b=0,r=1 branches.
//
//  The visible obstruction state is
//
//      <t, affine next-digit constants c in F_3^3,
//          left-obstruction residuals o in F_3^2>.
//
//  For certified counting to a fixed target K, the implementation carries
//  the exact finite tail
//
//      x mod 3^(K-k),  H(x)/3^k mod 3^(K-k),
//
//  once k >= K-k.  This makes the grouped transition exact: the next digit
//  is evaluated only on the affine coset d0 + span(kernel_rows), never by
//  expanding all raw descendants.
//
//  Typical run:
//
//      magma -b max_k:=10 \
//          code/agent_Z35_compressed_automaton.m \
//          > results/Z35_compressed_k10.log
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned max_k then
    max_k := 10;
elif Type(max_k) eq MonStgElt then
    max_k := StringToInteger(max_k);
end if;

if not assigned use_finite_tail then
    use_finite_tail := 0;
elif Type(use_finite_tail) eq MonStgElt then
    use_finite_tail := StringToInteger(use_finite_tail);
end if;

if not assigned seed_depth then
    seed_depth := (max_k ge 8 and use_finite_tail eq 0) select 6
                  else Ceiling(max_k/2);
elif Type(seed_depth) eq MonStgElt then
    seed_depth := StringToInteger(seed_depth);
end if;

if not assigned print_conflict_limit then
    print_conflict_limit := 8;
elif Type(print_conflict_limit) eq MonStgElt then
    print_conflict_limit := StringToInteger(print_conflict_limit);
end if;

if not assigned print_state_limit then
    print_state_limit := 18;
elif Type(print_state_limit) eq MonStgElt then
    print_state_limit := StringToInteger(print_state_limit);
end if;

Z := Integers();
F3 := GF(3);

R0<a,b,c0,c1,c2,r> := PolynomialRing(Z, 6);

G4 := c2^2 - b^2 - 5*r + 7;
G3 := 2*c1*c2 - 2*a*b - 21 + 10*r^2;
G2 := c1^2 + 2*c0*c2 - a^2 + 7*b + 35 - 10*r^3;
G1 := 2*c0*c1 + 7*a - 2*b - 35 + 5*r^4;
G0 := 4*c0^2 - 8*a + 35 - 4*r^5;
Gs := [G4, G3, G2, G1, G0];

DigitNames := ["dA", "dB", "dC0", "dC1", "dC2", "dR"];

CentralDirs := AssociativeArray();
CentralDirs[1] := [1, 1, 1, 0, 1, 0];
CentralDirs[2] := [1, 1, 1, 2, 1, 0];

ExpectedT1 := AssociativeArray();
ExpectedT1["1"] := 1;
ExpectedT1["2"] := 27;
ExpectedT1["3"] := 729;
ExpectedT1["4"] := 6561;
ExpectedT1["5"] := 177147;
ExpectedT1["6"] := 1594323;
ExpectedT1["7"] := 43046721;

ExpectedCentralDirect := AssociativeArray();
ExpectedCentralDirect["1"] := 1;
ExpectedCentralDirect["2"] := 27;
ExpectedCentralDirect["3"] := 729;
ExpectedCentralDirect["4"] := 6561;

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

function JoinStrings(items, sep)
    if #items eq 0 then
        return "";
    end if;
    out := items[1];
    for i in [2..#items] do
        out cat:= sep cat items[i];
    end for;
    return out;
end function;

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

function LinearEquationString(coeffs, c, names)
    terms := [];
    for i in [1..#coeffs] do
        ci := Z!coeffs[i];
        if ci eq 1 then
            Append(~terms, names[i]);
        elif ci eq 2 then
            Append(~terms, "2*" cat names[i]);
        end if;
    end for;
    cc := Z!c;
    if cc eq 1 then
        Append(~terms, "1");
    elif cc eq 2 then
        Append(~terms, "2");
    end if;
    if #terms eq 0 then
        return "0 = 0";
    end if;
    return JoinStrings(terms, " + ") cat " = 0";
end function;

function DesiredEquationRows(t)
    if t eq 1 then
        return [
            [0, 0, 0, 0, 1, 2],
            [0, 1, 0, 2, 0, 1],
            [1, 0, 2, 0, 0, 1]
        ];
    end if;
    return [
        [0, 0, 0, 0, 1, 1],
        [0, 1, 0, 1, 0, 1],
        [1, 0, 1, 0, 0, 1]
    ];
end function;

function RowPreimage(J, row)
    // Find p in F_3^5 with p*J = row.  The search space has size 3^5.
    target := [ F3!(row[i] mod 3) : i in [1..#row] ];
    for u1 in F3 do
    for u2 in F3 do
    for u3 in F3 do
    for u4 in F3 do
    for u5 in F3 do
        p := [u1, u2, u3, u4, u5];
        image := [ &+[ p[i]*J[i][j] : i in [1..5] ] : j in [1..6] ];
        if image eq target then
            return true, [ Z!p[i] : i in [1..5] ];
        end if;
    end for;
    end for;
    end for;
    end for;
    end for;
    return false, [];
end function;

function StateTransformRows(t, J, left_rows)
    desired := DesiredEquationRows(t);
    rows := [];
    for erow in desired do
        ok, prow := RowPreimage(J, erow);
        if not ok then
            error "could not express desired equation row in Jacobian rowspace";
        end if;
        Append(~rows, prow);
    end for;
    for lrow in left_rows do
        Append(~rows, [ Z!(lrow[i] mod 3) : i in [1..#lrow] ]);
    end for;
    U := Matrix(F3, rows);
    if Rank(U) ne 5 then
        error "state transform rows are not invertible";
    end if;
    return rows;
end function;

function ApplyStateTransform(Urows, scaled_mod3)
    return [
        (&+[ (Urows[i][j] mod 3)*(scaled_mod3[j] mod 3) : j in [1..5] ]) mod 3 :
        i in [1..5]
    ];
end function;

function FormalStateFromScaled(t, Urows, scaled)
    smod := [ (Z!scaled[i]) mod 3 : i in [1..5] ];
    u := ApplyStateTransform(Urows, smod);
    return <t, u[[1..3]], u[[4..5]]>;
end function;

function FormalKey(st)
    ccode := st[2][1] + 3*st[2][2] + 9*st[2][3];
    ocode := st[3][1] + 3*st[3][2];
    return (st[1] - 1)*243 + 9*ccode + ocode;
end function;

function FormalPrintable(st)
    return "<t=" cat Sprint(st[1]) cat ", c=" cat Sprint(st[2]) cat
           ", obs=" cat Sprint(st[3]) cat ">";
end function;

function IsZeroObs(st)
    return &and [ st[3][i] eq 0 : i in [1..#st[3]] ];
end function;

function ObsCode(st)
    return st[3][1] + 3*st[3][2];
end function;

function EncodeVec(vals, M)
    code := 0;
    place := 1;
    for v in vals do
        code +:= (Z!v mod M)*place;
        place *:= M;
    end for;
    return code;
end function;

function ExactStateKey(t, vmod, scaled, M)
    return (t - 1)*M^11 + EncodeVec(vmod, M)*M^5 + EncodeVec(scaled, M);
end function;

procedure AddExactState(~states, t, vmod, scaled, wt, M)
    key := ExactStateKey(t, vmod, scaled, M);
    if IsDefined(states, key) then
        old := states[key];
        states[key] := <old[1], old[2], old[3], old[4] + wt>;
    else
        states[key] := <t, vmod, scaled, wt>;
    end if;
end procedure;

function EvalJacobianMod(derivs, vmod, M)
    return [ [ (Z!Evaluate(derivs[i][j], vmod)) mod M : j in [1..6] ] :
             i in [1..#derivs] ];
end function;

function DotMod(row, vec, m)
    return (&+[ (row[i] mod m)*(vec[i] mod m) : i in [1..#row] ]) mod m;
end function;

function CosetKey(t, constants)
    return Sprint(t) cat "|" cat Sprint(constants);
end function;

function CosetDigitsRaw(t, constants, kernel_rows)
    erows := DesiredEquationRows(t);
    d0 := [];
    found := false;
    for u1 in [0..2] do
    for u2 in [0..2] do
    for u3 in [0..2] do
    for u4 in [0..2] do
    for u5 in [0..2] do
    for u6 in [0..2] do
        d := [u1, u2, u3, u4, u5, u6];
        if &and [ (DotMod(erows[i], d, 3) + constants[i]) mod 3 eq 0 :
                  i in [1..3] ] then
            d0 := d;
            found := true;
            break u6;
        end if;
    end for;
    end for;
    end for;
    end for;
    end for;
    end for;
    if not found then
        return [];
    end if;

    digits := [];
    seen := AssociativeArray();
    for a1 in [0..2] do
    for a2 in [0..2] do
    for a3 in [0..2] do
        d := [
            (d0[j] + a1*kernel_rows[1][j] + a2*kernel_rows[2][j] +
             a3*kernel_rows[3][j]) mod 3 : j in [1..6]
        ];
        dkey := Sprint(d);
        if not IsDefined(seen, dkey) then
            seen[dkey] := true;
            Append(~digits, d);
        end if;
    end for;
    end for;
    end for;
    Sort(~digits);
    return digits;
end function;

function PrecomputeCosetCache(t, kernel_rows)
    cache := AssociativeArray();
    for c1 in [0..2] do
    for c2 in [0..2] do
    for c3 in [0..2] do
        constants := [c1, c2, c3];
        cache[CosetKey(t, constants)] := CosetDigitsRaw(t, constants, kernel_rows);
    end for;
    end for;
    end for;
    return cache;
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

function LiftStoreToDepth(polys, start_vals, target_k)
    residues := [ [ Z!x : x in start_vals ] ];
    base_mod := 3;
    table := CorrectionTable(polys, start_vals);
    rows := [ <1, 3, 1, [ <27, 0> ], true> ];

    for kk in [2..target_k] do
        next_residues := [];
        hist := AssociativeArray();
        total := 0;
        for vals in residues do
            corrs := LiftCorrectionsFromTable(polys, vals, base_mod, table);
            Increment(~hist, Sprint(#corrs));
            total +:= #corrs;
            for d in corrs do
                Append(~next_residues, [ vals[i] + base_mod*d[i] : i in [1..6] ]);
            end for;
        end for;
        Append(~rows, <kk, 3^kk, total,
                       Sort([ <StringToInteger(k), hist[k]> : k in Keys(hist) ]),
                       true>);
        residues := next_residues;
        base_mod *:= 3;
    end for;
    return rows, residues;
end function;

function SeedExactStates(t, Hs, residues, k, target_k)
    rem := target_k - k;
    M := 3^rem;
    states := AssociativeArray();
    for vals in residues do
        scaled := [
            (ExactQuotient(Z!Evaluate(Hs[i], vals), 3^k)) mod M :
            i in [1..#Hs]
        ];
        vmod := [ (Z!vals[i]) mod M : i in [1..6] ];
        AddExactState(~states, t, vmod, scaled, 1, M);
    end for;
    return states;
end function;

function SignatureKey(hist)
    return Sprint(SortedHist(hist));
end function;

function ExactTransitionRows(t, Hs, derivs, Urows, kernel_rows, seed_states,
                             seed_k, target_k, conflicts, all_none_bad)
    states := seed_states;
    rows := [];
    counts := AssociativeArray();
    counts[Sprint(seed_k)] := &+[ states[key][4] : key in Keys(states) ];
    coset_cache := PrecomputeCosetCache(t, kernel_rows);
    projected_signature_seen := AssociativeArray();

    current_k := seed_k;
    while current_k lt target_k do
        rem := target_k - current_k;
        M := 3^rem;
        Mnext := 3^(rem - 1);
        next_states := AssociativeArray();
        formal_hist := AssociativeArray();
        corr_hist := AssociativeArray();
        child_good_hist := AssociativeArray();
        child_obs_hist := AssociativeArray();
        total_in := 0;
        total_out := 0;
        liftable_weight := 0;
        dead_weight := 0;
        projected_conflicts_this_row := 0;
        all_none_counter_this_row := 0;

        for key in Keys(states) do
            st := states[key];
            vmod := st[2];
            scaled := st[3];
            wt := st[4];
            total_in +:= wt;

            fstate := FormalStateFromScaled(t, Urows, scaled);
            fkey := FormalKey(fstate);
            IncrementBy(~formal_hist, fkey, wt);

            need_signature := #conflicts lt print_conflict_limit;
            child_formal_hist_unweighted := AssociativeArray();
            child_good_count := 0;
            corr_count := 0;

            if IsZeroObs(fstate) then
                digits := coset_cache[CosetKey(t, fstate[2])];
                J := EvalJacobianMod(derivs, vmod, M);
                for d in digits do
                    tmp := [
                        (scaled[i] + &+[ J[i][j]*d[j] : j in [1..6] ]) mod M :
                        i in [1..#scaled]
                    ];
                    if &and [ (tmp[i] mod 3) eq 0 : i in [1..#tmp] ] then
                        corr_count +:= 1;
                        scaled_next := [ ExactQuotient(tmp[i], 3) mod Mnext :
                                         i in [1..#tmp] ];
                        vnext := [ vmod[i] mod Mnext : i in [1..6] ];
                        AddExactState(~next_states, t, vnext, scaled_next, wt,
                                      Mnext);

                        child_state := FormalStateFromScaled(t, Urows, scaled_next);
                        if need_signature then
                            child_key := FormalKey(child_state);
                            Increment(~child_formal_hist_unweighted, child_key);
                        end if;
                        IncrementBy(~child_obs_hist, ObsCode(child_state), wt);
                        if IsZeroObs(child_state) then
                            child_good_count +:= 1;
                        end if;
                    end if;
                end for;
            end if;

            if corr_count eq 0 then
                dead_weight +:= wt;
            else
                liftable_weight +:= wt;
            end if;
            total_out +:= wt*corr_count;
            IncrementBy(~corr_hist, corr_count, wt);
            IncrementBy(~child_good_hist, child_good_count, wt);

            if corr_count gt 0 and child_good_count ne 0 and child_good_count ne corr_count then
                all_none_counter_this_row +:= 1;
                if #all_none_bad lt print_conflict_limit then
                    Append(~all_none_bad,
                           <current_k, fstate, vmod, scaled, child_good_count,
                            corr_count, SortedHist(child_formal_hist_unweighted)>);
                end if;
            end if;

            if need_signature then
                sig := SignatureKey(child_formal_hist_unweighted);
                proj_key := current_k*1000 + fkey;
                if IsDefined(projected_signature_seen, proj_key) then
                    old := projected_signature_seen[proj_key];
                    if old[1] ne sig then
                        projected_conflicts_this_row +:= 1;
                        Append(~conflicts,
                               <current_k, fstate, old[2], old[3], vmod, scaled,
                                old[4], child_good_count>);
                    end if;
                else
                    projected_signature_seen[proj_key] :=
                        <sig, vmod, scaled, child_good_count>;
                end if;
            end if;
        end for;

        current_k +:= 1;
        counts[Sprint(current_k)] := total_out;
        Append(~rows,
               <current_k, 3^current_k, total_out, #Keys(states),
                total_in, liftable_weight, dead_weight,
                SortedHist(corr_hist), SortedHist(child_good_hist),
                SortedHist(child_obs_hist), #Keys(formal_hist),
                #Keys(next_states), projected_conflicts_this_row,
                all_none_counter_this_row>);
        states := next_states;
    end while;

    return rows, counts, states, conflicts, all_none_bad;
end function;

function FastCosetLookahead(t, Hs, derivs, Urows, kernel_rows, residues,
                            current_k, conflicts)
    // Certifies counts at current_k+1 and current_k+2 by scanning every
    // stored parent, grouped by the exact mod-9 tail needed for this
    // two-step calculation, and only its 27 affine next digits.
    coset_cache := PrecomputeCosetCache(t, kernel_rows);
    M := 9;
    fast_states := AssociativeArray();
    for vals in residues do
        scaled := [
            (ExactQuotient(Z!Evaluate(Hs[i], vals), 3^current_k)) mod M :
            i in [1..#Hs]
        ];
        vmod := [ vals[i] mod M : i in [1..6] ];
        AddExactState(~fast_states, t, vmod, scaled, 1, M);
    end for;

    total_next := 0;
    total_next2 := 0;
    parent_count := #residues;
    liftable_parent_count := 0;
    dead_parent_count := 0;
    good_parent_count := 0;
    corr_hist := AssociativeArray();
    child_good_hist := AssociativeArray();
    child_obs_hist := AssociativeArray();
    all_none_counter := 0;
    good_parents := [];

    for key in Keys(fast_states) do
        st := fast_states[key];
        vals := st[2];
        scaled := st[3];
        wt := st[4];
        fstate := FormalStateFromScaled(t, Urows, scaled);
        corr_count := 0;
        child_good_count := 0;
        if IsZeroObs(fstate) then
            digits := coset_cache[CosetKey(t, fstate[2])];
            J := EvalJacobianMod(derivs, vals, M);
            for d in digits do
                tmp := [
                    (scaled[i] + &+[ J[i][j]*d[j] : j in [1..6] ]) mod M :
                    i in [1..#scaled]
                ];
                if &and [ (tmp[i] mod 3) eq 0 : i in [1..#tmp] ] then
                    corr_count +:= 1;
                    scaled_child := [ ExactQuotient(tmp[i], 3) mod 3 :
                                      i in [1..#tmp] ];
                    child_state := FormalStateFromScaled(t, Urows, scaled_child);
                    IncrementBy(~child_obs_hist, ObsCode(child_state), wt);
                    if IsZeroObs(child_state) then
                        child_good_count +:= 1;
                    end if;
                end if;
            end for;
        end if;

        total_next +:= wt*corr_count;
        total_next2 +:= wt*27*child_good_count;
        IncrementBy(~corr_hist, corr_count, wt);
        IncrementBy(~child_good_hist, child_good_count, wt);
        if corr_count eq 0 then
            dead_parent_count +:= wt;
        else
            liftable_parent_count +:= wt;
        end if;
        if child_good_count gt 0 then
            good_parent_count +:= wt;
        end if;
        if corr_count gt 0 and child_good_count ne 0 and child_good_count ne corr_count then
            all_none_counter +:= wt;
            if #conflicts lt print_conflict_limit then
                Append(~conflicts,
                       <current_k, fstate, vals, child_good_count, corr_count>);
            end if;
        end if;
    end for;

    row1 := <current_k + 1, 3^(current_k + 1), total_next,
             parent_count, liftable_parent_count, dead_parent_count,
             #Keys(fast_states), SortedHist(corr_hist),
             "certified_all_parent_cosets_grouped_mod9">;
    row2 := <current_k + 2, 3^(current_k + 2), total_next2,
             parent_count, good_parent_count, #Keys(fast_states),
             SortedHist(child_good_hist),
             SortedHist(child_obs_hist), all_none_counter,
             "certified_child_obstruction_scan_grouped_mod9">;
    return row1, row2, good_parents, conflicts;
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

function RepresentativeRows(polys, left_rows, parents, start_k, target_k,
                            start_vals)
    table := CorrectionTable(polys, start_vals);
    rows := [];
    current_parents := parents;
    for kk in [start_k + 2..target_k] do
        good_parents := [];
        obs_hist := AssociativeArray();
        no_rep := 0;
        for vals in current_parents do
            ok, rep := FirstLiftToDepth(polys, vals, start_k, kk - 1, table);
            if ok then
                obs := ObstructionVectorAtDepth(polys, left_rows, rep, kk - 1);
                Increment(~obs_hist, ObsCode(<0, [], obs>));
                if &and [ obs[i] eq 0 : i in [1..#obs] ] then
                    Append(~good_parents, vals);
                end if;
            else
                no_rep +:= 1;
            end if;
        end for;
        inferred_lifts := #good_parents * 27^(kk - start_k);
        Append(~rows, <kk, 3^kk, inferred_lifts, #current_parents,
                       #good_parents, no_rep, SortedHist(obs_hist),
                       "representative_not_certified">);
        current_parents := good_parents;
        if #current_parents eq 0 then
            break;
        end if;
    end for;
    return rows;
end function;

function DirectCountMap(rows)
    A := AssociativeArray();
    for row in rows do
        A[Sprint(row[1])] := row[3];
    end for;
    return A;
end function;

procedure PrintRows(rows)
    for row in rows do
        print row;
    end for;
end procedure;

print "Z35 compressed obstruction automaton";
print "max_k", max_k, "seed_depth", seed_depth,
      "use_finite_tail", use_finite_tail,
      "print_conflict_limit", print_conflict_limit;
print "local convention: H=0 mod 3^k; original G_i=0 mod 3^(k+1)";
if use_finite_tail ne 0 and seed_depth lt max_k - seed_depth then
    error "need seed_depth >= max_k - seed_depth for certified finite-tail propagation";
end if;

AllCounts := AssociativeArray();
for kk in [1..max_k] do
    AllCounts[Sprint(kk)] := 0;
end for;

AllBranchCounts := AssociativeArray();
AllBranchStatus := AssociativeArray();
AllNoneBad := [];
ProjectionConflicts := [];

for t in [1,2] do
    dir := CentralDirs[t];
    S, Hs, derivs, scales := BChartScaledEquations(t);
    J3 := JacobianMod3(Hs, dir);
    left_rows := IntegerRowsFromBasis(Basis(Nullspace(J3)));
    kernel_rows := IntegerRowsFromBasis(Basis(Nullspace(Transpose(J3))));
    Urows := StateTransformRows(t, J3, left_rows);
    erows := DesiredEquationRows(t);

    print "";
    print "CENTER", t, "dir", dir, "chart_scales", scales;
    print " jacobian_mod3_rank", Rank(J3);
    print " state_definition";
    print "  affine next-digit equations E_i(d)+c_i=0:";
    for i in [1..3] do
        print "   E", i, LinearEquationString(erows[i], 0, DigitNames);
    end for;
    print "  left_obstruction_rows_on_scaled_residuals", left_rows;
    print "  homogeneous_correction_kernel_rows", kernel_rows;
    print "  state_transform_rows_on_scaled_residual", Urows;
    initial_scaled := [
        (ExactQuotient(Z!Evaluate(Hs[i], dir), 3)) mod 3 : i in [1..#Hs]
    ];
    initial_state := FormalStateFromScaled(t, Urows, initial_scaled);
    print " initial_obstruction_state_at_k1", FormalPrintable(initial_state);
    print " initial_affine_next_digit_equations";
    for i in [1..3] do
        print "  ", LinearEquationString(erows[i], initial_state[2][i], DigitNames);
    end for;

    rows_direct, residues := LiftStoreToDepth(Hs, dir, seed_depth);
    print " direct_seed_rows";
    PrintRows(rows_direct);
    direct_counts := DirectCountMap(rows_direct);

    branch_counts := AssociativeArray();
    branch_status := AssociativeArray();
    for key in Keys(direct_counts) do
        branch_counts[key] := direct_counts[key];
        branch_status[key] := "direct_seed_certified";
    end for;

    if use_finite_tail ne 0 then
        seed_states := SeedExactStates(t, Hs, residues, seed_depth, max_k);
        seed_formal_hist := AssociativeArray();
        for key in Keys(seed_states) do
            st := seed_states[key];
            fstate := FormalStateFromScaled(t, Urows, st[3]);
            IncrementBy(~seed_formal_hist, FormalKey(fstate), st[4]);
        end for;
        print " seed_exact_state_classes", #Keys(seed_states),
              "seed_formal_state_classes", #Keys(seed_formal_hist);
        printed := 0;
        print " seed_formal_state_hist_sample";
        for item in SortedHist(seed_formal_hist) do
            if printed ge print_state_limit then
                break;
            end if;
            printed +:= 1;
            print "  ", item;
        end for;

        trans_rows, trans_counts, final_states, ProjectionConflicts, AllNoneBad :=
            ExactTransitionRows(t, Hs, derivs, Urows, kernel_rows, seed_states,
                                seed_depth, max_k, ProjectionConflicts,
                                AllNoneBad);
        print " certified_grouped_transition_rows";
        print " columns: next_k modulus total_lifts exact_state_classes input_weight liftable_weight dead_weight correction_count_hist child_good_count_hist child_obstruction_hist formal_state_classes next_exact_state_classes projected_state_conflicts all_or_none_counter_state_classes";
        PrintRows(trans_rows);

        for key in Keys(trans_counts) do
            branch_counts[key] := trans_counts[key];
            branch_status[key] := "finite_tail_certified";
        end for;
    else
        if max_k gt seed_depth then
            row1, row2, good_parents, AllNoneBad :=
                FastCosetLookahead(t, Hs, derivs, Urows, kernel_rows,
                                   residues, seed_depth, AllNoneBad);
            print " certified_fast_coset_rows";
            print " columns row1: k modulus lifts parent_count liftable_parent_count dead_parent_count mod9_state_classes correction_count_hist status";
            print " columns row2: k modulus lifts parent_count good_parent_count mod9_state_classes child_good_count_hist child_obstruction_code_hist all_or_none_counter status";
            print row1;
            print row2;
            branch_counts[Sprint(row1[1])] := row1[3];
            branch_status[Sprint(row1[1])] := "fast_coset_certified";
            if row2[1] le max_k then
                branch_counts[Sprint(row2[1])] := row2[3];
                branch_status[Sprint(row2[1])] := "fast_child_obstruction_certified";
            end if;

            if max_k gt seed_depth + 2 then
                rep_rows := RepresentativeRows(Hs, left_rows, residues,
                                               seed_depth, max_k, dir);
                print " representative_recurrence_rows";
                print " columns: k modulus inferred_lifts input_parent_count good_parent_count no_representative obstruction_code_hist status";
                PrintRows(rep_rows);
                for row in rep_rows do
                    if row[1] gt seed_depth + 2 then
                        branch_counts[Sprint(row[1])] := row[3];
                        branch_status[Sprint(row[1])] := row[8];
                    end if;
                end for;
            end if;
        end if;
    end if;
    AllBranchCounts[Sprint(t)] := branch_counts;
    AllBranchStatus[Sprint(t)] := branch_status;

    print " branch_certified_counts";
    print " columns: k modulus lifts status validation";
    for kk in [1..max_k] do
        key := Sprint(kk);
        val := IsDefined(branch_counts, key) select branch_counts[key] else 0;
        status := IsDefined(branch_status, key) select branch_status[key]
                                                else "not_computed";
        if t eq 1 and IsDefined(ExpectedT1, key) then
            status cat:= (val eq ExpectedT1[key]) select "_matches_recorded_t1_exact"
                                                   else "_MISMATCH_recorded_t1_exact";
        end if;
        if IsDefined(ExpectedCentralDirect, key) then
            status cat:= (val eq ExpectedCentralDirect[key]) select "_matches_direct_3^4"
                                                           else "_MISMATCH_direct_3^4";
        end if;
        print kk, 3^kk, val, status;
        AllCounts[key] +:= val;
    end for;
end for;

print "";
print "TWO_BRANCH_TOTALS";
print "columns: k modulus total_lifts nonzero_branches aggregate_status validation";
for kk in [1..max_k] do
    key := Sprint(kk);
    nonzero := 0;
    aggregate_status := "certified";
    for t in [1,2] do
        bc := AllBranchCounts[Sprint(t)];
        if IsDefined(bc, key) and bc[key] gt 0 then
            nonzero +:= 1;
        end if;
        bs := AllBranchStatus[Sprint(t)];
        if (not IsDefined(bs, key)) or bs[key] eq "representative_not_certified" then
            aggregate_status := "has_representative_not_certified";
        end if;
    end for;
    status := aggregate_status;
    if IsDefined(ExpectedCentralDirect, key) then
        status cat:= (AllCounts[key] eq 2*ExpectedCentralDirect[key]) select "_matches_direct_3^4"
                                                                  else "_MISMATCH_direct_3^4";
    end if;
    print kk, 3^kk, AllCounts[key], nonzero, status;
end for;

print "";
print "PROJECTED_STATE_MARKOV_TEST";
if use_finite_tail eq 0 then
    print " verdict",
          "not run in fast mode; use use_finite_tail:=1 for the exact projected-state Markov test";
elif #ProjectionConflicts eq 0 then
    print " verdict", "no projected-state transition conflicts encountered through target";
else
    print " verdict", "projected obstruction state is not Markov at the tested precision";
    print " conflict_count_recorded", #ProjectionConflicts;
    for c in ProjectionConflicts do
        print "  conflict", c;
    end for;
end if;

print "";
print "ALL_OR_NONE_CHILD_COSET_TEST";
if #AllNoneBad eq 0 then
    print " verdict",
          "no counterexample through target: every encountered liftable coset has 0 or 27 liftable child states";
else
    print " verdict", "counterexamples found to all-or-none child-coset rule";
    print " counterexample_count_recorded", #AllNoneBad;
    for c in AllNoneBad do
        print "  counterexample", c;
    end for;
end if;

print "";
print "verdict_marker",
      "certified rows are explicitly labeled; finite-tail projected-state test is optional";

quit;
