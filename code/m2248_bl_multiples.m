//////////////////////////////////////////////////////////////////////
//  Generate rational curves on the Bertin--Lecacheux fibration of
//  S = M_1(8,2,2,2), score their normalized [a,b,c,d] pullbacks, and
//  search for M(2,2,4,8) witnesses.
//
//  The Magma Weierstrass convention used here is
//      y^2 + (s^2-2s+2)xy = x(x-1)(x-s^2),
//  which is the convention for which P+(0,0) matches the section written
//  in NotesAndTodo.tex.  With the maps from the notes and v=-t, this
//  reproduces the existing BLuv(s) formula in code/m2248_pullback_bl.m.
//
//  Typical runs from torsion_jac:
//
//      magma -b bl_output_file:=data/m2248_bl_multiple_squareclasses_n5.txt \
//          max_multiple:=5 code/m2248_bl_multiples.m
//
//      magma -b bl_search_output_file:=data/m2248_bl_multiple_search_n8_h50.txt \
//          max_multiple:=8 height:=50 full_check:=false progress_interval:=500 \
//          code/m2248_bl_multiples.m
//////////////////////////////////////////////////////////////////////

load "code/m2248_section_multiples.m";

function BLEllipticData()
    K<S> := RationalFunctionField(Q);
    A := S^2 - 2*S + 2;
    E := EllipticCurve([
        A,
        -(1+S^2),
        K!0,
        S^2,
        K!0
    ]);

    P := E![1, 0];
    torsion := [
        <"O", E!0>,
        <"T", E![0, 0]>
    ];

    return E, P, torsion, S;
end function;

function BLPointToUV(Pt, S)
    if Pt[3] eq 0 then
        return false, S, S;
    end if;

    xE := Pt[1]/Pt[3];
    yE := Pt[2]/Pt[3];

    den_xS := xE^2 - xE - yE;
    den_yS := yE + xE - xE^2;
    den_tS := S*(xE - 1);
    if den_xS eq 0 or den_yS eq 0 or den_tS eq 0 then
        return false, S, S;
    end if;

    xS := S + yE*S/den_xS;
    yS := yE*S/den_yS;
    tS := yE/den_tS;

    den8 := tS*xS + tS*yS + (tS-1)^2;
    if den8 eq 0 then
        return false, S, S;
    end if;

    x8 := -(tS-1)^2*tS^2/den8;
    y8 := xS*(tS-1)^2*tS^3/den8;

    xp := x8 + (tS^2 - tS);
    yp := y8 + (-tS^2 + 2*tS)*x8 + (-tS^4 + 2*tS^3 - tS^2);

    v := -tS;
    den_u := xp^2 + (v^3-v)*xp;
    if den_u eq 0 then
        return false, S, S;
    end if;

    u := ((v^2+v)*yp + ((v^3-v)*xp - v^5 - v^4 + v^3 + v^2))/den_u;

    return true, u, v;
end function;

function BLSectionUV(n, torsion_label)
    E, P, torsion, S := BLEllipticData();
    torsion_point := E!0;
    found := false;

    for data in torsion do
        if data[1] eq torsion_label then
            torsion_point := data[2];
            found := true;
        end if;
    end for;
    if not found then
        error "Unknown BL torsion label.";
    end if;

    return BLPointToUV(n*P + torsion_point, S);
end function;

function BLQuadraticOverParent(u, v)
    K := Parent(u);
    R<X> := PolynomialRing(K);
    return
        -(u^2 + u*v + v^2 + u + v + 1)*X^2
        + (v+1)*(u+1)*(u+v)*X
        - u*v*(u+v+1);
end function;

function BLNormalizedTuplesFromUVSymbolic(u, v)
    K := Parent(u);
    q := BLQuadraticOverParent(u, v);
    qroots := Roots(q);

    if #qroots ne 2 then
        return [];
    end if;
    if qroots[1][2] ne 1 or qroots[2][2] ne 1 then
        return [];
    end if;

    roots := [ K!1, u, v, -u-v-1, qroots[1][1], qroots[2][1] ];
    if not HasDistinctEntries(roots) then
        return [];
    end if;

    out := [];
    seen := {};

    for i_inf in [1..6] do
        for i_zero in [1..6] do
            if i_inf eq i_zero then
                continue;
            end if;

            remaining := [ i : i in [1..6] | i ne i_inf and i ne i_zero ];
            lambdas := [
                (roots[i] - roots[i_zero])/(roots[i] - roots[i_inf])
                : i in remaining
            ];

            if &or [ lam eq 0 : lam in lambdas ] then
                continue;
            end if;

            ref := lambdas[1];
            tuple := [ K!1 ];
            ok := true;

            for j in [2..#lambdas] do
                okj, rj := IsSquare(lambdas[j]/ref);
                if not okj then
                    ok := false;
                    break;
                end if;
                Append(~tuple, rj);
            end for;

            if not ok or not IsUsableTuple(tuple) then
                continue;
            end if;

            key := Sprint(tuple);
            if key notin seen then
                Include(~seen, key);
                Append(~out, <[i_inf, i_zero], tuple>);
            end if;
        end for;
    end for;

    return out;
end function;

function BLNormalizedTuplesFromUVQ(u, v)
    q := BLQuadratic(u, v);
    qroots := Roots(q);

    if #qroots ne 2 then
        return [];
    end if;
    if qroots[1][2] ne 1 or qroots[2][2] ne 1 then
        return [];
    end if;

    roots := [ Q!1, u, v, -u-v-1, qroots[1][1], qroots[2][1] ];
    if not HasDistinctEntries(roots) then
        return [];
    end if;

    out := [];
    seen := {};

    for i_inf in [1..6] do
        for i_zero in [1..6] do
            if i_inf eq i_zero then
                continue;
            end if;

            remaining := [ i : i in [1..6] | i ne i_inf and i ne i_zero ];
            lambdas := [
                (roots[i] - roots[i_zero])/(roots[i] - roots[i_inf])
                : i in remaining
            ];

            if &or [ lam eq 0 : lam in lambdas ] then
                continue;
            end if;

            ref := lambdas[1];
            tuple := [ Q!1 ];
            ok := true;

            for j in [2..#lambdas] do
                okj, rj := IsSquare(lambdas[j]/ref);
                if not okj then
                    ok := false;
                    break;
                end if;
                Append(~tuple, rj);
            end for;

            if not ok or not IsUsableTuple(tuple) then
                continue;
            end if;

            key := Sprint(tuple);
            if key notin seen then
                Include(~seen, key);
                Append(~out, <[i_inf, i_zero], tuple>);
            end if;
        end for;
    end for;

    return out;
end function;

function BLValidateKnownSection()
    E, P, torsion, S := BLEllipticData();
    ok, u, v := BLPointToUV(P + torsion[2][2], S);
    if not ok then
        return false;
    end if;

    u0, v0 := BLuv(S);
    return u eq u0 and v eq v0;
end function;

procedure WriteBLMultipleSquareclassScores(output_file, max_multiple)
    out := Open(output_file, "w");
    E, P, torsion, S := BLEllipticData();

    print "BL known-section validation:", BLValidateKnownSection();
    print "Writing BL multiple-section squareclass scores to", output_file;
    print "Maximum multiple:", max_multiple;

    for n in [1..max_multiple] do
        for data in torsion do
            label := data[1];
            ok, u, v := BLPointToUV(n*P + data[2], S);
            if not ok then
                print "Skipping BL", n, label, "because a map is undefined generically.";
                fprintf out, "SECTION n=%o | torsion=%o | skipped=true\n", n, label;
                continue;
            end if;

            tuples := BLNormalizedTuplesFromUVSymbolic(u, v);
            fprintf out,
                "SECTION n=%o | torsion=%o | skipped=false | normalized_tuples=%o\n",
                n,
                label,
                #tuples;

            print "BL section", n, label, "normalized symbolic tuples:", #tuples;

            for tuple_idx in [1..#tuples] do
                branch_choice := tuples[tuple_idx][1];
                vals := tuples[tuple_idx][2];
                a2244_count, intermediate_count, full_count, witnesses :=
                    M2248SymbolicCountsForTuple(vals, false);
                fprintf out,
                    "  TUPLE tuple_idx=%o | branch_choice=%o | A2244_orderings=%o | intermediate_signs=%o | full_signs=%o\n",
                    tuple_idx,
                    branch_choice,
                    a2244_count,
                    intermediate_count,
                    full_count;

                rows := SquareclassScoreRows(vals);
                for row in rows do
                    fprintf out,
                        "n=%o | torsion=%o | tuple_idx=%o | branch_choice=%o | total_degree=%o | max_degree=%o | distinct_classes=%o | bad_constant=%o | perm=%o | degrees=%o | squareclasses=%o\n",
                        n,
                        label,
                        tuple_idx,
                        branch_choice,
                        row[1],
                        row[2],
                        row[3],
                        row[4],
                        row[5],
                        row[6],
                        row[7];
                end for;
            end for;
        end for;
    end for;

    delete out;
    print "Done writing BL multiple-section scores.";
end procedure;

procedure SearchBLMultipleSpecializations(output_file, max_multiple, height, full_check, progress_interval)
    params := RationalParametersOfHeight(height);
    out := Open(output_file, "w");
    E, P, torsion, S := BLEllipticData();

    total_specs := 0;
    total_normalized := 0;
    total_sources := 0;
    total_witnesses := 0;

    print "BL known-section validation:", BLValidateKnownSection();
    print "Searching BL multiple sections";
    print "Maximum multiple:", max_multiple;
    print "Parameters:", #params, "height <=", height;
    print "Full check:", full_check;
    print "Output:", output_file;

    for n in [1..max_multiple] do
        for data in torsion do
            label := data[1];
            ok, ufun, vfun := BLPointToUV(n*P + data[2], S);
            if not ok then
                print "Skipping BL", n, label, "because a map is undefined generically.";
                continue;
            end if;

            section_normalized := 0;
            section_sources := 0;
            section_witnesses := 0;

            for idx in [1..#params] do
                r := params[idx];
                total_specs +:= 1;

                ok_u, u := EvaluateRationalFunctionAtQ(ufun, r);
                ok_v, v := EvaluateRationalFunctionAtQ(vfun, r);
                if not (ok_u and ok_v) then
                    continue;
                end if;

                tuples := [];
                try
                    tuples := BLNormalizedTuplesFromUVQ(u, v);
                catch e
                    tuples := [];
                end try;

                section_normalized +:= #tuples;
                total_normalized +:= #tuples;

                for tuple_data in tuples do
                    branch_choice := tuple_data[1];
                    tuple_rat := tuple_data[2];
                    tuple := PrimitiveIntegerTuple(tuple_rat);

                    witnesses := M2248WitnessesForTupleAllPermutations(tuple, full_check);
                    if #witnesses eq 0 then
                        continue;
                    end if;

                    section_sources +:= 1;
                    section_witnesses +:= #witnesses;
                    total_sources +:= 1;
                    total_witnesses +:= #witnesses;

                    for W in witnesses do
                        fprintf out,
                            "n=%o | torsion=%o | param=%o | branch_choice=%o | raw=%o | source=%o | ordered=%o | perm=%o | eps=[%o,%o,%o] | full=%o | intermediate=%o | F=[%o,%o,%o,%o]\n",
                            n,
                            label,
                            r,
                            branch_choice,
                            tuple,
                            W`source_tuple,
                            W`tuple,
                            W`permutation,
                            W`eps_rho,
                            W`eps_sigma,
                            W`eps_tau,
                            W`full_pass,
                            W`intermediate,
                            W`F1,
                            W`F2,
                            W`F3,
                            W`F4;
                    end for;
                end for;

                if progress_interval gt 0 and total_specs mod progress_interval eq 0 then
                    print "Processed", total_specs, "BL specializations; normalized:",
                          total_normalized, "sources:", total_sources,
                          "witnesses:", total_witnesses;
                end if;
            end for;

            print "BL section", n, label,
                  "normalized:", section_normalized,
                  "sources:", section_sources,
                  "witnesses:", section_witnesses;
        end for;
    end for;

    delete out;
    print "BL normalized tuples:", total_normalized;
    print "BL total sources:", total_sources;
    print "BL total witnesses:", total_witnesses;
    print "Done.";
end procedure;

if assigned max_multiple then
    if Type(max_multiple) eq MonStgElt then
        max_multiple := StringToInteger(max_multiple);
    end if;
else
    max_multiple := 5;
end if;

if assigned bl_output_file then
    WriteBLMultipleSquareclassScores(bl_output_file, max_multiple);
end if;

if assigned bl_search_output_file then
    if not assigned height then
        height := 50;
    end if;
    if not assigned full_check then
        full_check := false;
    end if;
    if not assigned progress_interval then
        progress_interval := 500;
    end if;

    if Type(height) eq MonStgElt then
        height := StringToInteger(height);
    end if;
    if Type(full_check) eq MonStgElt then
        full_check := full_check eq "true" or full_check eq "True" or full_check eq "1";
    end if;
    if Type(progress_interval) eq MonStgElt then
        progress_interval := StringToInteger(progress_interval);
    end if;

    SearchBLMultipleSpecializations(
        bl_search_output_file,
        max_multiple,
        height,
        full_check,
        progress_interval
    );
end if;
