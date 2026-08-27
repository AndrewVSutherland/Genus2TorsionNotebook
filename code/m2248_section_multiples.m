//////////////////////////////////////////////////////////////////////
//  Generate rational curves on M(2,2,2,8) from multiples of the
//  section P=(1,1) on the elliptic fibration A_t and score their
//  pullbacks to the M(2,2,4,8) cover.
//
//  This complements code/m2248_pullback_bl.m, which handles the named
//  section families and the Bertin--Lecacheux branch model.
//
//  Typical runs from torsion_jac:
//
//      magma -b multiples_output_file:=data/m2248_multiple_squareclasses_n8.txt \
//          max_multiple:=8 code/m2248_section_multiples.m
//
//      magma -b search_output_file:=data/m2248_multiple_search_n8_h30.txt \
//          max_multiple:=8 height:=30 full_check:=false progress_interval:=1000 \
//          code/m2248_section_multiples.m
//////////////////////////////////////////////////////////////////////

load "code/m2248_pullback_bl.m";

function ATEllipticData()
    K<T> := RationalFunctionField(Q);
    E := EllipticCurve([
        K!0,
        -((T+1) + (T+1)/T),
        K!0,
        (T+1)^2/T,
        K!0
    ]);

    P := E![1, 1];
    torsion := [
        <"O", E!0>,
        <"T1", E![0, 0]>,
        <"T2", E![T+1, 0]>,
        <"T3", E![(T+1)/T, 0]>
    ];

    return E, P, torsion, T;
end function;

function ATPointToTuple(Pt, T)
    if Pt[3] eq 0 then
        return false, [];
    end if;

    X := Pt[1]/Pt[3];
    Y := Pt[2]/Pt[3];
    denom := T*X - (T+1)^2;
    if denom eq 0 then
        return false, [];
    end if;

    a := -T*(T+1)*(T*X^2 - (T^2+T+1)*X + 2*T*Y + (T+1)^2)/denom^2;
    b := -X*T/(T+1);

    return true, [a, b, Parent(T)!1, T];
end function;

function ATSectionTuple(n, torsion_label)
    E, P, torsion, T := ATEllipticData();
    torsion_point := E!0;
    found := false;

    for data in torsion do
        if data[1] eq torsion_label then
            torsion_point := data[2];
            found := true;
        end if;
    end for;
    if not found then
        error "Unknown torsion label.";
    end if;

    Pt := n*P + torsion_point;
    return ATPointToTuple(Pt, T);
end function;

function ConstantSquareclassIsBad(sc)
    if Degree(sc) ne 0 then
        return false;
    end if;
    ok, root := IsSquare(Coefficient(sc, 0));
    return not ok;
end function;

function SquareclassScoreRows(vals)
    rows := [];

    for perm in M2248PermutationIndices4 do
        ordered := [ vals[i] : i in perm ];
        squareclasses := A2244SquareclassesForOrderedTuple(ordered);
        degrees := [ Degree(sc) : sc in squareclasses ];
        bad_constant := &or [ ConstantSquareclassIsBad(sc) : sc in squareclasses ];
        Append(~rows, <
            &+degrees,
            Max(degrees),
            #Set([ Sprint(sc) : sc in squareclasses ]),
            bad_constant,
            perm,
            degrees,
            squareclasses
        >);
    end for;

    return rows;
end function;

function EvaluateRationalFunctionAtQ(f, r)
    num := Evaluate(Numerator(f), r);
    den := Evaluate(Denominator(f), r);
    if den eq 0 then
        return false, Q!0;
    end if;

    return true, Q!num / Q!den;
end function;

function SpecializeTupleAtQ(vals, r)
    out := [];
    for f in vals do
        ok, val := EvaluateRationalFunctionAtQ(f, r);
        if not ok then
            return false, [];
        end if;
        Append(~out, val);
    end for;

    return true, out;
end function;

procedure WriteMultipleSquareclassScores(output_file, max_multiple)
    out := Open(output_file, "w");
    E, P, torsion, T := ATEllipticData();

    print "Writing multiple-section squareclass scores to", output_file;
    print "Maximum multiple:", max_multiple;

    for n in [1..max_multiple] do
        for data in torsion do
            label := data[1];
            ok, vals := ATPointToTuple(n*P + data[2], T);
            if not ok then
                print "Skipping", n, label, "because the inverse map is undefined generically.";
                fprintf out, "SECTION n=%o | torsion=%o | skipped=true\n", n, label;
                continue;
            end if;

            a2244_count, intermediate_count, full_count, witnesses :=
                M2248SymbolicCountsForTuple(vals, false);
            fprintf out,
                "SECTION n=%o | torsion=%o | skipped=false | A2244_orderings=%o | intermediate_signs=%o | full_signs=%o\n",
                n,
                label,
                a2244_count,
                intermediate_count,
                full_count;

            rows := SquareclassScoreRows(vals);
            for row in rows do
                fprintf out,
                    "n=%o | torsion=%o | total_degree=%o | max_degree=%o | distinct_classes=%o | bad_constant=%o | perm=%o | degrees=%o | squareclasses=%o\n",
                    n,
                    label,
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

    delete out;
    print "Done writing multiple-section scores.";
end procedure;

procedure SearchMultipleSectionSpecializations(output_file, max_multiple, height, full_check, progress_interval)
    params := RationalParametersOfHeight(height);
    out := Open(output_file, "w");
    E, P, torsion, T := ATEllipticData();

    total_specs := 0;
    total_sources := 0;
    total_witnesses := 0;

    print "Searching multiple sections";
    print "Maximum multiple:", max_multiple;
    print "Parameters:", #params, "height <=", height;
    print "Full check:", full_check;
    print "Output:", output_file;

    for n in [1..max_multiple] do
        for data in torsion do
            label := data[1];
            ok, vals := ATPointToTuple(n*P + data[2], T);
            if not ok then
                print "Skipping", n, label, "because the inverse map is undefined generically.";
                continue;
            end if;

            section_sources := 0;
            section_witnesses := 0;

            for idx in [1..#params] do
                r := params[idx];
                total_specs +:= 1;

                ok_spec, vals_q := SpecializeTupleAtQ(vals, r);
                if not ok_spec or not IsUsableTuple(vals_q) then
                    continue;
                end if;

                tuple := PrimitiveIntegerTuple(vals_q);
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
                        "n=%o | torsion=%o | param=%o | raw=%o | source=%o | ordered=%o | perm=%o | eps=[%o,%o,%o] | full=%o | intermediate=%o | F=[%o,%o,%o,%o]\n",
                        n,
                        label,
                        r,
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

                if progress_interval gt 0 and total_specs mod progress_interval eq 0 then
                    print "Processed", total_specs, "specializations; sources:",
                          total_sources, "witnesses:", total_witnesses;
                end if;
            end for;

            print "Section", n, label, "sources:", section_sources,
                  "witnesses:", section_witnesses;
        end for;
    end for;

    delete out;
    print "Total sources:", total_sources;
    print "Total witnesses:", total_witnesses;
    print "Done.";
end procedure;

if assigned max_multiple then
    if Type(max_multiple) eq MonStgElt then
        max_multiple := StringToInteger(max_multiple);
    end if;
else
    max_multiple := 8;
end if;

if assigned multiples_output_file then
    WriteMultipleSquareclassScores(multiples_output_file, max_multiple);
end if;

if assigned search_output_file then
    if not assigned height then
        height := 30;
    end if;
    if not assigned full_check then
        full_check := false;
    end if;
    if not assigned progress_interval then
        progress_interval := 1000;
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

    SearchMultipleSectionSpecializations(
        search_output_file,
        max_multiple,
        height,
        full_check,
        progress_interval
    );
end if;
