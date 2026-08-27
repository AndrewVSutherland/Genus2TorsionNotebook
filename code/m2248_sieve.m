//////////////////////////////////////////////////////////////////////
//  Sieve for the extra halving cover M(2,2,4,8) -> A(2,2,4,4)
//
//  A tuple [a,b,c,d] represents the curve
//      y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2).
//
//  For tuples already satisfying the (2,2,4,4) square conditions, this
//  script tests the first intermediate double cover and, optionally, the
//  full four-square halving condition from m2248_equations_gpt55.tex.
//
//  Typical interactive use from the torsion_jac directory:
//
//      load "code/m2248_sieve.m";
//      SieveM2248TupleFile(
//          "data/tor2244_bank.txt",
//          "data/m2248_candidates.txt",
//          true, 0, 1000
//      );
//
//  Arguments to SieveM2248TupleFile are:
//      input file, output file, full_check, max_rows, progress_interval.
//
//  If full_check is false, only the intermediate condition
//      rho*sigma*tau in Q^{*2}
//  is tested.  If full_check is true, all four q_i square conditions are
//  tested.
//////////////////////////////////////////////////////////////////////

Q := Rationals();

M2248WitnessFormat := recformat<
    source_tuple,
    tuple,
    permutation,
    eps_rho,
    eps_sigma,
    eps_tau,
    rho,
    sigma,
    tau,
    intermediate,
    intermediate_root,
    F1,
    F2,
    F3,
    F4,
    q1,
    q2,
    q3,
    q4,
    full_pass
>;

M2248PermutationIndices4 := [
    [1,2,3,4], [1,2,4,3], [1,3,2,4], [1,3,4,2], [1,4,2,3], [1,4,3,2],
    [2,1,3,4], [2,1,4,3], [2,3,1,4], [2,3,4,1], [2,4,1,3], [2,4,3,1],
    [3,1,2,4], [3,1,4,2], [3,2,1,4], [3,2,4,1], [3,4,1,2], [3,4,2,1],
    [4,1,2,3], [4,1,3,2], [4,2,1,3], [4,2,3,1], [4,3,1,2], [4,3,2,1]
];

function ReadM2248TupleFile(filename)
    S := Read(filename);
    rows := Split(S, "\n");
    tuples := [];

    for rawrow in rows do
        n := #rawrow;
        if n gt 0 and rawrow[n] eq "\r" then
            n -:= 1;
        end if;

        if n ge 2 and rawrow[1] eq "[" and rawrow[n] eq "]" then
            body := rawrow[2..n-1];
            parts := Split(body, ",");
            tup := [ StringToInteger(part) : part in parts ];
            if #tup eq 4 then
                Append(~tuples, tup);
            end if;
        end if;
    end for;

    return tuples;
end function;

function SquareRootOrFail(n, label, tup)
    ok, r := IsSquare(n);
    if not ok then
        error Sprintf(
            "Tuple %o does not satisfy expected (2,2,4,4) square condition %o = %o",
            tup,
            label,
            n
        );
    end if;
    return r;
end function;

function M2248WitnessesForTuple(tup, full_check)
    if #tup ne 4 then
        error "Expected a tuple [a,b,c,d].";
    end if;

    a := Q!tup[1];
    b := Q!tup[2];
    c := Q!tup[3];
    d := Q!tup[4];

    if a eq 0 or b eq 0 or c eq 0 or d eq 0 then
        error "The tuple must have nonzero entries.";
    end if;

    A := a^2;
    B := b^2;
    C := c^2;
    D := d^2;

    if #(Set([A,B,C,D])) ne 4 then
        error "The tuple has repeated branch points.";
    end if;

    // These are the square roots appearing on A(2,2,4,4).
    // y10 and y20 are not used below, but checking them keeps this
    // function honest on arbitrary ordered tuples, not just tor2244.txt.
    y10 := SquareRootOrFail((A-C)*(A-D), "y1^2", tup);
    y20 := SquareRootOrFail((B-C)*(B-D), "y2^2", tup);
    w0 := SquareRootOrFail((A-C)*(B-C), "w^2", tup);
    t0 := SquareRootOrFail((A-D)*(B-D), "t^2", tup);

    witnesses := [];

    for eps_rho in [-1, 1] do
        for eps_sigma in [-1, 1] do
            for eps_tau in [-1, 1] do
                rho := eps_rho * a / b;
                sigma := eps_sigma * (A-C) / w0;
                tau := eps_tau * (A-D) / t0;

                intermediate := rho * sigma * tau;
                ok0, r0 := IsSquare(intermediate);
                if not ok0 then
                    continue;
                end if;

                F1 := (1 + rho) * (1 + sigma) * (1 + tau);
                F2 := rho * (1 + rho) * (rho + sigma) * (rho + tau);
                F3 := sigma * (1 + sigma) * (rho + sigma) * (sigma + tau);
                F4 := tau * (1 + tau) * (rho + tau) * (sigma + tau);

                ok1, q1 := IsSquare(F1);
                ok2, q2 := IsSquare(F2);
                ok3, q3 := IsSquare(F3);
                ok4, q4 := IsSquare(F4);
                full_pass := ok1 and ok2 and ok3 and ok4;

                if full_check and not full_pass then
                    continue;
                end if;

                Append(~witnesses, rec<M2248WitnessFormat |
                    source_tuple := tup,
                    tuple := tup,
                    permutation := [1,2,3,4],
                    eps_rho := eps_rho,
                    eps_sigma := eps_sigma,
                    eps_tau := eps_tau,
                    rho := rho,
                    sigma := sigma,
                    tau := tau,
                    intermediate := intermediate,
                    intermediate_root := r0,
                    F1 := F1,
                    F2 := F2,
                    F3 := F3,
                    F4 := F4,
                    q1 := ok1 select q1 else Q!0,
                    q2 := ok2 select q2 else Q!0,
                    q3 := ok3 select q3 else Q!0,
                    q4 := ok4 select q4 else Q!0,
                    full_pass := full_pass
                >);
            end for;
        end for;
    end for;

    return witnesses;
end function;

function M2248IntermediateWitnessesForTuple(tup)
    return M2248WitnessesForTuple(tup, false);
end function;

function M2248FullWitnessesForTuple(tup)
    return M2248WitnessesForTuple(tup, true);
end function;

function M2248OrderedTupleHasA2244Squares(tup)
    if #tup ne 4 then
        return false;
    end if;

    a := Q!tup[1];
    b := Q!tup[2];
    c := Q!tup[3];
    d := Q!tup[4];

    if a eq 0 or b eq 0 or c eq 0 or d eq 0 then
        return false;
    end if;

    A := a^2; B := b^2; C := c^2; D := d^2;
    if #(Set([A,B,C,D])) ne 4 then
        return false;
    end if;

    ok_y1, r_y1 := IsSquare((A-C)*(A-D));
    ok_y2, r_y2 := IsSquare((B-C)*(B-D));
    ok_w, r_w := IsSquare((A-C)*(B-C));
    ok_t, r_t := IsSquare((A-D)*(B-D));

    return ok_y1 and ok_y2 and ok_w and ok_t;
end function;

function M2248WitnessesForTupleAllPermutations(tup, full_check)
    if #tup ne 4 then
        error "Expected a tuple [a,b,c,d].";
    end if;

    witnesses := [];
    for perm in M2248PermutationIndices4 do
        ptup := [ tup[i] : i in perm ];
        if not M2248OrderedTupleHasA2244Squares(ptup) then
            continue;
        end if;

        perm_witnesses := M2248WitnessesForTuple(ptup, full_check);
        for W in perm_witnesses do
            W2 := W;
            W2`source_tuple := tup;
            W2`permutation := perm;
            Append(~witnesses, W2);
        end for;
    end for;

    return witnesses;
end function;

procedure WriteM2248Witness(out, W)
    fprintf out,
        "source=%o | ordered=%o | perm=%o | eps=[%o,%o,%o] | full=%o | rho=%o sigma=%o tau=%o | intermediate=%o | F=[%o,%o,%o,%o] | q=[%o,%o,%o,%o]\n",
        W`source_tuple,
        W`tuple,
        W`permutation,
        W`eps_rho,
        W`eps_sigma,
        W`eps_tau,
        W`full_pass,
        W`rho,
        W`sigma,
        W`tau,
        W`intermediate,
        W`F1,
        W`F2,
        W`F3,
        W`F4,
        W`q1,
        W`q2,
        W`q3,
        W`q4;
end procedure;

procedure SieveM2248TupleFile(input_file, output_file, full_check, max_rows, progress_interval)
    tuples := ReadM2248TupleFile(input_file);
    if max_rows gt 0 and max_rows lt #tuples then
        tuples := tuples[1..max_rows];
    end if;

    out := Open(output_file, "w");

    total_witnesses := 0;
    total_tuples := 0;

    print "Read", #tuples, "tuples from", input_file;
    print "Writing candidates to", output_file;
    print "Full check:", full_check;

    for i in [1..#tuples] do
        witnesses := M2248WitnessesForTuple(tuples[i], full_check);

        if #witnesses gt 0 then
            total_tuples +:= 1;
            total_witnesses +:= #witnesses;
            for W in witnesses do
                WriteM2248Witness(out, W);
            end for;
        end if;

        if progress_interval gt 0 and i mod progress_interval eq 0 then
            print "Processed", i, "of", #tuples,
                  "tuples; candidate tuples:", total_tuples,
                  "witnesses:", total_witnesses;
        end if;
    end for;

    delete out;

    print "Done.";
    print "Candidate tuples:", total_tuples;
    print "Witnesses:", total_witnesses;
end procedure;

procedure SieveM2248TupleFileAllPermutations(input_file, output_file, full_check, max_rows, progress_interval)
    tuples := ReadM2248TupleFile(input_file);
    if max_rows gt 0 and max_rows lt #tuples then
        tuples := tuples[1..max_rows];
    end if;

    out := Open(output_file, "w");

    total_witnesses := 0;
    total_tuples := 0;

    print "Read", #tuples, "tuples from", input_file;
    print "Writing all-permutation candidates to", output_file;
    print "Full check:", full_check;

    for i in [1..#tuples] do
        witnesses := M2248WitnessesForTupleAllPermutations(tuples[i], full_check);

        if #witnesses gt 0 then
            total_tuples +:= 1;
            total_witnesses +:= #witnesses;
            for W in witnesses do
                WriteM2248Witness(out, W);
            end for;
        end if;

        if progress_interval gt 0 and i mod progress_interval eq 0 then
            print "Processed", i, "of", #tuples,
                  "source tuples; candidate source tuples:", total_tuples,
                  "witnesses:", total_witnesses;
        end if;
    end for;

    delete out;

    print "Done.";
    print "Candidate source tuples:", total_tuples;
    print "Witnesses:", total_witnesses;
end procedure;

if assigned input_file then
    if not assigned output_file then
        output_file := "m2248_candidates.txt";
    end if;
    if not assigned full_check then
        full_check := true;
    end if;
    if not assigned max_rows then
        max_rows := 0;
    end if;
    if not assigned progress_interval then
        progress_interval := 1000;
    end if;
    if not assigned all_permutations then
        all_permutations := false;
    end if;

    if Type(full_check) eq MonStgElt then
        full_check := full_check eq "true" or full_check eq "True" or full_check eq "1";
    end if;
    if Type(max_rows) eq MonStgElt then
        max_rows := StringToInteger(max_rows);
    end if;
    if Type(progress_interval) eq MonStgElt then
        progress_interval := StringToInteger(progress_interval);
    end if;

    if Type(all_permutations) eq MonStgElt then
        all_permutations := all_permutations eq "true" or all_permutations eq "True" or all_permutations eq "1";
    end if;

    if all_permutations then
        SieveM2248TupleFileAllPermutations(
            input_file,
            output_file,
            full_check,
            max_rows,
            progress_interval
        );
    else
        SieveM2248TupleFile(
            input_file,
            output_file,
            full_check,
            max_rows,
            progress_interval
        );
    end if;
end if;
