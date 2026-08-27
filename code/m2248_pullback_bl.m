//////////////////////////////////////////////////////////////////////
//  Pull back the M(2,2,4,8) halving conditions to explicit rational
//  curves on the M_1(8,2,2,2) K3 surface.
//
//  This script depends on code/m2248_sieve.m.
//
//  It handles two sources of rational curves:
//    (1) the explicit M(2,2,2,8) section families in main.tex;
//    (2) the Bertin--Lecacheux section from NotesAndTodo.tex, converted
//        from the (u,v,q) branch model to normalized [a,b,c,d] tuples
//        when the quadratic q splits rationally.
//
//  Typical command-line runs from the torsion_jac directory:
//
//      magma -b output_file:=data/m2248_pullback_intermediate_h20.txt \
//          height:=20 full_check:=false progress_interval:=100 \
//          code/m2248_pullback_bl.m
//
//      magma -b output_file:=data/m2248_pullback_full_h40.txt \
//          height:=40 full_check:=true progress_interval:=250 \
//          code/m2248_pullback_bl.m
//
//      magma -b symbolic_output_file:=data/m2248_section_symbolic.txt \
//          full_check:=true code/m2248_pullback_bl.m
//
//      magma -b squareclass_output_file:=data/m2248_section_squareclasses.txt \
//          code/m2248_pullback_bl.m
//////////////////////////////////////////////////////////////////////

load "code/m2248_sieve.m";

Q := Rationals();
P<x> := PolynomialRing(Q);

M2248SectionFamilyNames := [
    "k3_section_P",
    "k3_section_P_plus_T1",
    "k3_section_P_plus_T2",
    "k3_section_P_plus_T3",
    "filip_projective"
];

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

function PrimitiveIntegerTuple(vals)
    dens := [ Denominator(v) : v in vals ];
    L := LCM(dens);
    ints := [ Integers()!(L*v) : v in vals ];
    g := GCD([ Abs(n) : n in ints | n ne 0 ]);
    if g ne 0 then
        ints := [ n div g : n in ints ];
    end if;

    return ints;
end function;

function SectionTuple(family, t)
    K := Parent(t);
    one := K!1;

    if family eq "k3_section_P" then
        return [
            -4*t^2*(t+one)/(t^2+t+one)^2,
            -t/(t+one),
            one,
            t
        ];
    elif family eq "k3_section_P_plus_T1" then
        return [
            -(t^2+t+one)^2/(4*(t^2+t)),
            -(t+one),
            one,
            t
        ];
    elif family eq "k3_section_P_plus_T2" then
        return [
            -(t^4 - 2*t^3 - t^2 + 2*t + one)/
                (t^4 + 2*t^3 - t^2 - 2*t + one),
            -1/t,
            one,
            t
        ];
    elif family eq "k3_section_P_plus_T3" then
        return [
            -t*(t+one)^2*(t-one)/((t^2 - 2*t - one)*(t^2 + one)),
            -t^2,
            one,
            t
        ];
    elif family eq "filip_projective" then
        return [
            -(t^2+t+one)^2,
            -4*t*(t+one)^2,
            4*t*(t+one),
            4*t^2*(t+one)
        ];
    else
        error "Unknown section family.";
    end if;
end function;

function IsUsableTuple(vals)
    if #vals ne 4 then
        return false;
    end if;
    if &or [ v eq 0 : v in vals ] then
        return false;
    end if;

    sqs := [ v^2 : v in vals ];
    return #(Set(sqs)) eq 4;
end function;

function BLuv(s)
    u := ((2*s^2 - 2*s + 1)*(s^3 - s^2 + 2*s - 1)) /
         ((s^2 - 1)*(s^3 - 2*s^2 + s - 1));
    v := s*(s^2 - 2*s + 2)/(s^2 - 1);
    return u, v;
end function;

function BLQuadratic(u, v)
    return
        -(u^2 + u*v + v^2 + u + v + 1)*x^2
        + (v+1)*(u+1)*(u+v)*x
        - u*v*(u+v+1);
end function;

function HasDistinctEntries(vals)
    return #(Set(vals)) eq #vals;
end function;

function BLNormalizedTuples(s)
    u, v := BLuv(s);
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


function M2248SymbolicCountsForTuple(vals, full_check)
    if #vals ne 4 then
        return 0, 0, 0, [];
    end if;

    K := Parent(vals[1]);
    if &or [ Parent(v) ne K : v in vals ] then
        vals := [ K!v : v in vals ];
    end if;
    if &or [ v eq 0 : v in vals ] then
        return 0, 0, 0, [];
    end if;

    a2244_count := 0;
    intermediate_count := 0;
    full_count := 0;
    witnesses := [];

    for perm in M2248PermutationIndices4 do
        a := vals[perm[1]];
        b := vals[perm[2]];
        c := vals[perm[3]];
        d := vals[perm[4]];

        A := a^2; B := b^2; C := c^2; D := d^2;
        if #(Set([A,B,C,D])) ne 4 then
            continue;
        end if;

        ok_y1, y10 := IsSquare((A-C)*(A-D));
        ok_y2, y20 := IsSquare((B-C)*(B-D));
        ok_w, w0 := IsSquare((A-C)*(B-C));
        ok_t, t0 := IsSquare((A-D)*(B-D));

        if not (ok_y1 and ok_y2 and ok_w and ok_t) then
            continue;
        end if;

        a2244_count +:= 1;

        for eps_rho in [-1, 1] do
            for eps_sigma in [-1, 1] do
                for eps_tau in [-1, 1] do
                    rho := (K!eps_rho) * a / b;
                    sigma := (K!eps_sigma) * (A-C) / w0;
                    tau := (K!eps_tau) * (A-D) / t0;

                    intermediate := rho * sigma * tau;
                    ok0, r0 := IsSquare(intermediate);
                    if not ok0 then
                        continue;
                    end if;

                    intermediate_count +:= 1;

                    F1 := (K!1 + rho) * (K!1 + sigma) * (K!1 + tau);
                    F2 := rho * (K!1 + rho) * (rho + sigma) * (rho + tau);
                    F3 := sigma * (K!1 + sigma) * (rho + sigma) * (sigma + tau);
                    F4 := tau * (K!1 + tau) * (rho + tau) * (sigma + tau);

                    ok1, q1 := IsSquare(F1);
                    ok2, q2 := IsSquare(F2);
                    ok3, q3 := IsSquare(F3);
                    ok4, q4 := IsSquare(F4);
                    full_pass := ok1 and ok2 and ok3 and ok4;

                    if full_pass then
                        full_count +:= 1;
                    end if;
                    if full_check and not full_pass then
                        continue;
                    end if;

                    Append(~witnesses, <
                        perm,
                        eps_rho,
                        eps_sigma,
                        eps_tau,
                        full_pass,
                        intermediate,
                        F1,
                        F2,
                        F3,
                        F4
                    >);
                end for;
            end for;
        end for;
    end for;

    return a2244_count, intermediate_count, full_count, witnesses;
end function;

procedure AnalyzeSectionPullbacksSymbolic(output_file, full_check)
    K<T> := RationalFunctionField(Q);
    out := Open(output_file, "w");

    total_a2244 := 0;
    total_intermediate := 0;
    total_full := 0;
    total_reported := 0;

    print "Symbolic section-family pullback analysis over Q(t)";
    print "Full check:", full_check;
    print "Output:", output_file;

    for family in M2248SectionFamilyNames do
        vals := SectionTuple(family, T);
        a2244_count, intermediate_count, full_count, witnesses :=
            M2248SymbolicCountsForTuple(vals, full_check);

        total_a2244 +:= a2244_count;
        total_intermediate +:= intermediate_count;
        total_full +:= full_count;
        total_reported +:= #witnesses;

        print "Symbolic family", family,
              "A2244 orderings:", a2244_count,
              "intermediate signs:", intermediate_count,
              "full signs:", full_count,
              "reported:", #witnesses;

        fprintf out,
            "family=%o | A2244_orderings=%o | intermediate_signs=%o | full_signs=%o | reported=%o\n",
            family,
            a2244_count,
            intermediate_count,
            full_count,
            #witnesses;

        for W in witnesses do
            fprintf out,
                "  perm=%o | eps=[%o,%o,%o] | full=%o | intermediate=%o | F=[%o,%o,%o,%o]\n",
                W[1],
                W[2],
                W[3],
                W[4],
                W[5],
                W[6],
                W[7],
                W[8],
                W[9],
                W[10];
        end for;
    end for;

    fprintf out,
        "TOTAL | A2244_orderings=%o | intermediate_signs=%o | full_signs=%o | reported=%o\n",
        total_a2244,
        total_intermediate,
        total_full,
        total_reported;

    delete out;

    print "Symbolic totals - A2244 orderings:", total_a2244,
          "intermediate signs:", total_intermediate,
          "full signs:", total_full,
          "reported:", total_reported;
end procedure;


function SquarefreeRationalConstant(c)
    c := Q!c;
    if c eq 0 then
        return Q!0;
    end if;

    out := c lt 0 select Q!-1 else Q!1;
    num := Abs(Numerator(c));
    den := Denominator(c);

    for fac in Factorization(num) do
        if fac[2] mod 2 eq 1 then
            out *:= Q!fac[1];
        end if;
    end for;
    for fac in Factorization(den) do
        if fac[2] mod 2 eq 1 then
            out /:= Q!fac[1];
        end if;
    end for;

    return out;
end function;

function SquareclassRepresentative(f)
    num := Numerator(f);
    den := Denominator(f);
    poly := num * den;
    R := Parent(poly);

    if poly eq 0 then
        return R!0;
    end if;

    rep := R!SquarefreeRationalConstant(LeadingCoefficient(poly));
    for fac in Factorization(poly) do
        if fac[2] mod 2 eq 1 then
            rep *:= fac[1];
        end if;
    end for;

    return rep;
end function;

function A2244SquareclassesForOrderedTuple(vals)
    a := vals[1];
    b := vals[2];
    c := vals[3];
    d := vals[4];

    A := a^2; B := b^2; C := c^2; D := d^2;
    conditions := [
        (A-C)*(A-D),
        (B-C)*(B-D),
        (A-C)*(B-C),
        (A-D)*(B-D)
    ];

    return [ SquareclassRepresentative(f) : f in conditions ];
end function;

procedure AnalyzeSectionA2244Squareclasses(output_file)
    K<T> := RationalFunctionField(Q);
    out := Open(output_file, "w");

    print "Writing section-family A2244 squareclass diagnostics to", output_file;

    for family in M2248SectionFamilyNames do
        vals := SectionTuple(family, T);
        for perm in M2248PermutationIndices4 do
            ordered := [ vals[i] : i in perm ];
            squareclasses := A2244SquareclassesForOrderedTuple(ordered);
            degrees := [ Degree(sc) : sc in squareclasses ];
            total_degree := &+degrees;
            max_degree := Max(degrees);
            distinct_classes := #Set([ Sprint(sc) : sc in squareclasses ]);
            bad_constant := false;
            for sc in squareclasses do
                if Degree(sc) eq 0 then
                    ok_const, root_const := IsSquare(Coefficient(sc, 0));
                    if not ok_const then
                        bad_constant := true;
                    end if;
                end if;
            end for;

            fprintf out,
                "total_degree=%o | max_degree=%o | distinct_classes=%o | bad_constant=%o | family=%o | perm=%o | degrees=%o | squareclasses=%o\n",
                total_degree,
                max_degree,
                distinct_classes,
                bad_constant,
                family,
                perm,
                degrees,
                squareclasses;
        end for;
    end for;

    delete out;
    print "Done writing squareclass diagnostics.";
end procedure;

procedure WritePullbackWitness(out, source, param, raw_tuple, W)
    fprintf out,
        "family=%o | param=%o | raw=%o | source=%o | ordered=%o | perm=%o | eps=[%o,%o,%o] | full=%o | intermediate=%o | F=[%o,%o,%o,%o] | q=[%o,%o,%o,%o]\n",
        source,
        param,
        raw_tuple,
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
        W`F4,
        W`q1,
        W`q2,
        W`q3,
        W`q4;
end procedure;

procedure SearchSectionPullbacks(out, params, full_check, progress_interval)
    total_specs := 0;
    total_sources := 0;
    total_witnesses := 0;

    for family in M2248SectionFamilyNames do
        fam_sources := 0;
        fam_witnesses := 0;

        for idx in [1..#params] do
            t := params[idx];
            total_specs +:= 1;

            ok := false;
            vals := [];
            try
                vals := SectionTuple(family, t);
                ok := IsUsableTuple(vals);
            catch e
                ok := false;
            end try;

            if not ok then
                continue;
            end if;

            tuple := PrimitiveIntegerTuple(vals);
            witnesses := M2248WitnessesForTupleAllPermutations(tuple, full_check);
            if #witnesses eq 0 then
                continue;
            end if;

            total_sources +:= 1;
            total_witnesses +:= #witnesses;
            fam_sources +:= 1;
            fam_witnesses +:= #witnesses;

            for W in witnesses do
                WritePullbackWitness(out, family, t, tuple, W);
            end for;

            if progress_interval gt 0 and total_specs mod progress_interval eq 0 then
                print "Processed", total_specs, "section specializations;",
                      "sources:", total_sources,
                      "witnesses:", total_witnesses;
            end if;
        end for;

        print "Section family", family, "candidate sources:", fam_sources,
              "witnesses:", fam_witnesses;
    end for;

    print "Section pullback total sources:", total_sources;
    print "Section pullback total witnesses:", total_witnesses;
end procedure;

procedure SearchBLPullbacks(out, params, full_check, progress_interval)
    total_sources := 0;
    total_witnesses := 0;
    total_tuples := 0;

    for idx in [1..#params] do
        s := params[idx];
        tuples := [];

        try
            tuples := BLNormalizedTuples(s);
        catch e
            tuples := [];
        end try;

        total_tuples +:= #tuples;

        for data in tuples do
            branch_choice := data[1];
            tuple_rat := data[2];
            tuple := PrimitiveIntegerTuple(tuple_rat);

            witnesses := M2248WitnessesForTupleAllPermutations(tuple, full_check);
            if #witnesses eq 0 then
                continue;
            end if;

            total_sources +:= 1;
            total_witnesses +:= #witnesses;

            raw := <branch_choice, tuple>;
            for W in witnesses do
                WritePullbackWitness(out, "bertin_lecacheux_section", s, raw, W);
            end for;
        end for;

        if progress_interval gt 0 and idx mod progress_interval eq 0 then
            print "Processed", idx, "BL parameters;",
                  "normalized tuples:", total_tuples,
                  "sources:", total_sources,
                  "witnesses:", total_witnesses;
        end if;
    end for;

    print "BL normalized tuples:", total_tuples;
    print "BL pullback total sources:", total_sources;
    print "BL pullback total witnesses:", total_witnesses;
end procedure;

procedure SearchM2248Pullbacks(output_file, height, full_check, progress_interval)
    params := RationalParametersOfHeight(height);
    out := Open(output_file, "w");

    print "Searching", #params, "rational parameters of height <=", height;
    print "Full check:", full_check;
    print "Output:", output_file;

    SearchSectionPullbacks(out, params, full_check, progress_interval);
    SearchBLPullbacks(out, params, full_check, progress_interval);

    delete out;
    print "Done.";
end procedure;

if assigned squareclass_output_file then
    AnalyzeSectionA2244Squareclasses(squareclass_output_file);
end if;

if assigned symbolic_output_file then
    if not assigned full_check then
        full_check := false;
    end if;
    if Type(full_check) eq MonStgElt then
        full_check := full_check eq "true" or full_check eq "True" or full_check eq "1";
    end if;

    AnalyzeSectionPullbacksSymbolic(symbolic_output_file, full_check);
end if;

if assigned output_file then
    if not assigned height then
        height := 20;
    end if;
    if not assigned full_check then
        full_check := false;
    end if;
    if not assigned progress_interval then
        progress_interval := 100;
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

    SearchM2248Pullbacks(output_file, height, full_check, progress_interval);
end if;
