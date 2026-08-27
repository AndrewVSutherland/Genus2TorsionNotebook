//////////////////////////////////////////////////////////////////////
//  Check Howe--Poonen--Leprevost split (2,2,4,8) candidates against
//  the M(2,2,4,8) cover equations after all choices of the two
//  Weierstrass points moved to 0 and infinity.
//
//  This addresses the TODO in paper/NotesAndTodo.tex: a raw tuple
//  [a,b,c,d] for y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2) need not already
//  be in the normalized M2248 chart.  One must test all ordered choices
//  of zero and infinity among the six branch points.
//
//  Typical run from torsion_jac:
//      magma code/m2248_hpl_normalization_check.m
//////////////////////////////////////////////////////////////////////

load "code/m2248_sieve.m";

Q := Rationals();

function FiniteBranchRoots(vals, use_squares)
    roots := [ Q!0 ];
    for v in vals do
        if use_squares then
            Append(~roots, -Q!v^2);
        else
            Append(~roots, -Q!v);
        end if;
    end for;
    return roots;
end function;

function BranchLabel(i)
    if i eq 1 then
        return "0";
    elif i eq 2 then
        return "a";
    elif i eq 3 then
        return "b";
    elif i eq 4 then
        return "c";
    elif i eq 5 then
        return "d";
    elif i eq 6 then
        return "inf";
    end if;
    return "?";
end function;

function MobiusValue(idx, zero_idx, inf_idx, roots)
    // Branch index 6 denotes the original infinity.
    if idx eq zero_idx or idx eq inf_idx then
        error "MobiusValue called on zero or infinity branch.";
    end if;

    if zero_idx eq 6 then
        // T(x) = 1/(x-w), sending original infinity to 0 and w to infinity.
        w := roots[inf_idx];
        return Q!1/(roots[idx] - w);
    elif inf_idx eq 6 then
        // T(x) = x-z, sending z to 0 and original infinity to infinity.
        z := roots[zero_idx];
        return roots[idx] - z;
    else
        z := roots[zero_idx];
        w := roots[inf_idx];
        if idx eq 6 then
            // Limit of (x-z)/(x-w) at original infinity.
            return Q!1;
        end if;
        return (roots[idx] - z)/(roots[idx] - w);
    end if;
end function;

function NormalizedSquareTuple(vals, use_squares, zero_idx, inf_idx)
    roots := FiniteBranchRoots(vals, use_squares);
    remaining := [ i : i in [1..6] | i ne zero_idx and i ne inf_idx ];
    lambdas := [ MobiusValue(i, zero_idx, inf_idx, roots) : i in remaining ];

    if &or [ l eq 0 : l in lambdas ] then
        return false, [], [];
    end if;

    base := lambdas[1];
    tuple := [ Q!1 ];
    for j in [2..#lambdas] do
        ratio := lambdas[j]/base;
        ok, rt := IsSquare(ratio);
        if not ok then
            return false, [], lambdas;
        end if;
        Append(~tuple, rt);
    end for;

    return true, tuple, lambdas;
end function;

procedure CheckCandidate(name, vals, use_squares)
    square_models := 0;
    intermediate_sources := 0;
    full_sources := 0;
    intermediate_witnesses := 0;
    full_witnesses := 0;

    print "CANDIDATE", name, "use_squares", use_squares;

    for zero_idx in [1..6] do
        for inf_idx in [1..6] do
            if zero_idx eq inf_idx then
                continue;
            end if;

            ok, tuple, lambdas := NormalizedSquareTuple(vals, use_squares, zero_idx, inf_idx);
            if not ok then
                continue;
            end if;
            square_models +:= 1;

            iw := M2248WitnessesForTupleAllPermutations(tuple, false);
            fw := M2248WitnessesForTupleAllPermutations(tuple, true);
            intermediate_witnesses +:= #iw;
            full_witnesses +:= #fw;
            if #iw gt 0 then
                intermediate_sources +:= 1;
            end if;
            if #fw gt 0 then
                full_sources +:= 1;
            end if;

            print "NORMALIZED", "zero", BranchLabel(zero_idx),
                  "infinity", BranchLabel(inf_idx),
                  "tuple", tuple,
                  "intermediate", #iw,
                  "full", #fw;
            for W in fw do
                print "  FULL", "ordered", W`tuple, "perm", W`permutation,
                      "eps", [W`eps_rho, W`eps_sigma, W`eps_tau],
                      "rho", W`rho, "sigma", W`sigma, "tau", W`tau;
            end for;
        end for;
    end for;

    print "SUMMARY", name, "use_squares", use_squares,
          "square_normalizations", square_models,
          "intermediate_sources", intermediate_sources,
          "intermediate_witnesses", intermediate_witnesses,
          "full_sources", full_sources,
          "full_witnesses", full_witnesses;
end procedure;

remark_vals := [
    20615879288622173904,
    18727056172448817625,
    17984089341487641600,
    16336390342285800000
];

meeting_vals := [
    124204233895590989520,
    88395624606340961712,
    70848937725658001712,
    64357770918797390375
];

CheckCandidate("remark_1091", remark_vals, true);
CheckCandidate("remark_1091", remark_vals, false);
CheckCandidate("meeting_294", meeting_vals, true);
CheckCandidate("meeting_294", meeting_vals, false);
