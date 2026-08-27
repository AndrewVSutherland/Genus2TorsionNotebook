//////////////////////////////////////////////////////////////////////
//  Z/5 x Z/5 full Mumford norm, linear-B branch b2=0.
//
//  Starting point:
//
//      f = h^2 - K*x^5,              h = 1 + h1*x + h2*x^2,
//      A^2 - B^2*f = U^5,
//      U = x^2 + s*x + t,
//      B = b0 + b1*x                 (b2=0),
//      A = x^5 + a4*x^4 + ... + a0.
//
//  The coefficients a4..a0 are forced by the x^9..x^5 terms.  This
//  script studies the five residual equations on the open chart
//
//      K*disc(f)*disc(U)*Res(B,U)*b0*b1 != 0.
//
//  Typical bounded run from the repository root:
//
//      magma -b prime_bound:=7 count_prime_bound:=7 \
//          do_primary:=true code/agent_z5x5_b2zero_elim.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned prime_bound then
    prime_bound := 7;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;

if not assigned count_prime_bound then
    count_prime_bound := 7;
elif Type(count_prime_bound) eq MonStgElt then
    count_prime_bound := StringToInteger(count_prime_bound);
end if;

if not assigned sample_limit then
    sample_limit := 5;
elif Type(sample_limit) eq MonStgElt then
    sample_limit := StringToInteger(sample_limit);
end if;

if not assigned do_primary then
    do_primary := false;
elif Type(do_primary) eq MonStgElt then
    do_primary := (do_primary eq "true") or (do_primary eq "True");
end if;

if not assigned do_saturation then
    do_saturation := false;
elif Type(do_saturation) eq MonStgElt then
    do_saturation := (do_saturation eq "true") or (do_saturation eq "True");
end if;

if not assigned do_dimension then
    do_dimension := false;
elif Type(do_dimension) eq MonStgElt then
    do_dimension := (do_dimension eq "true") or (do_dimension eq "True");
end if;

if not assigned do_slice_saturation then
    do_slice_saturation := false;
elif Type(do_slice_saturation) eq MonStgElt then
    do_slice_saturation := (do_slice_saturation eq "true") or (do_slice_saturation eq "True");
end if;

function ForcedAData(h, kap, U, B)
    P := Parent(h);
    x := P.1;
    f := h^2 - kap*x^5;
    A := x^5;
    E := A^2 - B^2*f - U^5;
    avec := [];

    for deg in [9,8,7,6,5] do
        a := -Coefficient(E, deg)/2;
        Append(~avec, a);
        A +:= a*x^(deg - 5);
        E := A^2 - B^2*f - U^5;
    end for;

    return A, E, f, avec;
end function;

function FullResidualData(k)
    R<h1,h2,kap,s,t,b0,b1,b2> := PolynomialRing(k, 8, "grevlex");
    Fr := FieldOfFractions(R);
    P<x> := PolynomialRing(Fr);

    h := 1 + h1*x + h2*x^2;
    U := x^2 + s*x + t;
    B := b0 + b1*x + b2*x^2;
    A, E, f, avec := ForcedAData(h, kap, U, B);
    res := [R!Numerator(Coefficient(E, i)) : i in [0..4]];
    return R, res, A, f, U, B;
end function;

function B2ZeroData(k)
    R<h1,h2,kap,s,t,b0,b1> := PolynomialRing(k, 7, "grevlex");
    Fr := FieldOfFractions(R);
    P<x> := PolynomialRing(Fr);

    h := 1 + h1*x + h2*x^2;
    U := x^2 + s*x + t;
    B := b0 + b1*x;
    A, E, f, avec := ForcedAData(h, kap, U, B);
    res := [R!Numerator(Coefficient(E, i)) : i in [0..4]];

    disc_f := R!Numerator(Discriminant(f));
    disc_U := R!Numerator(Discriminant(U));
    res_BU := R!Numerator(Resultant(B, U));
    boundary_factors := [kap, b1, b0, disc_U, res_BU, disc_f];
    boundary_names := ["K", "b1", "b0", "discU", "resBU", "discF"];
    boundary := &*boundary_factors;

    return R, res, boundary_factors, boundary_names, boundary, A, f, U, B;
end function;

function RelationPairs(J, D0, D2)
    rels := [];
    for a in [0..4] do
        for b in [0..4] do
            if a eq 0 and b eq 0 then
                continue;
            end if;
            if a*D0 + b*D2 eq J!0 then
                Append(~rels, <a,b>);
            end if;
        end for;
    end for;
    return rels;
end function;

function CheckClasses(f, h, U, V)
    P := Parent(f);
    x := P.1;
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    D0 := J![x, Evaluate(h, 0)];
    D2 := J![U, V];

    if D0 eq J!0 or D2 eq J!0 then
        return false, [], "zero class";
    end if;
    if 5*D0 ne J!0 then
        return false, [], "5*D0 != 0";
    end if;
    if 5*D2 ne J!0 then
        return false, [], "5*D2 != 0";
    end if;

    rels := RelationPairs(J, D0, D2);
    return (#rels eq 0), rels, "checked";
end function;

function B2ZeroModel(k, vals)
    P<x> := PolynomialRing(k);
    h := 1 + vals[1]*x + vals[2]*x^2;
    kap := vals[3];
    U := x^2 + vals[4]*x + vals[5];
    B := vals[6] + vals[7]*x;
    A, E, f, avec := ForcedAData(h, kap, U, B);

    if &or[Coefficient(E, i) ne 0 : i in [0..4]] then
        return false, A, B, f, U, P!0, [], "residual nonzero";
    end if;
    if kap eq 0 or Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false, A, B, f, U, P!0, [], "bad contact curve";
    end if;
    if vals[7] eq 0 then
        return false, A, B, f, U, P!0, [], "constant B";
    end if;
    if Discriminant(U) eq 0 or GCD(U, f) ne 1 then
        return false, A, B, f, U, P!0, [], "bad U";
    end if;
    if GCD(B, U) ne 1 then
        return false, A, B, f, U, P!0, [], "B not invertible mod U";
    end if;

    Binv := InverseMod(B, U);
    V := (-A*Binv) mod U;
    if (V^2 - f) mod U ne 0 then
        return false, A, B, f, U, V, [], "V recovery failed";
    end if;

    independent, rels, reason := CheckClasses(f, h, U, V);
    if not independent then
        return false, A, B, f, U, V, rels, reason;
    end if;
    return true, A, B, f, U, V, rels, "independent";
end function;

function ResidualJacobianRank(res, vals)
    R := Parent(res[1]);
    F := BaseRing(R);
    entries := [];
    for i in [1..#res] do
        for j in [1..7] do
            Append(~entries, Evaluate(Derivative(res[i], j), vals));
        end for;
    end for;
    M := Matrix(F, #res, 7, entries);
    return Rank(M);
end function;

function HasUnitGenerator(I)
    for g in Basis(I) do
        if g ne 0 and TotalDegree(g) eq 0 then
            return true;
        end if;
    end for;
    return false;
end function;

procedure PrintEquationSizes(label, res)
    printf "%o residual sizes\n", label;
    for i in [1..#res] do
        printf "  E%o degree=%o total_degree=%o terms=%o\n",
            i - 1, Degree(res[i]), TotalDegree(res[i]), #Terms(res[i]);
    end for;
end procedure;

procedure PrintIdealReport(label, I, primary)
    printf "  %o basis_len=%o unit=%o\n", label, #Basis(I), HasUnitGenerator(I);
    if do_dimension then
        try
            dim, component_degrees := Dimension(I);
            printf "  %o dimension=%o component_degrees=%o\n",
                label, dim, component_degrees;
        catch e
            printf "  %o dimension_failed=%o\n", label, e`Object;
        end try;
        try
            printf "  %o degree=%o\n", label, Degree(I);
        catch e
            printf "  %o degree_unavailable=%o\n", label, e`Object;
        end try;
    else
        printf "  %o dimension_skipped\n", label;
    end if;

    if not primary then
        return;
    end if;

    try
        comps := PrimaryDecomposition(I);
        printf "  %o primary_components=%o\n", label, #comps;
        for j in [1..#comps] do
            C := comps[j];
            try
                cdim, cdegrees := Dimension(C);
                printf "    component_%o dimension=%o component_degrees=%o",
                    j, cdim, cdegrees;
            catch e
                printf "    component_%o dimension_failed=%o", j, e`Object;
            end try;
            try
                printf " degree=%o", Degree(C);
            catch e
                printf " degree_unavailable";
            end try;
            printf " basis_degrees=%o\n", [TotalDegree(g) : g in Basis(C)];
        end for;
    catch e
        printf "  %o primary_decomposition_failed=%o\n", label, e`Object;
    end try;
end procedure;

function CountOpenPointsBrute(F, res, boundary_factors, limit)
    elts := [a : a in F];
    raw_count := 0;
    open_count := 0;
    samples := [];

    for h1v in elts do
    for h2v in elts do
    for kapv in elts do
    for sv in elts do
    for tv in elts do
    for b0v in elts do
    for b1v in elts do
        vals := [h1v,h2v,kapv,sv,tv,b0v,b1v];
        residual_zero := true;
        for e in res do
            if Evaluate(e, vals) ne 0 then
                residual_zero := false;
                break;
            end if;
        end for;
        if not residual_zero then
            continue;
        end if;

        raw_count +:= 1;
        open := true;
        for b in boundary_factors do
            if Evaluate(b, vals) eq 0 then
                open := false;
                break;
            end if;
        end for;
        if open then
            open_count +:= 1;
            if #samples lt limit then
                Append(~samples, vals);
            end if;
        end if;
    end for;
    end for;
    end for;
    end for;
    end for;
    end for;
    end for;

    return raw_count, open_count, samples;
end function;

function CountSlicePointsBrute(F, res, boundary_factors, h1fix, h2fix, limit)
    elts := [a : a in F];
    raw_count := 0;
    open_count := 0;
    samples := [];

    for kapv in elts do
    for sv in elts do
    for tv in elts do
    for b0v in elts do
    for b1v in elts do
        vals := [h1fix,h2fix,kapv,sv,tv,b0v,b1v];
        residual_zero := true;
        for e in res do
            if Evaluate(e, vals) ne 0 then
                residual_zero := false;
                break;
            end if;
        end for;
        if not residual_zero then
            continue;
        end if;

        raw_count +:= 1;
        open := true;
        for b in boundary_factors do
            if Evaluate(b, vals) eq 0 then
                open := false;
                break;
            end if;
        end for;
        if open then
            open_count +:= 1;
            if #samples lt limit then
                Append(~samples, vals);
            end if;
        end if;
    end for;
    end for;
    end for;
    end for;
    end for;

    return raw_count, open_count, samples;
end function;

procedure PrintSampleModels(p, F, res, samples)
    Z := Integers();
    printf "  sample_count_printed=%o\n", #samples;
    for i in [1..#samples] do
        vals := samples[i];
        jrank := ResidualJacobianRank(res, vals);
        ok, A, B, f, U, V, rels, reason := B2ZeroModel(F, vals);
        printf "    sample_%o vals=<h1=%o,h2=%o,K=%o,s=%o,t=%o,b0=%o,b1=%o> jacobian_rank=%o class=%o",
            i, Z!vals[1], Z!vals[2], Z!vals[3], Z!vals[4], Z!vals[5],
            Z!vals[6], Z!vals[7], jrank, reason;
        if #rels gt 0 then
            printf " rels=%o", rels;
        end if;
        print "";
        printf "      f=%o\n", f;
        printf "      U=%o\n", U;
        printf "      V=%o\n", V;
        printf "      A=%o\n", A;
        printf "      B=%o\n", B;
        if ok then
            printf "      #J(F_%o)=%o independent=true\n", p,
                #Jacobian(HyperellipticCurve(f));
        end if;
    end for;
end procedure;

procedure SliceReport(p, F, R, res, boundary_factors, samples)
    if #samples eq 0 then
        return;
    end if;

    h1fix := samples[1][1];
    h2fix := samples[1][2];
    printf "  slice_chart h1=%o h2=%o\n", Integers()!h1fix, Integers()!h2fix;
    raw_slice, open_slice, slice_samples :=
        CountSlicePointsBrute(F, res, boundary_factors, h1fix, h2fix, sample_limit);
    printf "  slice_F_%o raw_residual_points=%o open_points=%o\n",
        p, raw_slice, open_slice;
    PrintSampleModels(p, F, res, slice_samples);

    if not do_slice_saturation then
        print "  slice_saturation_skipped pass -b do_slice_saturation:=true";
        return;
    end if;

    S<kapS,sS,tS,b0S,b1S> := PolynomialRing(F, 5, "grevlex");
    phi := hom<R -> S | h1fix, h2fix, kapS, sS, tS, b0S, b1S>;
    sres := [phi(e) : e in res];
    sboundary_factors := [phi(b) : b in boundary_factors];
    sboundary := &*sboundary_factors;
    Islice := ideal<S | sres>;
    PrintIdealReport("slice_raw", Islice, false);
    Jslice := Saturation(Islice, ideal<S | sboundary>);
    PrintIdealReport("slice_sat_product", Jslice, do_primary);
end procedure;

procedure SymbolicReport()
    Q := Rationals();
    Rfull, full_res, Afull, ffull, Ufull, Bfull := FullResidualData(Q);
    PrintEquationSizes("# full b2-variable", full_res);

    Rlin, lin_res, boundary_factors, boundary_names, boundary, A, f, U, B :=
        B2ZeroData(Q);
    PrintEquationSizes("# b2=0", lin_res);
    print "# b2=0 open factors";
    for i in [1..#boundary_factors] do
        printf "  %o degree=%o total_degree=%o terms=%o : %o\n",
            boundary_names[i], Degree(boundary_factors[i]),
            TotalDegree(boundary_factors[i]), #Terms(boundary_factors[i]),
            boundary_factors[i];
    end for;
    printf "  product degree=%o total_degree=%o terms=%o\n",
        Degree(boundary), TotalDegree(boundary), #Terms(boundary);
end procedure;

procedure SaturationPrimeReport(p)
    F := GF(p);
    R, res, boundary_factors, boundary_names, boundary, A, f, U, B :=
        B2ZeroData(F);
    printf "\n# prime p=%o b2=0 modular saturation\n", p;

    if p le count_prime_bound then
        raw_count, open_count, samples :=
            CountOpenPointsBrute(F, res, boundary_factors, sample_limit);
        printf "  brute_F_%o raw_residual_points=%o open_points=%o\n",
            p, raw_count, open_count;
        PrintSampleModels(p, F, res, samples);
        SliceReport(p, F, R, res, boundary_factors, samples);
    else
        printf "  brute_count_skipped p=%o count_prime_bound=%o\n",
            p, count_prime_bound;
    end if;

    if not do_saturation then
        print "  saturation_skipped pass -b do_saturation:=true to compute ideals";
        return;
    end if;

    I := ideal<R | res>;
    PrintIdealReport("raw", I, false);

    J := I;
    for i in [1..#boundary_factors] do
        if boundary_factors[i] eq 0 then
            printf "  skip_saturation_%o factor_is_zero_mod_%o\n",
                boundary_names[i], p;
            continue;
        end if;
        J := Saturation(J, ideal<R | boundary_factors[i]>);
        PrintIdealReport("sat_" cat boundary_names[i], J, false);
        if HasUnitGenerator(J) then
            printf "  saturation became unit after %o\n", boundary_names[i];
            break;
        end if;
    end for;

    Jprod := Saturation(I, ideal<R | boundary>);
    PrintIdealReport("sat_product", Jprod, do_primary);
end procedure;

SymbolicReport();

print "";
print "# modular b2=0 reports";
for p in PrimesUpTo(prime_bound) do
    if p in {2,5} then
        continue;
    end if;
    SaturationPrimeReport(p);
end for;

quit;
