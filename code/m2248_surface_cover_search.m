//////////////////////////////////////////////////////////////////////
//  Direct search on the K3 surface S = M(2,2,2,8)
//
//      s2(a,b,c,d)^2 = 4 a b c d
//
//  for points that lift to the M(2,2,4,8) cover.  This is the direct
//  finite-cover search, as opposed to restricting to a few rational
//  curves on S.
//
//  A signed tuple [a,b,c,d] on S certifies the (2,2,2,8) structure for
//  the curve y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2).  Since the curve only
//  depends on absolute values and the cover code checks all label
//  permutations/signs, candidates are canonicalized as sorted positive
//  primitive curve tuples.
//
//  Typical runs from torsion_jac:
//
//      magma -b output_file:=data/m2248_surface_full_B100.txt \
//          max_abs:=100 full_check:=true primes:="3,5,7,11,13,17,19" \
//          progress_interval:=10000 code/m2248_surface_cover_search.m
//
//      magma -b output_file:=data/m2248_surface_intermediate_B100.txt \
//          max_abs:=100 full_check:=false code/m2248_surface_cover_search.m
//////////////////////////////////////////////////////////////////////

load "code/m2248_sieve.m";

Z := Integers();
Q := Rationals();

function GCDNonzero(vals)
    nz := [ Abs(Z!v) : v in vals | v ne 0 ];
    if #nz eq 0 then
        return Z!0;
    end if;
    return GCD(nz);
end function;

function PrimitiveIntegerTupleFromRationals(vals)
    dens := [ Denominator(Q!v) : v in vals ];
    L := LCM(dens);
    ints := [ Z!(L*(Q!v)) : v in vals ];
    g := GCDNonzero(ints);
    if g ne 0 then
        ints := [ n div g : n in ints ];
    end if;
    return ints;
end function;

function SurfaceEquationValue(tup)
    a := Z!tup[1]; b := Z!tup[2]; c := Z!tup[3]; d := Z!tup[4];
    s2 := a*b + a*c + a*d + b*c + b*d + c*d;
    s4 := a*b*c*d;
    return s2^2 - 4*s4;
end function;

function CanonicalCurveTuple(tup)
    vals := Sort([ Abs(Z!x) : x in tup ]);
    return vals;
end function;

function IsUsableCurveTuple(tup)
    if #tup ne 4 then
        return false;
    end if;
    if &or [ x eq 0 : x in tup ] then
        return false;
    end if;
    return #(Set(tup)) eq 4;
end function;

procedure MaybeRecordSurfaceTuple(~seen, ~records, signed_tuple, max_abs)
    if SurfaceEquationValue(signed_tuple) ne 0 then
        return;
    end if;

    curve_tuple := CanonicalCurveTuple(signed_tuple);
    if not IsUsableCurveTuple(curve_tuple) then
        return;
    end if;
    if Max(curve_tuple) gt max_abs then
        return;
    end if;

    key := Sprint(curve_tuple);
    if key notin seen then
        Include(~seen, key);
        Append(~records, <signed_tuple, curve_tuple>);
    end if;
end procedure;

function EnumerateSurfaceCurveTuples(max_abs)
    records := [];
    seen := {};

    for a in [-max_abs..max_abs] do
        if a eq 0 then
            continue;
        end if;
        for b in [-max_abs..max_abs] do
            if b eq 0 then
                continue;
            end if;
            for c in [-max_abs..max_abs] do
                if c eq 0 then
                    continue;
                end if;

                u := a*b + a*c + b*c;
                v := a + b + c;
                abc := a*b*c;

                if v eq 0 then
                    den := 4*abc;
                    if den eq 0 then
                        continue;
                    end if;
                    d := (Q!u^2)/(Q!den);
                    signed_tuple := PrimitiveIntegerTupleFromRationals([Q!a, Q!b, Q!c, d]);
                    MaybeRecordSurfaceTuple(~seen, ~records, signed_tuple, max_abs);
                else
                    disc := 16*abc*(abc - u*v);
                    ok_disc, sqrt_disc := IsSquare(disc);
                    if not ok_disc then
                        continue;
                    end if;

                    for eps in [-1, 1] do
                        d := (Q!(4*abc - 2*u*v + eps*sqrt_disc))/(Q!(2*v^2));
                        signed_tuple := PrimitiveIntegerTupleFromRationals([Q!a, Q!b, Q!c, d]);
                        MaybeRecordSurfaceTuple(~seen, ~records, signed_tuple, max_abs);
                    end for;
                end if;
            end for;
        end for;
    end for;

    return records;
end function;

function SquareRootsInField(x)
    F := Parent(x);
    if x eq F!0 then
        return [ F!0 ];
    end if;

    ok, r := IsSquare(x);
    if not ok then
        return [];
    end if;
    if r eq -r then
        return [ r ];
    end if;
    return [ r, -r ];
end function;

function A2244LocalPossibleAtPrime(tup, p)
    F := GF(p);

    for perm in M2248PermutationIndices4 do
        a := F!tup[perm[1]]; b := F!tup[perm[2]];
        c := F!tup[perm[3]]; d := F!tup[perm[4]];
        A := a^2; B := b^2; C := c^2; D := d^2;

        conditions := [
            (A-C)*(A-D),
            (B-C)*(B-D),
            (A-C)*(B-C),
            (A-D)*(B-D)
        ];

        if &and [ IsSquare(x) : x in conditions ] then
            return true;
        end if;
    end for;

    return false;
end function;

function A2244LocalPossible(tup, primes)
    return &and [ A2244LocalPossibleAtPrime(tup, p) : p in primes ];
end function;

function CoverLocalPossibleAtPrime(tup, p, full_check)
    F := GF(p);

    // Treat bad-reduction primes as inconclusive, not obstructing.
    if &or [ F!x eq 0 : x in tup ] then
        return true;
    end if;

    for perm in M2248PermutationIndices4 do
        a := F!tup[perm[1]]; b := F!tup[perm[2]];
        c := F!tup[perm[3]]; d := F!tup[perm[4]];
        A := a^2; B := b^2; C := c^2; D := d^2;

        if #(Set([A,B,C,D])) ne 4 then
            // Bad reduction for this labeling; do not let it falsely obstruct.
            return true;
        end if;
        if b eq 0 then
            return true;
        end if;

        y1_roots := SquareRootsInField((A-C)*(A-D));
        y2_roots := SquareRootsInField((B-C)*(B-D));
        w_roots := SquareRootsInField((A-C)*(B-C));
        t_roots := SquareRootsInField((A-D)*(B-D));
        if #y1_roots eq 0 or #y2_roots eq 0 or #w_roots eq 0 or #t_roots eq 0 then
            continue;
        end if;

        for w0 in w_roots do
            if w0 eq 0 then
                return true;
            end if;
            for t0 in t_roots do
                if t0 eq 0 then
                    return true;
                end if;

                sigma := (A-C)/w0;
                tau := (A-D)/t0;

                for eps_rho in [F!-1, F!1] do
                    rho := eps_rho*a/b;

                    intermediate := rho*sigma*tau;
                    if not IsSquare(intermediate) then
                        continue;
                    end if;

                    if not full_check then
                        return true;
                    end if;

                    F1 := (1 + rho)*(1 + sigma)*(1 + tau);
                    F2 := rho*(1 + rho)*(rho + sigma)*(rho + tau);
                    F3 := sigma*(1 + sigma)*(rho + sigma)*(sigma + tau);
                    F4 := tau*(1 + tau)*(rho + tau)*(sigma + tau);
                    if IsSquare(F1) and IsSquare(F2) and IsSquare(F3) and IsSquare(F4) then
                        return true;
                    end if;
                end for;
            end for;
        end for;
    end for;

    return false;
end function;

function CoverLocalPossible(tup, primes, full_check)
    return &and [ CoverLocalPossibleAtPrime(tup, p, full_check) : p in primes ];
end function;

function FirstA2244ObstructionPrime(tup, primes)
    for p in primes do
        if not A2244LocalPossibleAtPrime(tup, p) then
            return p;
        end if;
    end for;
    return 0;
end function;

function FirstCoverObstructionPrime(tup, primes, full_check)
    for p in primes do
        if not CoverLocalPossibleAtPrime(tup, p, full_check) then
            return p;
        end if;
    end for;
    return 0;
end function;

function ParsePrimeList(prime_string)
    if Type(prime_string) ne MonStgElt then
        return prime_string;
    end if;

    parts := Split(prime_string, ",");
    return [ StringToInteger(part) : part in parts | #part gt 0 ];
end function;

procedure WriteSurfaceWitness(out, signed_tuple, curve_tuple, W)
    fprintf out,
        "signed=%o | curve=%o | ordered=%o | perm=%o | eps=[%o,%o,%o] | full=%o | intermediate=%o | F=[%o,%o,%o,%o] | q=[%o,%o,%o,%o]\n",
        signed_tuple,
        curve_tuple,
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

procedure SearchSurfaceCover(output_file, max_abs, full_check, primes, progress_interval)
    out := Open(output_file, "w");

    print "Enumerating primitive curve tuples from S with height <=", max_abs;
    records := EnumerateSurfaceCurveTuples(max_abs);
    print "Unique curve tuples on S:", #records;
    print "Full check:", full_check;
    print "Local primes:", primes;
    print "Output:", output_file;

    local_a2244 := 0;
    local_cover := 0;
    exact_sources := 0;
    exact_witnesses := 0;

    for idx in [1..#records] do
        signed_tuple := records[idx][1];
        curve_tuple := records[idx][2];

        a2244_obstruction := FirstA2244ObstructionPrime(curve_tuple, primes);
        a2244_ok := a2244_obstruction eq 0;
        cover_obstruction := 0;
        cover_ok := false;
        witnesses := [];

        if a2244_ok then
            local_a2244 +:= 1;
            cover_obstruction := FirstCoverObstructionPrime(curve_tuple, primes, full_check);
            cover_ok := cover_obstruction eq 0;
        end if;

        if cover_ok then
            local_cover +:= 1;
            witnesses := M2248WitnessesForTupleAllPermutations(curve_tuple, full_check);
            if #witnesses gt 0 then
                exact_sources +:= 1;
                exact_witnesses +:= #witnesses;
            end if;
        end if;

        fprintf out,
            "CANDIDATE | signed=%o | curve=%o | local_A2244=%o | A2244_obstruction=%o | local_cover=%o | cover_obstruction=%o | exact_witnesses=%o\n",
            signed_tuple,
            curve_tuple,
            a2244_ok,
            a2244_obstruction,
            cover_ok,
            cover_obstruction,
            #witnesses;

        for W in witnesses do
            WriteSurfaceWitness(out, signed_tuple, curve_tuple, W);
        end for;

        if progress_interval gt 0 and idx mod progress_interval eq 0 then
            print "Processed", idx, "of", #records,
                  "local A2244:", local_a2244,
                  "local cover:", local_cover,
                  "exact sources:", exact_sources,
                  "witnesses:", exact_witnesses;
        end if;
    end for;

    fprintf out,
        "SUMMARY | max_abs=%o | full_check=%o | primes=%o | surface_tuples=%o | local_A2244=%o | local_cover=%o | exact_sources=%o | exact_witnesses=%o\n",
        max_abs,
        full_check,
        primes,
        #records,
        local_a2244,
        local_cover,
        exact_sources,
        exact_witnesses;

    delete out;

    print "Surface tuples:", #records;
    print "After local A2244 sieve:", local_a2244;
    print "After local cover sieve:", local_cover;
    print "Exact sources:", exact_sources;
    print "Exact witnesses:", exact_witnesses;
    print "Done.";
end procedure;


procedure SieveSurfaceTupleFile(surface_input_file, output_file, full_check, primes, progress_interval)
    tuples := ReadM2248TupleFile(surface_input_file);
    out := Open(output_file, "w");

    print "Reading surface curve tuples from", surface_input_file;
    print "Tuples:", #tuples;
    print "Full check:", full_check;
    print "Local primes:", primes;
    print "Output:", output_file;

    local_a2244 := 0;
    local_cover := 0;
    exact_sources := 0;
    exact_witnesses := 0;

    for idx in [1..#tuples] do
        curve_tuple := tuples[idx];

        a2244_obstruction := FirstA2244ObstructionPrime(curve_tuple, primes);
        a2244_ok := a2244_obstruction eq 0;
        cover_obstruction := 0;
        cover_ok := false;
        witnesses := [];

        if a2244_ok then
            local_a2244 +:= 1;
            cover_obstruction := FirstCoverObstructionPrime(curve_tuple, primes, full_check);
            cover_ok := cover_obstruction eq 0;
        end if;

        if cover_ok then
            local_cover +:= 1;
            witnesses := M2248WitnessesForTupleAllPermutations(curve_tuple, full_check);
            if #witnesses gt 0 then
                exact_sources +:= 1;
                exact_witnesses +:= #witnesses;
            end if;
        end if;

        fprintf out,
            "CANDIDATE | curve=%o | local_A2244=%o | A2244_obstruction=%o | local_cover=%o | cover_obstruction=%o | exact_witnesses=%o\n",
            curve_tuple,
            a2244_ok,
            a2244_obstruction,
            cover_ok,
            cover_obstruction,
            #witnesses;

        for W in witnesses do
            WriteSurfaceWitness(out, curve_tuple, curve_tuple, W);
        end for;

        if progress_interval gt 0 and idx mod progress_interval eq 0 then
            print "Processed", idx, "of", #tuples,
                  "local A2244:", local_a2244,
                  "local cover:", local_cover,
                  "exact sources:", exact_sources,
                  "witnesses:", exact_witnesses;
        end if;
    end for;

    fprintf out,
        "SUMMARY | input_file=%o | full_check=%o | primes=%o | surface_tuples=%o | local_A2244=%o | local_cover=%o | exact_sources=%o | exact_witnesses=%o\n",
        surface_input_file,
        full_check,
        primes,
        #tuples,
        local_a2244,
        local_cover,
        exact_sources,
        exact_witnesses;

    delete out;

    print "Surface tuples:", #tuples;
    print "After local A2244 sieve:", local_a2244;
    print "After local cover sieve:", local_cover;
    print "Exact sources:", exact_sources;
    print "Exact witnesses:", exact_witnesses;
    print "Done.";
end procedure;

if assigned output_file then
    if not assigned max_abs then
        max_abs := 100;
    end if;
    if not assigned full_check then
        full_check := true;
    end if;
    if not assigned primes then
        primes := [3,5,7,11,13,17,19,23,29,31];
    else
        primes := ParsePrimeList(primes);
    end if;
    if not assigned progress_interval then
        progress_interval := 10000;
    end if;

    if Type(max_abs) eq MonStgElt then
        max_abs := StringToInteger(max_abs);
    end if;
    if Type(full_check) eq MonStgElt then
        full_check := full_check eq "true" or full_check eq "True" or full_check eq "1";
    end if;
    if Type(progress_interval) eq MonStgElt then
        progress_interval := StringToInteger(progress_interval);
    end if;

    if assigned surface_input_file then
        SieveSurfaceTupleFile(surface_input_file, output_file, full_check, primes, progress_interval);
    else
        SearchSurfaceCover(output_file, max_abs, full_check, primes, progress_interval);
    end if;
end if;
