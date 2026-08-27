//////////////////////////////////////////////////////////////////////
//  Focused mod-l^2 tangent diagnostic for traced first-cover points.
//
//  Input files are the compact trace_file outputs produced by
//  agent_m18_416_live_stratum_search.m with trace_firstposs:=true:
//
//      R W st0 st1 gate_R gate_W plus minus tangent
//
//  For each traced rational pair and each requested good aux prime ell,
//  this checks the tangent-candidate congruence used by TangentCandidates
//  modulo ell and modulo ell^2.  The intended use is to tell whether the
//  remaining rational first-cover possibilities already have a local
//  tangent obstruction at a good prime, and if so whether that obstruction
//  appears only one level deeper.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned primes then prime_list := [37,41];
elif Type(primes) eq MonStgElt then
    prime_list := (#primes eq 0) select []
        else [StringToInteger(s) : s in Split(primes, ",")];
else prime_list := primes; end if;

if not assigned trace_files then
    trace_file_list := [
        "data/agent_m18_416_trace_firstposs_aux37_41_h100_part0.txt",
        "data/agent_m18_416_trace_firstposs_aux37_41_h100_part1.txt",
        "data/agent_m18_416_trace_firstposs_aux37_41_h100_part2.txt",
        "data/agent_m18_416_trace_firstposs_aux37_41_h100_part3.txt",
        "data/agent_m18_416_trace_firstposs_aux37_41_h100_part4.txt",
        "data/agent_m18_416_trace_firstposs_aux37_41_h100_part5.txt"
    ];
elif Type(trace_files) eq MonStgElt then
    trace_file_list := [s : s in Split(trace_files, ",") | #s gt 0];
else
    trace_file_list := trace_files;
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

TraceRec := recformat<
    R, W, st0, st1, gateR, gateW, plus, minus, tangent, source
>;

function ReadTraceFiles(files)
    rows := [];
    for filename in files do
        lines := Split(Read(filename), "\n");
        for raw in lines do
            s := raw;
            if #s gt 0 and s[#s] eq "\r" then s := s[1..#s-1]; end if;
            if #s eq 0 or s[1] eq "#" then continue; end if;
            parts := [part : part in Split(s, " ") | #part gt 0];
            if #parts lt 9 then continue; end if;
            Append(~rows, rec<TraceRec |
                R := StringToRational(parts[1]),
                W := StringToRational(parts[2]),
                st0 := StringToInteger(parts[3]),
                st1 := StringToInteger(parts[4]),
                gateR := StringToInteger(parts[5]),
                gateW := StringToInteger(parts[6]),
                plus := StringToInteger(parts[7]),
                minus := StringToInteger(parts[8]),
                tangent := StringToInteger(parts[9]),
                source := filename
            >);
        end for;
    end for;
    return rows;
end function;

function BoundaryLabelsFinite(F, R, w)
    labels := [];
    if R eq 0 then Append(~labels, "R"); end if;
    if w eq 0 then Append(~labels, "w"); end if;
    if w - 1 eq 0 then Append(~labels, "w-1"); end if;
    if w + 1 eq 0 then Append(~labels, "w+1"); end if;
    if R - 1 eq 0 then Append(~labels, "R-1"); end if;
    if R + 1 eq 0 then Append(~labels, "R+1"); end if;
    if R - w eq 0 then Append(~labels, "R-w"); end if;
    if R + w eq 0 then Append(~labels, "R+w"); end if;
    if R*w - 3*R + 3*w - 1 eq 0 then Append(~labels, "Lplus"); end if;
    if R*w + 3*R + 3*w + 1 eq 0 then Append(~labels, "Lminus"); end if;
    if 2*R^2 - R*w^2 + R - 2*w^2 eq 0 then Append(~labels, "Q"); end if;
    if R^4 - 2*R^3 + R^2*w^2 - R^2 + 2*R*w^2 - w^2 eq 0 then
        Append(~labels, "Quartic");
    end if;
    return labels;
end function;

function ResidueInt(q, m)
    den := (Z!Denominator(q)) mod m;
    if GCD(den, m) ne 1 then return false, 0; end if;
    num := (Z!Numerator(q)) mod m;
    return true, (num * InverseMod(den, m)) mod m;
end function;

function ResidueStatus(R, W, ell)
    okR, rr := ResidueInt(R, ell);
    okW, ww := ResidueInt(W, ell);
    if not (okR and okW) then return "denom", -1, -1, []; end if;
    F := GF(ell);
    labels := BoundaryLabelsFinite(F, F!rr, F!ww);
    if #labels gt 0 then return "boundary", rr, ww, labels; end if;
    return "open", rr, ww, [];
end function;

function FamilyPolynomial(R, w)
    t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
    A := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
    B := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);
    return x*A*B, t, A, B;
end function;

function RatToMod(q, N)
    den := (Z!Denominator(q)) mod N;
    if GCD(den, N) ne 1 then return false, Integers(N)!0; end if;
    num := (Z!Numerator(q)) mod N;
    return true, Integers(N)!((num * InverseMod(den, N)) mod N);
end function;

function SquareResidues(N)
    R := Integers(N);
    return {R!(i^2) : i in [0..N-1]};
end function;

function TangentCountsModN(RatR, RatW, N)
    Rmod := Integers(N);
    f, t, A, B := FamilyPolynomial(RatR, RatW);
    h := ExactQuotient(f, x);
    c1 := Coefficient(h, 1);
    c2 := Coefficient(h, 2);
    c3 := Coefficient(h, 3);
    c4 := Coefficient(h, 4);
    ok1, c1m := RatToMod(c1, N);
    ok2, c2m := RatToMod(c2, N);
    ok3, c3m := RatToMod(c3, N);
    ok4, c4m := RatToMod(c4, N);
    if not (ok1 and ok2 and ok3 and ok4) then
        return "denom", 0, 0, [0,0], [0,0];
    end if;

    Vrats := [RatR^2*RatW, -RatR^2*RatW];
    squares := SquareResidues(N);
    root_counts := [0, 0];
    tangent_counts := [0, 0];
    for vi in [1..2] do
        okV, Vm := RatToMod(Vrats[vi], N);
        if not okV then return "denom", 0, 0, [0,0], [0,0]; end if;
        for uu in [0..N-1] do
            U := Rmod!uu;
            M2 := c3m - 2*c4m*U;
            N2 := c1m - 2*c4m*U*Vm;
            Fval := 4*M2*N2 - (c2m - c4m*(U^2 + 2*Vm))^2;
            if Fval eq 0 then
                root_counts[vi] +:= 1;
                if M2 in squares and N2 in squares then
                    tangent_counts[vi] +:= 1;
                end if;
            end if;
        end for;
    end for;
    return "ok", &+root_counts, &+tangent_counts, root_counts, tangent_counts;
end function;

procedure IncrementBy(~A, key, n)
    if IsDefined(A, key) then A[key] +:= n; else A[key] := n; end if;
end procedure;

rows := ReadTraceFiles(trace_file_list);
printf "M18_416_FIRSTPOSS_L2_DIAG primes=%o trace_files=%o rows=%o\n",
    prime_list, trace_file_list, #rows;

summary := AssociativeArray();

for idx in [1..#rows] do
    row := rows[idx];
    f, t, A, B := FamilyPolynomial(row`R, row`W);
    printf "TRACE_PAIR idx=%o R=%o W=%o t=%o st=<%o,%o> gate=<%o,%o> plus=%o minus=%o rational_tangent=%o\n",
        idx, row`R, row`W, t, row`st0, row`st1, row`gateR, row`gateW,
        row`plus, row`minus, row`tangent;
    for ell in prime_list do
        rstatus, rr, ww, labels := ResidueStatus(row`R, row`W, ell);
        if rstatus ne "open" then
            printf "L2_PRIME idx=%o ell=%o residue=%o labels=%o verdict=SKIP_NONOPEN\n",
                idx, ell, rstatus, labels;
            IncrementBy(~summary, <ell, rstatus>, 1);
            continue;
        end if;
        s1, roots1, tang1, roots_by_v1, tang_by_v1 :=
            TangentCountsModN(row`R, row`W, ell);
        s2, roots2, tang2, roots_by_v2, tang_by_v2 :=
            TangentCountsModN(row`R, row`W, ell^2);
        if s1 ne "ok" or s2 ne "ok" then
            verdict := "DENOM_IN_TANGENT_MODEL";
        elif tang1 eq 0 then
            verdict := "NO_TANGENT_MOD_L";
        elif tang2 eq 0 then
            verdict := "DIES_MOD_L2";
        else
            verdict := "SURVIVES_MOD_L2";
        end if;
        printf "L2_PRIME idx=%o ell=%o residue=open r=%o w=%o modl_roots=%o modl_tangent=%o modl_roots_byV=%o modl_tangent_byV=%o modl2_roots=%o modl2_tangent=%o modl2_roots_byV=%o modl2_tangent_byV=%o verdict=%o\n",
            idx, ell, rr, ww, roots1, tang1, roots_by_v1, tang_by_v1,
            roots2, tang2, roots_by_v2, tang_by_v2, verdict;
        IncrementBy(~summary, <ell, verdict>, 1);
    end for;
end for;

print "L2_SUMMARY";
for key in Sort([k : k in Keys(summary)]) do
    printf "  ell=%o status=%o count=%o\n", key[1], key[2], summary[key];
end for;

quit;
