//////////////////////////////////////////////////////////////////////
//  Aux-prime sieve scanner for the M_1(8,4) [4,16] smooth p=7 strata.
//
//  This probes candidate good-reduction primes against the existing
//  mod-343 smooth-strata gate.  It does not run the expensive exact
//  halving chain; it reports how many height-bounded rational pairs
//  survive the current base aux primes, then how much each candidate
//  prime or candidate-prime pair cuts that survivor set.
//
//  Usage:
//    magma -b height:=30 base_primes:="11,13" \
//      candidate_primes:="17,19,23,29,31,37,41" \
//      gate_mod:=343 \
//      gate_file:=data/agent_m18_416_smooth_strata_mod343_gate.txt \
//      code/agent_m18_416_aux_prime_sieve_scan.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned height then height := 30;
elif Type(height) eq MonStgElt then height := StringToInteger(height); end if;
if not assigned base_primes then base_prime_list := [11,13];
elif Type(base_primes) eq MonStgElt then
    base_prime_list := (#base_primes eq 0) select []
        else [StringToInteger(s) : s in Split(base_primes, ",")];
else base_prime_list := base_primes; end if;
if not assigned candidate_primes then candidate_prime_list := [17,19,23,29,31,37,41];
elif Type(candidate_primes) eq MonStgElt then
    candidate_prime_list := (#candidate_primes eq 0) select []
        else [StringToInteger(s) : s in Split(candidate_primes, ",")];
else candidate_prime_list := candidate_primes; end if;
if not assigned gate_mod then gate_mod := 343;
elif Type(gate_mod) eq MonStgElt then gate_mod := StringToInteger(gate_mod); end if;
if not assigned gate_file then
    gate_file := "data/agent_m18_416_smooth_strata_mod343_gate.txt";
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

function ReadGateFile(filename, modulus)
    gate := {};
    st_by_residue := AssociativeArray();
    lines := Split(Read(filename), "\n");
    for raw in lines do
        s := raw;
        if #s gt 0 and s[#s] eq "\r" then s := s[1..#s-1]; end if;
        if #s eq 0 or s[1] eq "#" then continue; end if;
        parts := [part : part in Split(s, " ") | #part gt 0];
        if #parts lt 2 then continue; end if;
        pr := <StringToInteger(parts[1]) mod modulus,
               StringToInteger(parts[2]) mod modulus>;
        Include(~gate, pr);
        if #parts ge 4 then
            st_by_residue[pr] := <StringToInteger(parts[3]), StringToInteger(parts[4])>;
        else
            st_by_residue[pr] := <pr[1] mod 7, pr[2] mod 7>;
        end if;
    end for;
    return gate, st_by_residue;
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

function GoodHyperellipticPolynomial(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function FamilyPolynomialFinite(F, R, w)
    PF<X> := PolynomialRing(F);
    t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
    A := X^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*X + R^4;
    B := (R + 2 + 4*t)*X^2 + (R^2 + 4*R + 1 + 8*t)*X + (2*R^2 + R + 4*t);
    return X*A*B, t, A, B;
end function;

function YRFinite(F, R, w)
    Qfac := R^2 - (F!1/2)*R*w^2 + (F!1/2)*R - w^2;
    return -2*R*(R-1)^2*Qfac/(w^2-1);
end function;

function DivisibleByNFinite(J, D, n)
    G, phi := AbelianGroup(J);
    a := D @@ phi;
    coords := Eltseq(a);
    invs := Invariants(G);
    for i in [1..#coords] do
        g := GCD(n, invs[i]);
        if coords[i] mod g ne 0 then return false; end if;
    end for;
    return true;
end function;

function OpenTarget416Residues(pp)
    F := GF(pp);
    PF<X> := PolynomialRing(F);
    allowed := {};
    good_open := 0;
    first_count := 0;
    pr_count := 0;
    target_count := 0;

    for R in F do for w in F do
        if #BoundaryLabelsFinite(F, R, w) ne 0 then continue; end if;
        f, t, A, B := FamilyPolynomialFinite(F, R, w);
        if not GoodHyperellipticPolynomial(f) then continue; end if;
        yR := YRFinite(F, R, w);
        if yR^2 ne Evaluate(f, -R) then continue; end if;
        C := HyperellipticCurve(f);
        J := Jacobian(C);
        Tx := J![X, F!0];
        PR := J![X + R, yR];
        has_first := DivisibleByNFinite(J, Tx, 2);
        has_pr2 := DivisibleByNFinite(J, PR, 2);
        good_open +:= 1;
        if has_first then first_count +:= 1; end if;
        if has_pr2 then pr_count +:= 1; end if;
        if has_first and has_pr2 then
            target_count +:= 1;
            Include(~allowed, <Z!R, Z!w>);
        end if;
    end for; end for;
    return allowed, good_open, first_count, pr_count, target_count;
end function;

function RationalParametersOfHeight(B)
    vals := [];
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) eq 1 then Append(~vals, Q!num/den); end if;
        end for;
    end for;
    return Sort(Setseq(Seqset(vals)));
end function;

function ResidueInt(q, m)
    den := (Z!Denominator(q)) mod m;
    if GCD(den, m) ne 1 then return false, 0; end if;
    num := (Z!Numerator(q)) mod m;
    return true, (num * InverseMod(den, m)) mod m;
end function;

function IsSquareQ(q)
    q := Q!q;
    if q lt 0 then return false; end if;
    return IsSquare(Numerator(q)) and IsSquare(Denominator(q));
end function;

function PlusDisc(R, w)
    return -4*(w-R)*(R-1)^2*(R+1)*(R+w)*(w+1)*(R*w - 3*R + 3*w - 1);
end function;

function MinusDisc(R, w)
    return 4*(w-1)*(R+1)*(R*w + 3*R + 3*w + 1)
           *(R^4 - 2*R^3 + R^2*w^2 - R^2 + 2*R*w^2 - w^2);
end function;

function FirstCoverPossible(R, w)
    return IsSquareQ(PlusDisc(R,w)) or IsSquareQ(MinusDisc(R,w));
end function;

function AuxStatus(R, w, pp, F, target)
    okR, rp := ResidueInt(R, pp);
    okW, wp := ResidueInt(w, pp);
    if not (okR and okW) then return true, "den"; end if;
    if #BoundaryLabelsFinite(F, F!rp, F!wp) ne 0 then return true, "boundary"; end if;
    if <rp,wp> notin target then return false, "kill"; end if;
    return true, "target";
end function;

procedure IncrementBy(~A, key, n)
    if IsDefined(A, key) then A[key] +:= n; else A[key] := n; end if;
end procedure;

function AssocCount(A, key)
    return IsDefined(A, key) select A[key] else 0;
end function;

all_primes := Sort(Setseq(Seqset(base_prime_list cat candidate_prime_list)));
target_by_prime := AssociativeArray();
field_by_prime := AssociativeArray();

print "M18_416_AUX_PRIME_SIEVE_SCAN";
printf "height=%o base_primes=%o candidate_primes=%o gate_mod=%o gate_file=%o\n",
    height, base_prime_list, candidate_prime_list, gate_mod, gate_file;

for pp in all_primes do
    field_by_prime[pp] := GF(pp);
    allowed, good_open, first_count, pr_count, target_count := OpenTarget416Residues(pp);
    target_by_prime[pp] := allowed;
    printf "AUX_PRIME_PROFILE p=%o good_open=%o first_Tx_half=%o PR_half=%o target416=%o\n",
        pp, good_open, first_count, pr_count, target_count;
end for;

gate, st_by_residue := ReadGateFile(gate_file, gate_mod);
allowed_seq := Sort(Setseq(gate));
params := RationalParametersOfHeight(height);
buckets := AssociativeArray();
for r in params do
    ok, rv := ResidueInt(r, gate_mod);
    if not ok then continue; end if;
    if IsDefined(buckets, rv) then Append(~buckets[rv], r);
    else buckets[rv] := [r]; end if;
end for;

gate_budget := 0;
gate_residue_by_st := AssociativeArray();
gate_budget_by_st := AssociativeArray();
for pr in allowed_seq do
    st := IsDefined(st_by_residue, pr) select st_by_residue[pr] else <pr[1] mod 7, pr[2] mod 7>;
    IncrementBy(~gate_residue_by_st, st, 1);
    rc := IsDefined(buckets, pr[1]) select #buckets[pr[1]] else 0;
    wc := IsDefined(buckets, pr[2]) select #buckets[pr[2]] else 0;
    b := rc*wc;
    gate_budget +:= b;
    IncrementBy(~gate_budget_by_st, st, b);
end for;
printf "GATE height=%o parameters=%o residues=%o budget=%o\n",
    height, #params, #allowed_seq, gate_budget;

base_raw := 0;
base_cand := 0;
base_auxkill := 0;
base_firstposs := 0;
base_cand_by_st := AssociativeArray();
base_first_by_st := AssociativeArray();

single_cand := [0 : i in [1..#candidate_prime_list]];
single_kill := [0 : i in [1..#candidate_prime_list]];
single_unknown := [0 : i in [1..#candidate_prime_list]];
single_active_pass := [0 : i in [1..#candidate_prime_list]];
single_firstposs := [0 : i in [1..#candidate_prime_list]];

combos := [];
for i in [1..#candidate_prime_list] do
    for j in [i+1..#candidate_prime_list] do
        Append(~combos, <i,j>);
    end for;
end for;
combo_cand := [0 : i in [1..#combos]];
combo_firstposs := [0 : i in [1..#combos]];

for pr in allowed_seq do
    if not (IsDefined(buckets, pr[1]) and IsDefined(buckets, pr[2])) then continue; end if;
    st := IsDefined(st_by_residue, pr) select st_by_residue[pr] else <pr[1] mod 7, pr[2] mod 7>;
    for R in buckets[pr[1]] do
        for w in buckets[pr[2]] do
            if R eq 0 or w in {Q!-1, Q!0, Q!1} then continue; end if;
            base_raw +:= 1;
            okbase := true;
            for pp in base_prime_list do
                ok, status := AuxStatus(R, w, pp, field_by_prime[pp], target_by_prime[pp]);
                if not ok then okbase := false; break; end if;
            end for;
            if not okbase then
                base_auxkill +:= 1;
                continue;
            end if;

            base_cand +:= 1;
            IncrementBy(~base_cand_by_st, st, 1);
            first_ok := FirstCoverPossible(R, w);
            if first_ok then
                base_firstposs +:= 1;
                IncrementBy(~base_first_by_st, st, 1);
            end if;

            pass_flags := [];
            for i in [1..#candidate_prime_list] do
                pp := candidate_prime_list[i];
                ok, status := AuxStatus(R, w, pp, field_by_prime[pp], target_by_prime[pp]);
                Append(~pass_flags, ok);
                if ok then
                    single_cand[i] +:= 1;
                    if first_ok then single_firstposs[i] +:= 1; end if;
                    if status eq "target" then
                        single_active_pass[i] +:= 1;
                    else
                        single_unknown[i] +:= 1;
                    end if;
                else
                    single_kill[i] +:= 1;
                end if;
            end for;

            for ci in [1..#combos] do
                ij := combos[ci];
                if pass_flags[ij[1]] and pass_flags[ij[2]] then
                    combo_cand[ci] +:= 1;
                    if first_ok then combo_firstposs[ci] +:= 1; end if;
                end if;
            end for;
        end for;
    end for;
end for;

printf "SCAN_BASE raw=%o cand=%o auxkill=%o firstposs=%o\n",
    base_raw, base_cand, base_auxkill, base_firstposs;

for st in Sort(Setseq(Keys(gate_residue_by_st) join Keys(base_cand_by_st))) do
    printf "SCAN_STRATUM <%o,%o> gate_res=%o gate_budget=%o base_cand=%o base_firstposs=%o\n",
        st[1], st[2],
        AssocCount(gate_residue_by_st, st),
        AssocCount(gate_budget_by_st, st),
        AssocCount(base_cand_by_st, st),
        AssocCount(base_first_by_st, st);
end for;

for i in [1..#candidate_prime_list] do
    pp := candidate_prime_list[i];
    printf "SCAN_SINGLE p=%o cand=%o killed=%o unknown_pass=%o target_pass=%o firstposs=%o\n",
        pp, single_cand[i], single_kill[i], single_unknown[i],
        single_active_pass[i], single_firstposs[i];
end for;

for ci in [1..#combos] do
    ij := combos[ci];
    pp := candidate_prime_list[ij[1]];
    qq := candidate_prime_list[ij[2]];
    printf "SCAN_PAIR p=%o q=%o cand=%o firstposs=%o\n",
        pp, qq, combo_cand[ci], combo_firstposs[ci];
end for;

quit;
