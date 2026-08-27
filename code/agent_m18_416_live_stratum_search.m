
//////////////////////////////////////////////////////////////////////
//  [4,16] refined search on the LIVE p=7 blowup strata of M_1(8,4).
//
//  Follows code/agent_m18_416_p7_blowup.m and
//  code/agent_m18_416_p7_blowup_level2.m: the mod-7 boundary closure
//  (21 residues) refines to 449 surviving (R,w) residues mod 49, with
//  the fully-live stratum (R,w) == (0,0) mod 7 (i.e. 7|R, 7|w) as the
//  prime target.  This script replaces the mod-7 CRT gate of
//  agent_m18_416_search_crt.m by the mod-49 survivor gate and pushes
//  the rational height well beyond the old 20/30.
//
//  Exact chain (verbatim from agent_m18_416_search_crt.m):
//    FirstCoverPossible -> TangentCandidates -> exact Tx=[x,0] halving
//    -> exact P_R=(-R,Y_R) halving -> TorsionSubgroup.
//
//  Usage:
//    magma -b height:=60 aux_primes:="11,13" NParts:=6 Part:=0 \
//        agent_m18_416_live_stratum_search.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned height then height := 60;
elif Type(height) eq MonStgElt then height := StringToInteger(height); end if;
if not assigned aux_primes then aux_prime_list := [11,13];
elif Type(aux_primes) eq MonStgElt then
    aux_prime_list := (#aux_primes eq 0) select []
        else [StringToInteger(s) : s in Split(aux_primes, ",")];
else aux_prime_list := aux_primes; end if;
if not assigned tangent_sieve_primes then tangent_sieve_prime_list := [];
elif Type(tangent_sieve_primes) eq MonStgElt then
    tangent_sieve_prime_list := (#tangent_sieve_primes eq 0) select []
        else [StringToInteger(s) : s in Split(tangent_sieve_primes, ",")];
else tangent_sieve_prime_list := tangent_sieve_primes; end if;
if not assigned max_hits then max_hits := 10;
elif Type(max_hits) eq MonStgElt then max_hits := StringToInteger(max_hits); end if;
if not assigned progress_interval then progress_interval := 2000;
elif Type(progress_interval) eq MonStgElt then progress_interval := StringToInteger(progress_interval); end if;
if not assigned NParts then NParts := 1;
elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;
if not assigned Part then Part := 0;
elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned smooth_strata_only then
    smooth_strata_only := false;
elif Type(smooth_strata_only) eq MonStgElt then
    smooth_strata_only := smooth_strata_only in {"true", "True", "1", "yes"};
end if;
if not assigned trace_firstposs then
    trace_firstposs := false;
elif Type(trace_firstposs) eq MonStgElt then
    trace_firstposs := trace_firstposs in {"true", "True", "1", "yes"};
end if;
if not assigned trace_file then
    trace_file := "";
end if;
if not assigned gate_mod then
    if assigned gate_file then
        gate_mod := 343;
    else
        gate_mod := 49;
    end if;
elif Type(gate_mod) eq MonStgElt then
    gate_mod := StringToInteger(gate_mod);
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);
p := 7;
Fp := GF(p);

//------------------------------------------------------------------
// mod-49 survivor gate, recomputed inline (level-1 cokernel blowup)
//------------------------------------------------------------------
Rng7<R7,w7,a7,b7,c7,d7,e7> := PolynomialRing(Q, 7, "grevlex");
KF7 := FieldOfFractions(Rng7); PX7<x7> := PolynomialRing(KF7);
tt := (2*R7^2 + (1-w7^2)*R7 - 2*w7^2)/(4*(w7^2-1));
AA := x7^2 + (R7^3 + 4*R7^2*tt + R7 - 8*R7*tt + 4*tt)*x7 + R7^4;
BB := (R7 + 2 + 4*tt)*x7^2 + (R7^2 + 4*R7 + 1 + 8*tt)*x7 + (2*R7^2 + R7 + 4*tt);
ff := x7*AA*BB; cc4 := R7 + 2 + 4*tt;
qq := x7^2 + a7*x7 + b7; ll := c7*x7^2 + d7*x7 + e7;
FF := ff - ll^2 - cc4*(x7+R7)*qq^2;
E416 := [];
for i in [0..4] do
    ci := Numerator(Coefficient(FF, i));
    den := LCM([Denominator(co) : co in Coefficients(ci)]); ci := den*ci;
    g := GCD([Z!co : co in Coefficients(ci)]);
    Append(~E416, Rng7!(ci/g));
end for;
DerAll := [[Derivative(E416[i], j) : j in [1..7]] : i in [1..5]];
auxidx := [3,4,5,6,7];

function AuxMod7(R0, w0)
    sols := [];
    for aa in [0..6] do for bb in [0..6] do for cv in [0..6] do
      for dv in [0..6] do for ev in [0..6] do
        v := [R0,w0,aa,bb,cv,dv,ev]; ok := true;
        for i in [1..5] do
            if (Z!Evaluate(E416[i], v)) mod p ne 0 then ok := false; break; end if;
        end for;
        if ok then Append(~sols, [aa,bb,cv,dv,ev]); end if;
    end for; end for; end for; end for; end for;
    return sols;
end function;

boundaryPolys := [R7, w7, w7-1, w7+1, R7-1, R7+1, R7-w7, R7+w7,
  R7*w7-3*R7+3*w7-1, R7*w7+3*R7+3*w7+1,
  2*R7^2-R7*w7^2+R7-2*w7^2, R7^4-2*R7^3+R7^2*w7^2-R7^2+2*R7*w7^2-w7^2];

smoothStrata7 := {<3,3>, <3,4>, <4,0>, <5,0>, <5,2>, <5,5>};
survivors49 := {};   // set of <R mod 49, w mod 49>
if assigned gate_file then
    print "external gate_file supplied; skipping inline mod-49 survivor gate";
else
    print "computing mod-49 survivor gate (level-1 blowup) ...";
    for R0 in [0..6] do for w0 in [0..6] do
        isbd := false;
        for bp in boundaryPolys do
            if (Z!Evaluate(bp, [R0,w0,0,0,0,0,0])) mod p eq 0 then isbd := true; break; end if;
        end for;
        if not isbd then continue; end if;
        aux7 := AuxMod7(R0, w0);
        if #aux7 eq 0 then continue; end if;
        live := {};
        for s in aux7 do
            base := [R0, w0] cat s;
            Jaux := Matrix(Fp, 5, 5,
                [Fp!(Z!Evaluate(DerAll[i][auxidx[j]], base)) : i in [1..5], j in [1..5]]);
            Pk := KernelMatrix(Jaux);
            nr := Nrows(Pk);
            if nr eq 0 then
                live := live join {<rr,oo> : rr in [0..6], oo in [0..6]};
                continue;
            end if;
            c0 := [Fp!((Z!Evaluate(E416[i], base)) div p) : i in [1..5]];
            gR := [Fp!(Z!Evaluate(DerAll[i][1], base)) : i in [1..5]];
            gw := [Fp!(Z!Evaluate(DerAll[i][2], base)) : i in [1..5]];
            Pc  := [&+[Pk[r][j]*c0[j] : j in [1..5]] : r in [1..nr]];
            PgR := [&+[Pk[r][j]*gR[j] : j in [1..5]] : r in [1..nr]];
            Pgw := [&+[Pk[r][j]*gw[j] : j in [1..5]] : r in [1..nr]];
            for rho in [0..6] do for om in [0..6] do
                ok := true;
                for r in [1..nr] do
                    if Pc[r] + PgR[r]*(Fp!rho) + Pgw[r]*(Fp!om) ne 0 then ok := false; break; end if;
                end for;
                if ok then Include(~live, <rho,om>); end if;
            end for; end for;
        end for;
        for pr in live do
            Include(~survivors49, <R0 + 7*pr[1], w0 + 7*pr[2]>);
        end for;
    end for; end for;
    printf "mod-49 survivor gate: %o residues\n", #survivors49;

    if smooth_strata_only then
        smooth49 := {pr : pr in survivors49 | <pr[1] mod 7, pr[2] mod 7> in smoothStrata7};
        printf "smooth-strata-only gate: %o residues over %o mod-7 strata\n",
            #smooth49, #smoothStrata7;
        for st in Sort(Setseq(smoothStrata7)) do
            cnt := #{pr : pr in smooth49 | <pr[1] mod 7, pr[2] mod 7> eq st};
            printf "  smooth stratum <%o,%o>: mod49 residues %o\n", st[1], st[2], cnt;
        end for;
        survivors49 := smooth49;
    end if;
end if;

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
            st_by_residue[pr] := <pr[1] mod p, pr[2] mod p>;
        end if;
    end for;
    return gate, st_by_residue;
end function;

procedure IncrementBy(~A, key, n)
    if IsDefined(A, key) then A[key] +:= n; else A[key] := n; end if;
end procedure;

function AssocCount(A, key)
    return IsDefined(A, key) select A[key] else 0;
end function;

function GateStratum(pr, st_by_residue)
    if IsDefined(st_by_residue, pr) then return st_by_residue[pr]; end if;
    return <pr[1] mod p, pr[2] mod p>;
end function;

//------------------------------------------------------------------
// exact-search machinery (from agent_m18_416_search_crt.m)
//------------------------------------------------------------------
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
        if DivisibleByNFinite(J, Tx, 2) and DivisibleByNFinite(J, PR, 2) then
            Include(~allowed, <Z!R,Z!w>);
        end if;
    end for; end for;
    return allowed;
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

function RatToField(q, F, ell)
    den := (Z!Denominator(q)) mod ell;
    if GCD(den, ell) ne 1 then return false, F!0; end if;
    num := (Z!Numerator(q)) mod ell;
    return true, F!((num * InverseMod(den, ell)) mod ell);
end function;

function IsSquareQ(q)
    q := Q!q;
    if q lt 0 then return false; end if;
    return IsSquare(Numerator(q)) and IsSquare(Denominator(q));
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do L := LCM(L, Denominator(Coefficient(f, i))); end for;
    return P!(L^2*f), L;
end function;

function FamilyPolynomial(R, w)
    t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
    A := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
    B := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);
    return x*A*B, t, A, B;
end function;

function YRational(R, w)
    Qfac := R^2 - (Q!1/2)*R*w^2 + (Q!1/2)*R - w^2;
    return -2*R*(R-1)^2*Qfac/(w^2-1);
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

function FirstCoverSigns(R, w)
    return IsSquareQ(PlusDisc(R,w)), IsSquareQ(MinusDisc(R,w));
end function;

function TangentCandidates(f, R, w)
    h := ExactQuotient(f, x);
    c1 := Coefficient(h, 1); c2 := Coefficient(h, 2);
    c3 := Coefficient(h, 3); c4 := Coefficient(h, 4);
    out := [];
    for V in [R^2*w, -R^2*w] do
        F := 4*(c3 - 2*c4*x)*(c1 - 2*c4*x*V) - (c2 - c4*(x^2 + 2*V))^2;
        for rt in Roots(F) do
            U0 := rt[1];
            M2 := c3 - 2*c4*U0;
            N2 := c1 - 2*c4*U0*V;
            if IsSquareQ(M2) and IsSquareQ(N2) then
                Append(~out, <U0, V, M2, N2>);
            end if;
        end for;
    end for;
    return out;
end function;

function TangentCongruenceStatus(R, w, ell)
    okR, rr := ResidueInt(R, ell);
    okW, ww := ResidueInt(w, ell);
    if not (okR and okW) then return "skip_denom", 0, 0; end if;

    F := GF(ell);
    if #BoundaryLabelsFinite(F, F!rr, F!ww) ne 0 then
        return "skip_boundary", 0, 0;
    end if;

    f, t, A, B := FamilyPolynomial(R, w);
    h := ExactQuotient(f, x);
    c1 := Coefficient(h, 1); c2 := Coefficient(h, 2);
    c3 := Coefficient(h, 3); c4 := Coefficient(h, 4);
    ok1, c1m := RatToField(c1, F, ell);
    ok2, c2m := RatToField(c2, F, ell);
    ok3, c3m := RatToField(c3, F, ell);
    ok4, c4m := RatToField(c4, F, ell);
    if not (ok1 and ok2 and ok3 and ok4) then
        return "skip_denom", 0, 0;
    end if;

    roots := 0;
    tangents := 0;
    for Vrat in [R^2*w, -R^2*w] do
        okV, Vm := RatToField(Vrat, F, ell);
        if not okV then return "skip_denom", 0, 0; end if;
        for U in F do
            M2 := c3m - 2*c4m*U;
            N2 := c1m - 2*c4m*U*Vm;
            Fval := 4*M2*N2 - (c2m - c4m*(U^2 + 2*Vm))^2;
            if Fval eq 0 then
                roots +:= 1;
                if IsSquare(M2) and IsSquare(N2) then
                    tangents +:= 1;
                end if;
            end if;
        end for;
    end for;
    if tangents eq 0 then return "kill", roots, tangents; end if;
    return "pass", roots, tangents;
end function;

function TangentSieveKillPrime(R, w, prime_list)
    for ell in prime_list do
        status, roots, tangents := TangentCongruenceStatus(R, w, ell);
        if status eq "kill" then return true, ell, roots; end if;
    end for;
    return false, 0, 0;
end function;

function Has416(invs)
    vals := Reverse(Sort([Valuation(n, 2) : n in invs]));
    return #vals ge 2 and vals[1] ge 4 and vals[2] ge 2;
end function;

//------------------------------------------------------------------
// main loop
//------------------------------------------------------------------
print "M_1(8,4) [4,16] live-stratum gated search";
printf "height=%o aux_primes=%o tangent_sieve_primes=%o Part=%o/%o smooth_strata_only=%o gate_mod=%o\n",
    height, aux_prime_list, tangent_sieve_prime_list, Part, NParts,
    smooth_strata_only, gate_mod;
if assigned gate_file then
    printf "gate_file=%o\n", gate_file;
end if;
trace_has_file := #trace_file gt 0;
trace_out := 0;
if trace_has_file then
    trace_out := Open(trace_file, "w");
    fprintf trace_out, "# columns: R W st0 st1 gate_R gate_W plus minus tangent\n";
    printf "trace_file=%o\n", trace_file;
end if;

open_target_by_prime := AssociativeArray();
for pp in aux_prime_list do
    open_target_by_prime[pp] := OpenTarget416Residues(pp);
    printf "aux_prime %o: open target416 residues %o\n", pp, #open_target_by_prime[pp];
end for;

params := RationalParametersOfHeight(height);
// bucket parameters by the active residue modulus (skip nonunit denominators)
buckets := AssociativeArray();
for r in params do
    ok, rv := ResidueInt(r, gate_mod);
    if not ok then continue; end if;
    if IsDefined(buckets, rv) then Append(~buckets[rv], r);
    else buckets[rv] := [r]; end if;
end for;

st_by_residue := AssociativeArray();
if assigned gate_file then
    allowed_set, st_by_residue := ReadGateFile(gate_file, gate_mod);
else
    if gate_mod ne 49 then
        error "gate_mod other than 49 requires gate_file";
    end if;
    allowed_set := survivors49;
    for pr in allowed_set do
        st_by_residue[pr] := <pr[1] mod p, pr[2] mod p>;
    end for;
end if;

allowed_seq := Sort(Setseq(allowed_set));
budget := 0;
gate_residue_by_st := AssociativeArray();
gate_budget_by_st := AssociativeArray();
for pr in allowed_seq do
    st := GateStratum(pr, st_by_residue);
    IncrementBy(~gate_residue_by_st, st, 1);
    rc := IsDefined(buckets, pr[1]) select #buckets[pr[1]] else 0;
    wc := IsDefined(buckets, pr[2]) select #buckets[pr[2]] else 0;
    b := rc*wc;
    budget +:= b;
    IncrementBy(~gate_budget_by_st, st, b);
end for;
printf "parameters=%o  gate_residues=%o  mod%o-gated pair budget=%o\n",
    #params, #allowed_seq, gate_mod, budget;

candidates := 0; family_smooth := 0; first_possible := 0;
tangent_bases := 0; exact_tests := 0; first_verified := 0;
pr_halved := 0; hits := 0; auxkill := 0; tangent_sieve_kill := 0;
torsion_counts := AssociativeArray();
work_residue_by_st := AssociativeArray();
work_budget_by_st := AssociativeArray();
auxkill_by_st := AssociativeArray();
cand_by_st := AssociativeArray();
smooth_by_st := AssociativeArray();
first_by_st := AssociativeArray();
tangent_sieve_kill_by_st := AssociativeArray();
tangent_by_st := AssociativeArray();
exact_by_st := AssociativeArray();
firstver_by_st := AssociativeArray();
prfail_by_st := AssociativeArray();
prhalf_by_st := AssociativeArray();
hit_by_st := AssociativeArray();
stop := false;
pairidx := 0;

for pr in allowed_seq do
    if stop then break; end if;
    pairidx +:= 1;
    if (pairidx mod NParts) ne Part then continue; end if;
    rr := pr[1]; ww := pr[2];
    if not (IsDefined(buckets, rr) and IsDefined(buckets, ww)) then continue; end if;
    st := GateStratum(pr, st_by_residue);
    IncrementBy(~work_residue_by_st, st, 1);
    IncrementBy(~work_budget_by_st, st, #buckets[rr] * #buckets[ww]);
    for R in buckets[rr] do
        if stop then break; end if;
        for w in buckets[ww] do
            if stop then break; end if;
            if R eq 0 or w in {Q!-1, Q!0, Q!1} then continue; end if;

            // aux good-open filters
            okaux := true;
            for pp in aux_prime_list do
                okR, rp := ResidueInt(R, pp);
                okW, wp := ResidueInt(w, pp);
                if not (okR and okW) then continue; end if;
                F := GF(pp);
                if #BoundaryLabelsFinite(F, F!rp, F!wp) ne 0 then continue; end if;
                if <rp,wp> notin open_target_by_prime[pp] then okaux := false; break; end if;
            end for;
            if not okaux then
                auxkill +:= 1;
                IncrementBy(~auxkill_by_st, st, 1);
                continue;
            end if;

            candidates +:= 1;
            IncrementBy(~cand_by_st, st, 1);
            if progress_interval gt 0 and candidates mod progress_interval eq 0 then
                printf "progress cand=%o smooth=%o firstposs=%o tansievekill=%o tangent=%o exact=%o firstver=%o prhalf=%o hits=%o\n",
                    candidates, family_smooth, first_possible, tangent_sieve_kill,
                    tangent_bases, exact_tests, first_verified,
                    pr_halved, hits;
            end if;

            f, t, A, B := FamilyPolynomial(R, w);
            if not GoodHyperellipticPolynomial(f) then continue; end if;
            family_smooth +:= 1;
            IncrementBy(~smooth_by_st, st, 1);
            plus_ok, minus_ok := FirstCoverSigns(R, w);
            if not (plus_ok or minus_ok) then continue; end if;
            first_possible +:= 1;
            IncrementBy(~first_by_st, st, 1);
            killed_by_tangent_sieve, kill_ell, kill_roots :=
                TangentSieveKillPrime(R, w, tangent_sieve_prime_list);
            if killed_by_tangent_sieve and not trace_firstposs then
                tangent_sieve_kill +:= 1;
                IncrementBy(~tangent_sieve_kill_by_st, st, 1);
                continue;
            end if;

            cands := TangentCandidates(f, R, w);
            if trace_firstposs then
                printf "FIRSTPOSS R=%o w=%o t=%o stratum=<%o,%o> mod%o=<%o,%o> plus=%o minus=%o tangent=%o tansievekill=%o kill_ell=%o\n",
                    R, w, t, st[1], st[2], gate_mod, rr, ww, plus_ok,
                    minus_ok, #cands, killed_by_tangent_sieve, kill_ell;
                if trace_has_file then
                    fprintf trace_out, "%o %o %o %o %o %o %o %o %o\n",
                        R, w, st[1], st[2], rr, ww,
                        plus_ok select 1 else 0,
                        minus_ok select 1 else 0,
                        #cands;
                end if;
            end if;
            if killed_by_tangent_sieve then
                tangent_sieve_kill +:= 1;
                IncrementBy(~tangent_sieve_kill_by_st, st, 1);
                continue;
            end if;
            if #cands eq 0 then continue; end if;
            tangent_bases +:= 1;
            IncrementBy(~tangent_by_st, st, 1);

            exact_tests +:= 1;
            IncrementBy(~exact_by_st, st, 1);
            fI, L := IntegralModelPolynomial(f);
            C := HyperellipticCurve(fI);
            J := Jacobian(C);
            Tx := J![x, Q!0];
            tx_div, tx_half := IsDivisibleBy(Tx, 2);
            if not tx_div then continue; end if;
            first_verified +:= 1;
            IncrementBy(~firstver_by_st, st, 1);

            yR := YRational(R, w);
            yRI := L*yR;
            if Evaluate(fI, -R) ne yRI^2 then continue; end if;
            PR := J![x + R, yRI];
            pr_div, pr_half := IsDivisibleBy(PR, 2);
            if not pr_div then
                IncrementBy(~prfail_by_st, st, 1);
                continue;
            end if;
            pr_halved +:= 1;
            IncrementBy(~prhalf_by_st, st, 1);
            G, phi := TorsionSubgroup(J);
            invs := Invariants(G);
            key := Sprint(invs);
            if IsDefined(torsion_counts, key) then torsion_counts[key] +:= 1;
            else torsion_counts[key] := 1; end if;
            printf "PR_HALF R=%o w=%o t=%o mod%o=<%o,%o> torsion=%o\n",
                R, w, t, gate_mod, rr, ww, invs;
            printf "  f = %o\n", fI;
            if Has416(invs) then
                hits +:= 1;
                IncrementBy(~hit_by_st, st, 1);
                printf "HIT_416 R=%o w=%o torsion=%o\n", R, w, invs;
                if hits ge max_hits and max_hits gt 0 then stop := true; end if;
            end if;
        end for;
    end for;
end for;

printf "DONE cand=%o auxkill=%o smooth=%o firstposs=%o tansievekill=%o tangent=%o exact=%o firstver=%o prhalf=%o hits=%o\n",
    candidates, auxkill, family_smooth, first_possible, tangent_sieve_kill,
    tangent_bases, exact_tests, first_verified, pr_halved, hits;
print "TORSION_COUNTS";
for key in Sort([k : k in Keys(torsion_counts)]) do
    printf "  %o : %o\n", key, torsion_counts[key];
end for;
strata_keys := {};
for st in Keys(gate_residue_by_st) do Include(~strata_keys, st); end for;
for st in Keys(work_residue_by_st) do Include(~strata_keys, st); end for;
for st in Keys(auxkill_by_st) do Include(~strata_keys, st); end for;
for st in Keys(cand_by_st) do Include(~strata_keys, st); end for;
for st in Keys(first_by_st) do Include(~strata_keys, st); end for;
for st in Keys(tangent_sieve_kill_by_st) do Include(~strata_keys, st); end for;
for st in Keys(prfail_by_st) do Include(~strata_keys, st); end for;
for st in Keys(prhalf_by_st) do Include(~strata_keys, st); end for;
print "P7_STRATUM_COUNTS";
for st in Sort(Setseq(strata_keys)) do
    printf "  <%o,%o> gate_res=%o gate_budget=%o work_res=%o work_budget=%o auxkill=%o cand=%o smooth=%o firstposs=%o tansievekill=%o tangent=%o exact=%o firstver=%o prfail=%o prhalf=%o hits=%o\n",
        st[1], st[2],
        AssocCount(gate_residue_by_st, st),
        AssocCount(gate_budget_by_st, st),
        AssocCount(work_residue_by_st, st),
        AssocCount(work_budget_by_st, st),
        AssocCount(auxkill_by_st, st),
        AssocCount(cand_by_st, st),
        AssocCount(smooth_by_st, st),
        AssocCount(first_by_st, st),
        AssocCount(tangent_sieve_kill_by_st, st),
        AssocCount(tangent_by_st, st),
        AssocCount(exact_by_st, st),
        AssocCount(firstver_by_st, st),
        AssocCount(prfail_by_st, st),
        AssocCount(prhalf_by_st, st),
        AssocCount(hit_by_st, st);
end for;
quit;
