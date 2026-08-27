//////////////////////////////////////////////////////////////////////
//  Denominator-aware p=7 boundary search for [4,16] in M_1(8,4).
//
//  This complements the strict affine CRT searches: rational parameters
//  with 7 dividing Denominator(R) or Denominator(w) reduce to the p=7
//  boundary of P^1_R x P^1_w and are therefore skipped by affine
//  residue vectors.
//
//  Charts covered here:
//      Rinf_w      R = 1/z, z = 0, w finite mod 7
//      r_Winf      w = 1/z, z = 0, R finite mod 7
//      Rinf_Winf   R = 1/zR, w = 1/zW, zR = zW = 0
//
//  For each chart, we first test the cleared FIRST_COVER and TARGET_416
//  closures over F_7.  The rational height search then keeps only pairs
//  landing on those live p=7 chart residues, applies optional good-prime
//  finite filters, and finally exact-checks Tx/2 and P_R/2 on J(Q).
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned height then
    height := 20;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

if not assigned aux_primes then
    aux_prime_list := [11,13,17];
elif Type(aux_primes) eq MonStgElt then
    if aux_primes eq "" or aux_primes eq "none" then
        aux_prime_list := [];
    else
        aux_prime_list := [StringToInteger(s) : s in Split(aux_primes, ",")];
    end if;
else
    aux_prime_list := aux_primes;
end if;

if not assigned max_exact then
    max_exact := 0; // 0 means no cap.
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;

if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;

if not assigned progress_interval then
    progress_interval := 250;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

procedure Increment(~A, key)
    if IsDefined(A, key) then
        A[key] +:= 1;
    else
        A[key] := 1;
    end if;
end procedure;

procedure AddSample(~A, key, val, max_samples)
    if not IsDefined(A, key) then
        A[key] := [];
    end if;
    if #A[key] lt max_samples then
        Append(~A[key], val);
    end if;
end procedure;

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            r := Q!num/den;
            key := Sprint(r);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, r);
            end if;
        end for;
    end for;
    return vals;
end function;

function GoodHyperellipticPolynomial(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function FamilyPolynomial(R, w)
    t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
    A := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
    B := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x
         + (2*R^2 + R + 4*t);
    return x*A*B, t, A, B;
end function;

function YRExact(R, w)
    Qfac := R^2 - (Q!1/2)*R*w^2 + (Q!1/2)*R - w^2;
    return -2*R*(R-1)^2*Qfac/(w^2-1);
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function IsSquareQ(q)
    q := Q!q;
    if q lt 0 then
        return false;
    end if;
    return IsSquare(Numerator(q)) and IsSquare(Denominator(q));
end function;

function SquareRootsQ(q)
    q := Q!q;
    if q eq 0 then
        return [Q!0];
    end if;
    if not IsSquareQ(q) then
        return [];
    end if;
    okn, rn := IsSquare(Numerator(q));
    okd, rd := IsSquare(Denominator(q));
    assert okn and okd;
    r := Q!rn/Q!rd;
    return [r, -r];
end function;

function TangentCandidates(f, R, w)
    h := ExactQuotient(f, x);
    c1 := Coefficient(h, 1);
    c2 := Coefficient(h, 2);
    c3 := Coefficient(h, 3);
    c4 := Coefficient(h, 4);

    out := [];
    PRing<Uvar> := PolynomialRing(Q);
    for V0 in [R^2*w, -R^2*w] do
        F := 4*(c3 - 2*c4*Uvar)*(c1 - 2*c4*Uvar*V0)
             - (c2 - c4*(Uvar^2 + 2*V0))^2;
        for rt in Roots(F) do
            U0 := rt[1];
            M2 := c3 - 2*c4*U0;
            N2 := c1 - 2*c4*U0*V0;
            if IsSquareQ(M2) and IsSquareQ(N2) then
                Append(~out, <U0, V0, M2, N2>);
            end if;
        end for;
    end for;
    return out;
end function;

function FirstCoverSolutionsExactFull(f)
    h := ExactQuotient(f, x);
    c0 := Coefficient(h, 0);
    c1 := Coefficient(h, 1);
    c2 := Coefficient(h, 2);
    c3 := Coefficient(h, 3);
    c4 := Coefficient(h, 4);

    out := [];
    PRing<Uvar> := PolynomialRing(Q);
    for V0 in SquareRootsQ(c0/c4) do
        F := 4*(c3 - 2*c4*Uvar)*(c1 - 2*c4*Uvar*V0)
             - (c2 - c4*(Uvar^2 + 2*V0))^2;
        for rt in Roots(F) do
            U0 := rt[1];
            M2 := c3 - 2*c4*U0;
            N2 := c1 - 2*c4*U0*V0;
            for M0 in SquareRootsQ(M2) do
                for N0 in SquareRootsQ(N2) do
                    q0 := x^2 + U0*x + V0;
                    identity := h - x*(M0*x+N0)^2 - c4*q0^2;
                    if identity eq 0 then
                        Append(~out, <U0,V0,M0,N0>);
                    end if;
                end for;
            end for;
        end for;
    end for;
    return out;
end function;

function TorsionOrder(invs)
    if #invs eq 0 then
        return 1;
    end if;
    return &*invs;
end function;

function Exponent(invs)
    if #invs eq 0 then
        return 1;
    end if;
    return invs[#invs];
end function;

function Has416Torsion(invs)
    vals := Reverse(Sort([Valuation(n, 2) : n in invs]));
    return #vals ge 2 and vals[1] ge 4 and vals[2] ge 2;
end function;

function ResidueInt(q, p)
    den := Z!Denominator(q);
    if den mod p eq 0 then
        return false, 0;
    end if;
    num := Z!Numerator(q);
    return true, ((num mod p) * InverseMod(den mod p, p)) mod p;
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

function FamilyPolynomialFinite(F, R, w)
    PF<X> := PolynomialRing(F);
    t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
    A := X^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*X + R^4;
    B := (R + 2 + 4*t)*X^2 + (R^2 + 4*R + 1 + 8*t)*X
         + (2*R^2 + R + 4*t);
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
        if coords[i] mod g ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function GoodPrimeAllowedTable(p)
    F := GF(p);
    allowed := [ true : i in [1..p^2] ];
    boundary := 0;
    good := 0;
    target416 := 0;
    killed := 0;

    for R in F do
        for w in F do
            idx := (Z!R)*p + (Z!w) + 1;
            labels := BoundaryLabelsFinite(F, R, w);
            if #labels ne 0 then
                boundary +:= 1;
                allowed[idx] := true;
                continue;
            end if;

            f, t, A, B := FamilyPolynomialFinite(F, R, w);
            if not GoodHyperellipticPolynomial(f) then
                boundary +:= 1;
                allowed[idx] := true;
                continue;
            end if;

            PF<X> := PolynomialRing(F);
            yR := YRFinite(F, R, w);
            C := HyperellipticCurve(f);
            J := Jacobian(C);
            Tx := J![X, F!0];
            PR := J![X + R, yR];
            has_first := DivisibleByNFinite(J, Tx, 2);
            has_pr2 := DivisibleByNFinite(J, PR, 2);

            good +:= 1;
            if has_first and has_pr2 then
                target416 +:= 1;
                allowed[idx] := true;
            else
                killed +:= 1;
                allowed[idx] := false;
            end if;
        end for;
    end for;

    return allowed, boundary, good, target416, killed;
end function;

function AllZero(eqs, vals)
    for pol in eqs do
        if Evaluate(pol, vals) ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function BuildChartEquations(F, chart)
    S<r,w,U,V,M,N,a,b,c,d,e> := PolynomialRing(F, 11, "grevlex");
    K := FieldOfFractions(S);
    PX<X> := PolynomialRing(K);

    if chart eq "Rinf_w" then
        Rexpr := 1/(K!r);
        wexpr := K!w;
    elif chart eq "r_Winf" then
        Rexpr := K!r;
        wexpr := 1/(K!w);
    elif chart eq "Rinf_Winf" then
        Rexpr := 1/(K!r);
        wexpr := 1/(K!w);
    else
        Rexpr := K!r;
        wexpr := K!w;
    end if;

    t := (2*Rexpr^2 + (1-wexpr^2)*Rexpr - 2*wexpr^2)/(4*(wexpr^2-1));
    A := X^2 + (Rexpr^3 + 4*Rexpr^2*t + Rexpr - 8*Rexpr*t + 4*t)*X
         + Rexpr^4;
    B := (Rexpr + 2 + 4*t)*X^2 + (Rexpr^2 + 4*Rexpr + 1 + 8*t)*X
         + (2*Rexpr^2 + Rexpr + 4*t);
    f := X*A*B;
    h := A*B;
    c4 := Coefficient(h, 4);

    a0 := X^2 + U*X + V;
    F0 := h - X*(M*X+N)^2 - c4*a0^2;
    E0 := [Numerator(Coefficient(F0, i)) : i in [0..4]];

    q416 := X^2 + a*X + b;
    ell416 := c*X^2 + d*X + e;
    F416 := f - ell416^2 - c4*(X+Rexpr)*q416^2;
    E416 := [Numerator(Coefficient(F416, i)) : i in [0..4]];

    return S, E0, E416;
end function;

function HasFirstClosure(E0, elems, base1, base2)
    z := elems[1];
    for UU in elems do
        for VV in elems do
            for MM in elems do
                for NN in elems do
                    vals := [base1, base2, UU, VV, MM, NN, z, z, z, z, z];
                    if AllZero(E0, vals) then
                        return true, <UU,VV,MM,NN>;
                    end if;
                end for;
            end for;
        end for;
    end for;
    return false, _;
end function;

function Has416Closure(E416, elems, base1, base2)
    z := elems[1];
    for aa in elems do
        for bb in elems do
            for cc in elems do
                for dd in elems do
                    for ee in elems do
                        vals := [base1, base2, z, z, z, z, aa, bb, cc, dd, ee];
                        if AllZero(E416, vals) then
                            return true, <aa,bb,cc,dd,ee>;
                        end if;
                    end for;
                end for;
            end for;
        end for;
    end for;
    return false, _;
end function;

function ChartKey(chart, i)
    if chart eq "Rinf_w" then
        return "Rinf_w:" cat IntegerToString(i);
    elif chart eq "r_Winf" then
        return "r_Winf:" cat IntegerToString(i);
    elif chart eq "Rinf_Winf" then
        return "Rinf_Winf";
    end if;
    return "unknown";
end function;

function P7BoundaryKey(R, w)
    okR, r7 := ResidueInt(R, 7);
    okW, w7 := ResidueInt(w, 7);
    if okR and okW then
        return false, "finite", r7, w7;
    elif (not okR) and okW then
        return true, "Rinf_w", 0, w7;
    elif okR and (not okW) then
        return true, "r_Winf", r7, 0;
    else
        return true, "Rinf_Winf", 0, 0;
    end if;
end function;

function AuxAllows(R, w, profiles)
    for prof in profiles do
        p := prof[1];
        allowed := prof[2];
        okR, rp := ResidueInt(R, p);
        okW, wp := ResidueInt(w, p);
        if not (okR and okW) then
            continue;
        end if;
        idx := rp*p + wp + 1;
        if not allowed[idx] then
            return false, p, rp, wp;
        end if;
    end for;
    return true, 0, 0, 0;
end function;

print "M_1(8,4) [4,16] denominator-aware p=7 boundary search";
print "height", height, "aux_primes", aux_prime_list,
      "max_exact", max_exact, "max_hits", max_hits,
      "progress_interval", progress_interval;
print "agent prefix", "agent_m18_416_den7";
print "";

F7 := GF(7);
elems7 := [u : u in F7];
p7_live := {};
p7_first_keys := {};
p7_samples := AssociativeArray();

for chart in ["Rinf_w", "r_Winf", "Rinf_Winf"] do
    S, E0, E416 := BuildChartEquations(F7, chart);
    bases := chart eq "Rinf_Winf" select [0] else [0..6];
    first_count := 0;
    target_count := 0;
    print "P7_CHART_BEGIN", chart;
    for i in bases do
        b1 := F7!0;
        b2 := F7!0;
        if chart eq "Rinf_w" then
            b1 := F7!0;
            b2 := F7!i;
        elif chart eq "r_Winf" then
            b1 := F7!i;
            b2 := F7!0;
        end if;

        has_first, first_sol := HasFirstClosure(E0, elems7, b1, b2);
        has_416, sol416 := Has416Closure(E416, elems7, b1, b2);
        key := ChartKey(chart, i);
        if has_first then
            first_count +:= 1;
            Include(~p7_first_keys, key);
        end if;
        if has_first and has_416 then
            target_count +:= 1;
            Include(~p7_live, key);
            p7_samples[key] := <first_sol, sol416>;
        end if;
        print "  base", key, "first_closure", has_first,
              "target416_closure", has_first and has_416;
        if has_first and has_416 then
            print "    sample_first", first_sol, "sample_416", sol416;
        end if;
    end for;
    print "P7_CHART_SUMMARY", chart, "bases", #bases,
          "first_closure", first_count, "target416_closure", target_count;
end for;
print "P7_LIVE_TARGET_KEYS", Sort(Setseq(p7_live));
print "";

aux_profiles := [* *];
for p in aux_prime_list do
    allowed, boundary, good, target416, killed := GoodPrimeAllowedTable(p);
    Append(~aux_profiles, <p, allowed>);
    print "AUX_PROFILE", p, "boundary_or_bad", boundary,
          "good_open", good, "good_target416", target416,
          "good_killed", killed;
end for;
print "";

params := RationalParametersOfHeight(height);
print "PARAMETERS", #params;

total_pairs := 0;
den7_pairs := 0;
p7_chart_filtered := 0;
p7_chart_pass := 0;
aux_filtered := 0;
aux_pass := 0;
smooth_tests := 0;
smooth := 0;
yr_fail := 0;
exact_tests := 0;
tx_halves := 0;
pr_halves := 0;
both_halves := 0;
torsion_tests := 0;
pr_divisibility_errors := 0;
hits := [];

chart_seen := AssociativeArray();
chart_pass_counts := AssociativeArray();
chart_exact_counts := AssociativeArray();
chart_tx_counts := AssociativeArray();
chart_pr_counts := AssociativeArray();
chart_both_counts := AssociativeArray();
chart_samples := AssociativeArray();
aux_kill_counts := AssociativeArray();

for R in params do
    for w in params do
        total_pairs +:= 1;
        is_den7, chart, rkey, wkey := P7BoundaryKey(R, w);
        if not is_den7 then
            continue;
        end if;
        den7_pairs +:= 1;

        key := chart eq "Rinf_w" select ChartKey(chart, wkey)
               else chart eq "r_Winf" select ChartKey(chart, rkey)
               else ChartKey(chart, 0);
        Increment(~chart_seen, key);

        if key notin p7_live then
            p7_chart_filtered +:= 1;
            continue;
        end if;
        p7_chart_pass +:= 1;
        Increment(~chart_pass_counts, key);
        AddSample(~chart_samples, key, <R,w>, 6);

        ok_aux, kill_p, kill_r, kill_w := AuxAllows(R, w, aux_profiles);
        if not ok_aux then
            aux_filtered +:= 1;
            kill_key := IntegerToString(kill_p) cat ":" cat IntegerToString(kill_r)
                        cat "," cat IntegerToString(kill_w);
            Increment(~aux_kill_counts, kill_key);
            continue;
        end if;
        aux_pass +:= 1;

        if w eq 1 or w eq -1 then
            continue;
        end if;
        smooth_tests +:= 1;
        f, t, A, B := FamilyPolynomial(R, w);
        if not GoodHyperellipticPolynomial(f) then
            continue;
        end if;
        smooth +:= 1;

        yR := YRExact(R, w);
        if yR^2 ne Evaluate(f, -R) then
            yr_fail +:= 1;
            print "WARNING_YR_FAIL", "R", R, "w", w, "t", t;
            continue;
        end if;

        tx_candidates := FirstCoverSolutionsExactFull(f);
        ok_tx := #tx_candidates gt 0;

        fI, L := IntegralModelPolynomial(f);
        C := HyperellipticCurve(fI);
        J := Jacobian(C);
        PR := J![x + R, Q!(L*yR)];

        exact_tests +:= 1;
        Increment(~chart_exact_counts, key);
        if progress_interval gt 0 and exact_tests mod progress_interval eq 0 then
            print "PROGRESS", "exact_tests", exact_tests,
                  "tx_halves", tx_halves, "pr_halves", pr_halves,
                  "both_halves", both_halves, "hits", #hits;
        end if;

        ok_pr := false;
        try
            ok_pr, Hpr := IsDivisibleBy(PR, 2);
        catch err1
            try
                Cq := HyperellipticCurve(f);
                Jq := Jacobian(Cq);
                PRq := Jq![x + R, yR];
                ok_pr, Hprq := IsDivisibleBy(PRq, 2);
                print "WARNING_PR_DIV_FALLBACK_RATIONAL_MODEL",
                      "R", R, "w", w, "chart", key, "ok_pr", ok_pr;
            catch err2
                pr_divisibility_errors +:= 1;
                print "WARNING_PR_DIV_ERROR", "R", R, "w", w,
                      "chart", key, "t", t;
            end try;
        end try;
        if ok_tx then
            tx_halves +:= 1;
            Increment(~chart_tx_counts, key);
        end if;
        if ok_pr then
            pr_halves +:= 1;
            Increment(~chart_pr_counts, key);
            print "PR_HALF", "R", R, "w", w, "chart", key, "t", t,
                  "Tx_half", ok_tx;
        end if;
        if ok_tx and ok_pr then
            both_halves +:= 1;
            Increment(~chart_both_counts, key);
            torsion_tests +:= 1;
            Gtors, phi := TorsionSubgroup(J);
            invs := Invariants(Gtors);
            Append(~hits, <R,w,t,key,invs,fI>);
            print "HIT_416_HALVING", "R", R, "w", w, "chart", key,
                  "t", t, "torsion", invs, "order", TorsionOrder(invs),
                  "exponent", Exponent(invs), "has416_torsion", Has416Torsion(invs);
            print "  fI", fI;
            if #hits ge max_hits then
                break R;
            end if;
        end if;

        if max_exact gt 0 and exact_tests ge max_exact then
            break R;
        end if;
    end for;
end for;

print "";
print "DONE height", height;
print "total_pairs", total_pairs,
      "den7_pairs", den7_pairs,
      "p7_chart_pass", p7_chart_pass,
      "p7_chart_filtered", p7_chart_filtered,
      "aux_pass", aux_pass,
      "aux_filtered", aux_filtered,
      "smooth_tests", smooth_tests,
      "smooth", smooth,
      "yr_fail", yr_fail,
      "exact_tests", exact_tests,
      "Tx_halves", tx_halves,
      "PR_halves", pr_halves,
      "both_halves", both_halves,
      "torsion_tests", torsion_tests,
      "pr_divisibility_errors", pr_divisibility_errors,
      "hits", #hits;

print "CHART_COUNTS";
for key in Sort([k : k in Keys(chart_seen)]) do
    pass := IsDefined(chart_pass_counts, key) select chart_pass_counts[key] else 0;
    exact := IsDefined(chart_exact_counts, key) select chart_exact_counts[key] else 0;
    txc := IsDefined(chart_tx_counts, key) select chart_tx_counts[key] else 0;
    prc := IsDefined(chart_pr_counts, key) select chart_pr_counts[key] else 0;
    bothc := IsDefined(chart_both_counts, key) select chart_both_counts[key] else 0;
    print " ", key, "seen", chart_seen[key], "p7_pass", pass,
          "exact", exact, "Tx_half", txc, "PR_half", prc,
          "both", bothc;
    if IsDefined(chart_samples, key) then
        print "    samples", chart_samples[key];
    end if;
end for;

print "AUX_KILL_COUNTS";
for key in Sort([k : k in Keys(aux_kill_counts)]) do
    print " ", key, aux_kill_counts[key];
end for;

if #hits gt 0 then
    print "HITS";
    for h in hits do
        print h;
    end for;
end if;

quit;
