//////////////////////////////////////////////////////////////////////
//  FRONT 1B: p=7-boundary/CRT-guided rational search for [4,16]
//  inside the M_1(8,4) family.
//
//  The good affine [4,16] locus is empty over F_7.  This script
//  recomputes the cleared p=7 boundary closure for TARGET_416, uses
//  those boundary residues as a CRT gate for rational (R,w), applies
//  optional good-open filters at larger primes, then runs exact
//  Jacobian divisibility tests:
//      T_x = [x,0] divisible by 2, and P_R=(-R,Y_R) divisible by 2.
//
//  Typical small validation run:
//      magma -b height:=35 max_exact_tests:=30 aux_primes:="11,13" \
//          agent_m18_416_search_crt.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned height then
    height := 50;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

if not assigned aux_primes then
    aux_prime_list := [11,13,17];
elif Type(aux_primes) eq MonStgElt then
    aux_prime_list := (#aux_primes eq 0) select [] else
        [StringToInteger(s) : s in Split(aux_primes, ",")];
else
    aux_prime_list := aux_primes;
end if;

if not assigned max_exact_tests then
    max_exact_tests := 50;
elif Type(max_exact_tests) eq MonStgElt then
    max_exact_tests := StringToInteger(max_exact_tests);
end if;

if not assigned max_candidates then
    max_candidates := 0; // 0 means no cap after the CRT/local filters.
elif Type(max_candidates) eq MonStgElt then
    max_candidates := StringToInteger(max_candidates);
end if;

if not assigned max_hits then
    max_hits := 10;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;

if not assigned progress_interval then
    progress_interval := 50;
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

procedure AppendAssoc(~A, key, val)
    if IsDefined(A, key) then
        vals := A[key];
    else
        vals := [];
    end if;
    Append(~vals, val);
    A[key] := vals;
end procedure;

function PairKey(a, b)
    return IntegerToString(Z!a) cat "," cat IntegerToString(Z!b);
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

function BoundaryKey(labels)
    if #labels eq 0 then
        return "none";
    end if;
    return Join(Sort(labels), "&");
end function;

function GoodHyperellipticPolynomial(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
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

function BuildBoundaryClosureEquations416(F)
    Rng<R,w,U,V,M,N,a,b,c,d,e,rho,sigma> := PolynomialRing(F, 13, "grevlex");
    K := FieldOfFractions(Rng);
    PX<X> := PolynomialRing(K);

    t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
    A := X^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*X + R^4;
    B := (R + 2 + 4*t)*X^2 + (R^2 + 4*R + 1 + 8*t)*X
         + (2*R^2 + R + 4*t);
    f := X*A*B;
    h := A*B;
    c4 := Coefficient(h, 4);

    a0 := X^2 + U*X + V;
    F0 := h - X*(M*X+N)^2 - c4*a0^2;
    E0 := [Numerator(Coefficient(F0, i)) : i in [0..4]];

    q416 := X^2 + a*X + b;
    ell416 := c*X^2 + d*X + e;
    F416 := f - ell416^2 - c4*(X+R)*q416^2;
    E416 := [Numerator(Coefficient(F416, i)) : i in [0..4]];
    return E0, E416;
end function;

function AllZero(eqs, vals)
    for pol in eqs do
        if Evaluate(pol, vals) ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function HasFirstClosure(E0, elems, r0, w0)
    z := Parent(r0)!0;
    for UU in elems do
        for VV in elems do
            for MM in elems do
                for NN in elems do
                    vals := [r0,w0,UU,VV,MM,NN,z,z,z,z,z,z,z];
                    if AllZero(E0, vals) then
                        return true, <UU,VV,MM,NN>;
                    end if;
                end for;
            end for;
        end for;
    end for;
    return false, <>;
end function;

function Has416Closure(E416, elems, r0, w0)
    z := Parent(r0)!0;
    for aa in elems do
        for bb in elems do
            for cc in elems do
                for dd in elems do
                    for ee in elems do
                        vals := [r0,w0,z,z,z,z,aa,bb,cc,dd,ee,z,z];
                        if AllZero(E416, vals) then
                            return true, <aa,bb,cc,dd,ee>;
                        end if;
                    end for;
                end for;
            end for;
        end for;
    end for;
    return false, <>;
end function;

function BoundaryClosure416Residues(p)
    F := GF(p);
    elems := [a : a in F];
    E0, E416 := BuildBoundaryClosureEquations416(F);
    allowed := {};
    key_by_pair := AssociativeArray();
    samples_by_key := AssociativeArray();
    total_boundary := 0;
    first_closure := 0;
    target416_closure := 0;

    for r0 in F do
        for w0 in F do
            labels := BoundaryLabelsFinite(F, r0, w0);
            if #labels eq 0 then
                continue;
            end if;
            total_boundary +:= 1;
            has_first, first_sample := HasFirstClosure(E0, elems, r0, w0);
            if has_first then
                first_closure +:= 1;
            end if;
            ok416, sol416 := Has416Closure(E416, elems, r0, w0);
            if has_first and ok416 then
                target416_closure +:= 1;
                rr := Z!r0;
                ww := Z!w0;
                Include(~allowed, <rr,ww>);
                key := BoundaryKey(labels);
                key_by_pair[PairKey(rr,ww)] := key;
                AppendAssoc(~samples_by_key, key, <rr,ww,first_sample,sol416>);
            end if;
        end for;
    end for;
    return allowed, key_by_pair, samples_by_key, total_boundary,
           first_closure, target416_closure;
end function;

function OpenTarget416Residues(p)
    F := GF(p);
    PF<X> := PolynomialRing(F);
    allowed := {};
    good_open := 0;
    first_count := 0;
    pr_count := 0;
    target_count := 0;

    for R in F do
        for w in F do
            if #BoundaryLabelsFinite(F, R, w) ne 0 then
                continue;
            end if;
            f, t, A, B := FamilyPolynomialFinite(F, R, w);
            if not GoodHyperellipticPolynomial(f) then
                continue;
            end if;
            yR := YRFinite(F, R, w);
            if yR^2 ne Evaluate(f, -R) then
                continue;
            end if;

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
                Include(~allowed, <Z!R,Z!w>);
            end if;
        end for;
    end for;
    return allowed, good_open, first_count, pr_count, target_count;
end function;

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

function ResidueInt(q, p)
    den := (Z!Denominator(q)) mod p;
    if den eq 0 then
        return false, 0;
    end if;
    num := (Z!Numerator(q)) mod p;
    return true, (num * InverseMod(den, p)) mod p;
end function;

function IsSquareQ(q)
    q := Q!q;
    if q lt 0 then
        return false;
    end if;
    return IsSquare(Numerator(q)) and IsSquare(Denominator(q));
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function FamilyPolynomial(R, w)
    t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
    A := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
    B := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x
         + (2*R^2 + R + 4*t);
    return x*A*B, t, A, B;
end function;

function YRational(R, w)
    Qfac := R^2 - (Q!1/2)*R*w^2 + (Q!1/2)*R - w^2;
    return -2*R*(R-1)^2*Qfac/(w^2-1);
end function;

function PlusDisc(R, w)
    return -4*(w-R)*(R-1)^2*(R+1)*(R+w)*(w+1)
            *(R*w - 3*R + 3*w - 1);
end function;

function MinusDisc(R, w)
    return 4*(w-1)*(R+1)*(R*w + 3*R + 3*w + 1)
           *(R^4 - 2*R^3 + R^2*w^2 - R^2 + 2*R*w^2 - w^2);
end function;

function FirstCoverPossible(R, w)
    return IsSquareQ(PlusDisc(R,w)) or IsSquareQ(MinusDisc(R,w));
end function;

function TangentCandidates(f, R, w)
    h := ExactQuotient(f, x);
    c1 := Coefficient(h, 1);
    c2 := Coefficient(h, 2);
    c3 := Coefficient(h, 3);
    c4 := Coefficient(h, 4);

    out := [];
    PRing<U> := PolynomialRing(Q);

    for V in [R^2*w, -R^2*w] do
        F := 4*(c3 - 2*c4*U)*(c1 - 2*c4*U*V)
             - (c2 - c4*(U^2 + 2*V))^2;
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

function Has416(invs)
    vals := Sort([Valuation(n, 2) : n in invs]);
    vals := Reverse(vals);
    return #vals ge 2 and vals[1] ge 4 and vals[2] ge 2;
end function;

function PassesAuxiliaryFilters(R, w, aux_prime_list, open_target_by_prime)
    for p in aux_prime_list do
        okR, rp := ResidueInt(R, p);
        okW, wp := ResidueInt(w, p);
        if not (okR and okW) then
            continue;
        end if;
        F := GF(p);
        labels := BoundaryLabelsFinite(F, F!rp, F!wp);
        if #labels ne 0 then
            continue;
        end if;
        if <rp,wp> notin open_target_by_prime[p] then
            return false, p;
        end if;
    end for;
    return true, 0;
end function;

procedure PrintAssocCounts(label, A)
    print label;
    for key in Sort([k : k in Keys(A)]) do
        print " ", key, A[key];
    end for;
end procedure;

print "M_1(8,4) [4,16] p=7-boundary CRT search";
print "height", height, "aux_primes", aux_prime_list,
      "max_exact_tests", max_exact_tests,
      "max_candidates", max_candidates, "max_hits", max_hits;

p7_allowed, p7_key_by_pair, p7_samples_by_key, p7_boundary_total,
    p7_first_closure, p7_target416_closure := BoundaryClosure416Residues(7);

print "p7_boundary_closure", "boundary", p7_boundary_total,
      "first_closure", p7_first_closure,
      "target416_closure", p7_target416_closure,
      "allowed_pairs", #p7_allowed;
print "p7_allowed_pairs", Sort(Setseq(p7_allowed));
print "p7_allowed_strata";
for key in Sort([k : k in Keys(p7_samples_by_key)]) do
    print " ", key, p7_samples_by_key[key];
end for;

open_target_by_prime := AssociativeArray();
for p in aux_prime_list do
    allowed, good_open, first_count, pr_count, target_count :=
        OpenTarget416Residues(p);
    open_target_by_prime[p] := allowed;
    print "aux_prime_profile", p,
          "good_open", good_open,
          "first_Tx_half", first_count,
          "PR_half", pr_count,
          "target416", target_count;
end for;

params := RationalParametersOfHeight(height);
p7_buckets := AssociativeArray();
p7_denominator_bad := 0;
for r in params do
    ok, rv := ResidueInt(r, 7);
    if ok then
        AppendAssoc(~p7_buckets, rv, r);
    else
        p7_denominator_bad +:= 1;
    end if;
end for;

p7_affine_param_count := #params - p7_denominator_bad;
p7_crt_pair_budget := 0;
for pair in p7_allowed do
    rr := pair[1];
    ww := pair[2];
    rcount := IsDefined(p7_buckets, rr) select #p7_buckets[rr] else 0;
    wcount := IsDefined(p7_buckets, ww) select #p7_buckets[ww] else 0;
    p7_crt_pair_budget +:= rcount*wcount;
end for;

print "parameter_count", #params,
      "p7_affine_parameters", p7_affine_param_count,
      "p7_denominator_bad_parameters", p7_denominator_bad,
      "p7_crt_pair_budget", p7_crt_pair_budget;
print "NOTE denominator-divisible p=7 parameters are not enumerated in this affine CRT pass.";

candidates_after_aux := 0;
aux_kill_counts := AssociativeArray();
p7_stratum_candidate_counts := AssociativeArray();
family_smooth := 0;
first_possible := 0;
tangent_bases := 0;
exact_tests := 0;
first_verified := 0;
pr_halved := 0;
torsion_tests := 0;
hits := [];
torsion_counts := AssociativeArray();

stop := false;
allowed_seq := Sort(Setseq(p7_allowed));

for pair in allowed_seq do
    if stop then break; end if;
    rr := pair[1];
    ww := pair[2];
    if not (IsDefined(p7_buckets, rr) and IsDefined(p7_buckets, ww)) then
        continue;
    end if;
    p7_pair_key := PairKey(rr, ww);
    p7_stratum := p7_key_by_pair[p7_pair_key];

    for R in p7_buckets[rr] do
        if stop then break; end if;
        for w in p7_buckets[ww] do
            if stop then break; end if;

            ok_aux, kill_p := PassesAuxiliaryFilters(R, w, aux_prime_list,
                                                      open_target_by_prime);
            if not ok_aux then
                Increment(~aux_kill_counts, kill_p);
                continue;
            end if;

            candidates_after_aux +:= 1;
            Increment(~p7_stratum_candidate_counts, p7_stratum);
            if progress_interval gt 0 and
               candidates_after_aux mod progress_interval eq 0 then
                print "progress", "candidates_after_aux", candidates_after_aux,
                      "smooth", family_smooth,
                      "first_possible", first_possible,
                      "tangent_bases", tangent_bases,
                      "exact_tests", exact_tests,
                      "first_verified", first_verified,
                      "pr_halved", pr_halved,
                      "hits", #hits;
            end if;
            if max_candidates gt 0 and candidates_after_aux ge max_candidates then
                stop := true;
            end if;

            if R eq 0 or w in {Q!-1, Q!0, Q!1} then
                continue;
            end if;

            f, t, A, B := FamilyPolynomial(R, w);
            if not GoodHyperellipticPolynomial(f) then
                continue;
            end if;
            family_smooth +:= 1;

            if not FirstCoverPossible(R, w) then
                continue;
            end if;
            first_possible +:= 1;

            candidates := TangentCandidates(f, R, w);
            if #candidates eq 0 then
                continue;
            end if;
            tangent_bases +:= 1;

            exact_tests +:= 1;
            fI, L := IntegralModelPolynomial(f);
            C := HyperellipticCurve(fI);
            J := Jacobian(C);
            Tx := J![x, Q!0];
            tx_divisible, tx_half := IsDivisibleBy(Tx, 2);
            if not tx_divisible then
                print "WARNING tangent candidate failed exact Tx half",
                      "R", R, "w", w, "t", t;
                if max_exact_tests gt 0 and exact_tests ge max_exact_tests then
                    stop := true;
                end if;
                continue;
            end if;
            first_verified +:= 1;

            yR := YRational(R, w);
            yRI := L*yR;
            if Evaluate(fI, -R) ne yRI^2 then
                print "WARNING P_R formula failed over Q",
                      "R", R, "w", w, "t", t;
                if max_exact_tests gt 0 and exact_tests ge max_exact_tests then
                    stop := true;
                end if;
                continue;
            end if;

            PR := J![x + R, yRI];
            pr_divisible, pr_half := IsDivisibleBy(PR, 2);
            if pr_divisible then
                pr_halved +:= 1;
                torsion_tests +:= 1;
                G, phi := TorsionSubgroup(J);
                invs := Invariants(G);
                key := Sprint(invs);
                Increment(~torsion_counts, key);
                ord := TorsionOrder(invs);
                exp := Exponent(invs);
                print "PR_HALF", "R", R, "w", w, "t", t,
                      "p7_pair", pair, "p7_stratum", p7_stratum,
                      "tx_half_order", Order(tx_half),
                      "pr_half_order", Order(pr_half),
                      "torsion", invs, "order", ord, "exponent", exp;
                print "  f =", fI;
                if Has416(invs) then
                    Append(~hits, <R,w,t,invs,fI>);
                    print "HIT_416", "R", R, "w", w, "t", t,
                          "torsion", invs;
                end if;
            end if;

            if (#hits ge max_hits and max_hits gt 0) or
               (max_exact_tests gt 0 and exact_tests ge max_exact_tests) then
                stop := true;
            end if;
        end for;
    end for;
end for;

print "DONE";
print "height", height, "aux_primes", aux_prime_list;
print "candidates_after_aux", candidates_after_aux,
      "family_smooth", family_smooth,
      "first_possible", first_possible,
      "tangent_bases", tangent_bases,
      "exact_tests", exact_tests,
      "first_verified", first_verified,
      "pr_halved", pr_halved,
      "torsion_tests", torsion_tests,
      "hits", #hits;
PrintAssocCounts("AUX_KILL_COUNTS", aux_kill_counts);
PrintAssocCounts("P7_STRATUM_CANDIDATE_COUNTS", p7_stratum_candidate_counts);
PrintAssocCounts("TORSION_COUNTS", torsion_counts);

quit;
