//////////////////////////////////////////////////////////////////////
//  CRT-guided exact search for [8,8] on the reduced square subcover.
//
//  This searches the good-open rational square branch w=s^2.  For a
//  chosen list of finite primes, it first computes the reduced [8,8]
//  allowed residue pairs (R,s) mod p.  It then buckets rational
//  parameters of bounded height by their residue vector and only tests
//  exact pairs satisfying all chosen congruence conditions.
//
//  This deliberately excludes rational parameters whose denominators are
//  divisible by one of the CRT primes.  Those are boundary cases at that
//  prime and should be handled by a separate boundary search.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned height then
    height := 200;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

if not assigned primes then
    prime_list := [7,11,13,17,19,23];
elif Type(primes) eq MonStgElt then
    prime_list := [StringToInteger(s) : s in Split(primes, ",")];
else
    prime_list := primes;
end if;

if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;

if not assigned max_candidates then
    max_candidates := 0; // 0 means no cap.
elif Type(max_candidates) eq MonStgElt then
    max_candidates := StringToInteger(max_candidates);
end if;

if not assigned progress_interval then
    progress_interval := 1000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

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

function BoundaryBad(R, w)
    if R eq 0 or w eq 0 or w eq 1 or w eq -1 then
        return true;
    end if;
    if R eq 1 or R eq -1 or R eq w or R eq -w then
        return true;
    end if;
    if R*w - 3*R + 3*w - 1 eq 0 then
        return true;
    end if;
    if R*w + 3*R + 3*w + 1 eq 0 then
        return true;
    end if;
    if 2*R^2 - R*w^2 + R - 2*w^2 eq 0 then
        return true;
    end if;
    if R^4 - 2*R^3 + R^2*w^2 - R^2 + 2*R*w^2 - w^2 eq 0 then
        return true;
    end if;
    return false;
end function;

function ResidueInt(q, p)
    den := Z!Denominator(q);
    if den mod p eq 0 then
        return false, 0;
    end if;
    num := Z!Numerator(q);
    return true, ((num mod p) * InverseMod(den mod p, p)) mod p;
end function;

function ResidueVector(q, prime_list)
    vals := [];
    for p in prime_list do
        ok, rp := ResidueInt(q, p);
        if not ok then
            return false, [];
        end if;
        Append(~vals, rp);
    end for;
    return true, vals;
end function;

function Key(vals)
    out := "";
    for i in [1..#vals] do
        if i gt 1 then
            out cat:= ",";
        end if;
        out cat:= IntegerToString(vals[i]);
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

function Has88(invs)
    vals := Sort([Valuation(n, 2) : n in invs]);
    vals := Reverse(vals);
    return #vals ge 2 and vals[1] ge 3 and vals[2] ge 3;
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

function FamilyDataFinite(F, R, s)
    PF<X> := PolynomialRing(F);
    w := s^2;
    t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
    A := X^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*X + R^4;
    B := (R + 2 + 4*t)*X^2 + (R^2 + 4*R + 1 + 8*t)*X
         + (2*R^2 + R + 4*t);
    h := A*B;
    return X*h, h, Coefficient(h, 1), Coefficient(h, 2),
           Coefficient(h, 3), Coefficient(h, 4);
end function;

function SquareRootsFinite(F, a)
    PF<T> := PolynomialRing(F);
    return [rt[1] : rt in Roots(T^2 - a)];
end function;

function FirstCoverSolutionsFinite(F, R, s)
    f, h, c1, c2, c3, c4 := FamilyDataFinite(F, R, s);
    PF<U> := PolynomialRing(F);
    V := R^2*s^2;
    FU := 4*(c3 - 2*c4*U)*(c1 - 2*c4*U*V)
          - (c2 - c4*(U^2 + 2*V))^2;

    sols := [];
    for rt in Roots(FU) do
        U0 := rt[1];
        M2 := c3 - 2*c4*U0;
        N2 := c1 - 2*c4*U0*V;
        C := c2 - c4*(U0^2 + 2*V);
        for M0 in SquareRootsFinite(F, M2) do
            if M0 ne 0 then
                N0 := C/(2*M0);
                if N0^2 eq N2 then
                    Append(~sols, <U0, M0, N0, c4>);
                end if;
            else
                if C ne 0 then
                    continue;
                end if;
                for N0 in SquareRootsFinite(F, N2) do
                    Append(~sols, <U0, M0, N0, c4>);
                end for;
            end if;
        end for;
    end for;
    return sols;
end function;

function TauSolvesFinite(F, A, B, C)
    if A eq 0 then
        if B eq 0 then
            return C eq 0, F!0;
        end if;
        return true, -C/B;
    end if;
    PF<T> := PolynomialRing(F);
    rts := Roots(A*T^2 + B*T + C);
    if #rts eq 0 then
        return false, _;
    end if;
    return true, rts[1][1];
end function;

function ReducedSecondSolvableFinite(F, R, s, sol)
    U := sol[1];
    M := sol[2];
    N := sol[3];
    c := sol[4];
    V := R^2*s^2;

    for eta in [F!1, -F!1] do
        mu := eta*R*s;
        A := 2*mu - U;
        for z in F do
            if z eq 0 then
                continue;
            end if;

            D := A - c*z^2;
            Pz := c^2*z^4 - 4*M*c*z^3
                  + (4*M^2 + 2*U*c)*z^2
                  + (4*M*U - 8*N)*z + U^2 - 4*V;
            B := 2*z*(M*mu - N) - mu*A - c*mu*z^2;
            C := V*c*z^2;

            if D ne 0 then
                tau := -Pz/(4*D);
                if A*tau^2 + B*tau + C eq 0 then
                    return true;
                end if;
            elif Pz eq 0 then
                ok, tau := TauSolvesFinite(F, A, B, C);
                if ok then
                    return true;
                end if;
            end if;
        end for;
    end for;
    return false;
end function;

function AllowedSByRForPrime(p)
    F := GF(p);
    allowed := [* [Z | ] : i in [1..p] *];
    open_good := 0;
    first_pairs := 0;
    reduced_pairs := 0;

    for R in F do
        for s in F do
            w := s^2;
            if #BoundaryLabelsFinite(F, R, w) ne 0 then
                continue;
            end if;
            f, h, c1, c2, c3, c4 := FamilyDataFinite(F, R, s);
            if not GoodHyperellipticPolynomial(f) then
                continue;
            end if;
            open_good +:= 1;
            fs := FirstCoverSolutionsFinite(F, R, s);
            if #fs eq 0 then
                continue;
            end if;
            first_pairs +:= 1;
            if &or [ReducedSecondSolvableFinite(F, R, s, sol) : sol in fs] then
                reduced_pairs +:= 1;
                Append(~allowed[Z!R + 1], Z!s);
            end if;
        end for;
    end for;

    return allowed, open_good, first_pairs, reduced_pairs;
end function;

function CompatibleSKeys(rvec, allowed_by_prime)
    keys := [""];
    for i in [1..#rvec] do
        opts := allowed_by_prime[i][rvec[i] + 1];
        if #opts eq 0 then
            return [];
        end if;
        new_keys := [];
        for key in keys do
            for sres in opts do
                if key eq "" then
                    Append(~new_keys, IntegerToString(sres));
                else
                    Append(~new_keys, key cat "," cat IntegerToString(sres));
                end if;
            end for;
        end for;
        keys := new_keys;
    end for;
    return keys;
end function;

function FirstCoverSolutionsExact(f, R, s)
    V := R^2*s^2;
    h := ExactQuotient(f, x);
    c1 := Coefficient(h, 1);
    c2 := Coefficient(h, 2);
    c3 := Coefficient(h, 3);
    c4 := Coefficient(h, 4);

    PR<Uvar> := PolynomialRing(Q);
    FU := 4*(c3 - 2*c4*Uvar)*(c1 - 2*c4*Uvar*V)
          - (c2 - c4*(Uvar^2 + 2*V))^2;

    sols := [];
    for rt in Roots(FU) do
        U0 := rt[1];
        M2 := c3 - 2*c4*U0;
        N2 := c1 - 2*c4*U0*V;
        for M0 in SquareRootsQ(M2) do
            for N0 in SquareRootsQ(N2) do
                a0 := x^2 + U0*x + V;
                identity := h - x*(M0*x+N0)^2 - c4*a0^2;
                if identity eq 0 then
                    Append(~sols, <U0, V, M0, N0, c4>);
                end if;
            end for;
        end for;
    end for;
    return sols;
end function;

function TauRootsQ(A, B, C)
    QT<T> := PolynomialRing(Q);
    return [rt[1] : rt in Roots(A*T^2 + B*T + C)];
end function;

function ReducedSecondSolutionsExact(R, s, U, V, M, N, c)
    out := [];
    S<z> := PolynomialRing(Q);
    for eta in [Q!1, Q!-1] do
        mu := eta*R*s;
        A := 2*mu - U;
        D := A - c*z^2;
        Pz := c^2*z^4 - 4*M*c*z^3
              + (4*M^2 + 2*U*c)*z^2
              + (4*M*U - 8*N)*z + U^2 - 4*V;
        B := 2*z*(M*mu - N) - mu*A - c*mu*z^2;
        C := V*c*z^2;
        G := A*Pz^2 - 4*B*Pz*D + 16*C*D^2;

        if G eq 0 then
            print "WARNING zero exact z-eliminant",
                  "R", R, "s", s, "U", U, "M", M, "N", N, "eta", eta;
            continue;
        end if;

        for rt in Roots(G) do
            z0 := rt[1];
            if z0 eq 0 then
                continue;
            end if;
            Dz := Evaluate(D, z0);
            Pz0 := Evaluate(Pz, z0);
            Bz := Evaluate(B, z0);
            Cz := Evaluate(C, z0);
            taus := [];

            if Dz ne 0 then
                tau := -Pz0/(4*Dz);
                if A*tau^2 + Bz*tau + Cz eq 0 then
                    taus := [tau];
                end if;
            elif Pz0 eq 0 then
                taus := TauRootsQ(A, Bz, Cz);
            end if;

            for tau in taus do
                a_val := U/2 + M*z0 + tau - (c/2)*z0^2;
                b_val := eta*R*s*tau;
                Append(~out, <z0, tau, eta, a_val, b_val>);
            end for;
        end for;
    end for;
    return out;
end function;

function IrreducibleFrobeniusCertificate(fI)
    C := HyperellipticCurve(fI);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79] do
        try
            fp := ChangeRing(fI, GF(p));
            if not GoodHyperellipticPolynomial(fp) then
                continue;
            end if;
            Lp := LPolynomial(ChangeRing(C, GF(p)));
            fac := Factorization(Lp);
            if #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4 then
                return true, p, Lp;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false, 0, P!0;
end function;

print "CRT-guided reduced [8,8] exact search";
print "height", height, "primes", prime_list,
      "max_hits", max_hits, "max_candidates", max_candidates,
      "progress_interval", progress_interval;

allowed_by_prime := [* *];
for p in prime_list do
    allowed, open_good_p, first_pairs_p, reduced_pairs_p := AllowedSByRForPrime(p);
    Append(~allowed_by_prime, allowed);
    print "prime_profile", p,
          "open_good", open_good_p,
          "first_pairs", first_pairs_p,
          "reduced_pairs", reduced_pairs_p;
end for;

params := RationalParametersOfHeight(height);
s_buckets := AssociativeArray();
r_data := [];
denom_bad := 0;

for q in params do
    ok, vec := ResidueVector(q, prime_list);
    if not ok then
        denom_bad +:= 1;
        continue;
    end if;
    Append(~r_data, <q, vec>);
    key := Key(vec);
    if IsDefined(s_buckets, key) then
        lst := s_buckets[key];
        Append(~lst, q);
        s_buckets[key] := lst;
    else
        s_buckets[key] := [q];
    end if;
end for;

print "parameter_count", #params,
      "strict_good_parameters", #r_data,
      "denominator_bad_parameters", denom_bad,
      "s_residue_buckets", #Keys(s_buckets);

pair_candidates := 0;
open_good := 0;
first_bases := 0;
first_solutions := 0;
second_solutions := 0;
exact_verified := 0;
torsion_tests := 0;
hits := [];
torsion_counts := AssociativeArray();

for rd in r_data do
    if #hits ge max_hits then
        break;
    end if;
    R := rd[1];
    rvec := rd[2];
    skeys := CompatibleSKeys(rvec, allowed_by_prime);

    for skey in skeys do
        if #hits ge max_hits then
            break;
        end if;
        if not IsDefined(s_buckets, skey) then
            continue;
        end if;
        for s in s_buckets[skey] do
            if #hits ge max_hits then
                break;
            end if;
            pair_candidates +:= 1;
            if progress_interval gt 0 and pair_candidates mod progress_interval eq 0 then
                print "progress", pair_candidates,
                      "open_good", open_good,
                      "first_bases", first_bases,
                      "first_solutions", first_solutions,
                      "second_solutions", second_solutions,
                      "verified", exact_verified,
                      "hits", #hits;
            end if;
            if max_candidates gt 0 and pair_candidates gt max_candidates then
                break rd;
            end if;

            w := s^2;
            if BoundaryBad(R, w) then
                continue;
            end if;
            f, t, A, B := FamilyPolynomial(R, w);
            if not GoodHyperellipticPolynomial(f) then
                continue;
            end if;
            open_good +:= 1;

            first := FirstCoverSolutionsExact(f, R, s);
            if #first eq 0 then
                continue;
            end if;
            first_bases +:= 1;
            first_solutions +:= #first;

            fI, L := IntegralModelPolynomial(f);
            Cc := HyperellipticCurve(fI);
            J := Jacobian(Cc);
            Tx := J![x, Q!0];

            for fs in first do
                U0 := fs[1];
                V := fs[2];
                M0 := fs[3];
                N0 := fs[4];
                c4 := fs[5];
                seconds := ReducedSecondSolutionsExact(R, s, U0, V, M0, N0, c4);
                second_solutions +:= #seconds;
                if #seconds eq 0 then
                    continue;
                end if;

                a0 := x^2 + U0*x + V;
                v0 := (M0*U0 - N0)*x + M0*V;
                Hx := J![a0, Q!(L*v0)];
                ok_half, Q16 := IsDivisibleBy(Hx, 2);
                if not ok_half then
                    print "WARNING reduced equations did not verify exact halving",
                          "R", R, "s", s, "U", U0, "M", M0, "N", N0,
                          "seconds", seconds;
                    continue;
                end if;
                exact_verified +:= 1;
                if 2*Hx ne Tx then
                    print "WARNING Hx is not half of Tx", "R", R, "s", s;
                    continue;
                end if;

                torsion_tests +:= 1;
                Gtors, phi := TorsionSubgroup(J);
                invs := Invariants(Gtors);
                key := Sprint(invs);
                if IsDefined(torsion_counts, key) then
                    torsion_counts[key] +:= 1;
                else
                    torsion_counts[key] := 1;
                end if;

                if Has88(invs) then
                    simple, pcert, Lp := IrreducibleFrobeniusCertificate(fI);
                    Append(~hits, <R,s,w,t,invs,simple,pcert,Lp,fI,fs,seconds[1]>);
                    print "HIT", "R", R, "s", s, "w", w, "t", t,
                          "torsion", invs, "order", TorsionOrder(invs),
                          "exponent", Exponent(invs),
                          "simple", simple, "pcert", pcert;
                    print "  first", fs;
                    print "  second", seconds[1];
                    print "  fI", fI;
                else
                    print "NON_88_VERIFIED", "R", R, "s", s, "w", w,
                          "torsion", invs;
                end if;
            end for;
        end for;
    end for;
end for;

print "DONE height", height, "primes", prime_list;
print "pair_candidates", pair_candidates,
      "open_good", open_good,
      "first_bases", first_bases,
      "first_solutions", first_solutions,
      "second_solutions", second_solutions,
      "exact_verified", exact_verified,
      "torsion_tests", torsion_tests,
      "hits", #hits;
print "TORSION_COUNTS";
for key in Sort([k : k in Keys(torsion_counts)]) do
    print " ", key, torsion_counts[key];
end for;

quit;
