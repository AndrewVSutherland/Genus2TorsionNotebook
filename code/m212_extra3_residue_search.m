//////////////////////////////////////////////////////////////////////
//  Fast residue-table search for M(2,12) plus independent 3-torsion.
//
//  The M(2,12) chart is
//      a = (1-z^2)/(4*(r+1)).
//
//  The marked order-12 class supplies one rational 3-direction.  An
//  independent rational 3-torsion point forces J(F_p)[3] to have rank at
//  least 2 at every good prime p != 3.
//
//  This script precomputes, for each p, the good residues (z,r) with
//  3-rank >= 2 and the bad/boundary residues.  The rational search then
//  only does residue lookups before exact torsion.
//
//  Typical run:
//      magma -b height:=50 prime_bound:=43 max_exact:=200 \
//          code/m212_extra3_residue_search.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned height then
    height := 50;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned prime_bound then
    prime_bound := 43;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;
if not assigned max_exact then
    max_exact := 200;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned max_print then
    max_print := 50;
elif Type(max_print) eq MonStgElt then
    max_print := StringToInteger(max_print);
end if;
if not assigned progress_interval then
    progress_interval := 1000000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Q := Rationals();
Z := Integers();
Qx<x> := PolynomialRing(Q);

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            q := Q!num / Q!den;
            key := Sprint(q);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, q);
            end if;
        end for;
    end for;
    return vals;
end function;

function M212PolynomialFromZR(z, r)
    if r eq -1 or z^2 eq 1 then
        return false, Qx!0, Q!0;
    end if;
    a := (1-z^2)/(4*(r+1));
    if a eq 0 then
        return false, Qx!0, Q!0;
    end if;
    T := a*x^2 - x + r;
    h := (x-r)*(T+1);
    W := h^2 + 4*a*x^2*T*(T+1);
    return true, W, a;
end function;

function M212PolynomialFinite(K, z, r)
    P<X> := PolynomialRing(K);
    if r eq -K!1 or z^2 eq K!1 then
        return false, P!0, K!0;
    end if;
    a := (1-z^2)/(4*(r+1));
    if a eq 0 then
        return false, P!0, K!0;
    end if;
    T := a*X^2 - X + r;
    h := (X-r)*(T+1);
    W := h^2 + 4*a*X^2*T*(T+1);
    return true, W, a;
end function;

function GoodHyperellipticPolynomial(f)
    return Degree(f) in {5,6} and Discriminant(f) ne 0;
end function;

function ThreeRankFromInvariants(invs)
    return #[n : n in invs | (Z!n) mod 3 eq 0];
end function;

function ThreeRankFinitePolynomial(f)
    C := HyperellipticCurve(f);
    A, phi := AbelianGroup(Jacobian(C));
    invs := Invariants(A);
    return ThreeRankFromInvariants(invs), invs;
end function;

function PairKey(a, b)
    return <Z!a, Z!b>;
end function;

function ResidueTables(p)
    K := GF(p);
    allowed := {};
    bad := {};
    rank_counts := AssociativeArray();

    for zz in K do
        for rr in K do
            key := PairKey(Z!zz, Z!rr);
            ok, W, a := M212PolynomialFinite(K, zz, rr);
            if not ok or not GoodHyperellipticPolynomial(W) then
                Include(~bad, key);
                continue;
            end if;
            trank, invs := ThreeRankFinitePolynomial(W);
            if IsDefined(rank_counts, trank) then
                rank_counts[trank] +:= 1;
            else
                rank_counts[trank] := 1;
            end if;
            if trank ge 2 then
                Include(~allowed, key);
            end if;
        end for;
    end for;

    return allowed, bad, rank_counts;
end function;

function ResidueOfRational(q, p)
    num := Numerator(q);
    den := Denominator(q);
    if (den mod p) eq 0 then
        return false, 0;
    end if;
    K := GF(p);
    return true, Z!(K!num / K!den);
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return Qx!(L^2*f), L;
end function;

function TorsionOrder(invs)
    if #invs eq 0 then
        return 1;
    end if;
    return &*invs;
end function;

function TorsionThreeRank(invs)
    return ThreeRankFromInvariants(invs);
end function;

primes := [p : p in PrimesUpTo(prime_bound) | p notin {2,3}];
allowed_by_p := AssociativeArray();
bad_by_p := AssociativeArray();

print "M(2,12)+independent 3 residue-table search";
print "height", height, "prime_bound", prime_bound, "primes", primes,
      "max_exact", max_exact;

for p in primes do
    allowed, bad, rank_counts := ResidueTables(p);
    allowed_by_p[p] := allowed;
    bad_by_p[p] := bad;
    print "prime", p, "allowed_rank_ge2", #allowed, "bad", #bad,
          "rank_counts", Sort([<k, rank_counts[k]> : k in Keys(rank_counts)]);
end for;

params := RationalParametersOfHeight(height);
checked := 0;
smooth := 0;
residue_survivors := 0;
exact_tests := 0;
hits := [];
printed := 0;
seen_ar := {};
kill_counts := AssociativeArray();

for z in params do
    for r in params do
        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "progress", checked, "smooth", smooth,
                  "residue_survivors", residue_survivors,
                  "exact", exact_tests, "hits", #hits;
        end if;

        ok, W, a := M212PolynomialFromZR(z, r);
        if not ok or not GoodHyperellipticPolynomial(W) then
            continue;
        end if;
        smooth +:= 1;

        ar_key := Sprint(<a,r>);
        if ar_key in seen_ar then
            continue;
        end if;
        Include(~seen_ar, ar_key);

        killed := false;
        bad_primes := [];
        for p in primes do
            okz, rz := ResidueOfRational(z, p);
            okr, rr := ResidueOfRational(r, p);
            if not okz or not okr then
                Append(~bad_primes, p);
                continue;
            end if;
            key := PairKey(rz, rr);
            if key in bad_by_p[p] then
                Append(~bad_primes, p);
                continue;
            end if;
            if key notin allowed_by_p[p] then
                killed := true;
                if IsDefined(kill_counts, p) then
                    kill_counts[p] +:= 1;
                else
                    kill_counts[p] := 1;
                end if;
                break;
            end if;
        end for;
        if killed then
            continue;
        end if;

        residue_survivors +:= 1;
        if printed lt max_print then
            print "SURVIVOR", "z", z, "r", r, "a", a,
                  "bad_primes", bad_primes;
            printed +:= 1;
        end if;

        if exact_tests ge max_exact then
            continue;
        end if;

        WI, L := IntegralModelPolynomial(W);
        C := HyperellipticCurve(WI);
        G, phi := TorsionSubgroup(Jacobian(C));
        invs := Invariants(G);
        exact_tests +:= 1;
        print "EXACT", "z", z, "r", r, "a", a,
              "torsion", invs, "order", TorsionOrder(invs),
              "three_rank", TorsionThreeRank(invs);
        if TorsionThreeRank(invs) ge 2 then
            Append(~hits, <z,r,a,invs,WI>);
            print "EXTRA3_HIT", "z", z, "r", r, "a", a,
                  "torsion", invs, "W", WI;
        end if;
    end for;
end for;

print "DONE M(2,12)+independent 3 residue-table search";
print "checked", checked, "smooth", smooth,
      "unique_ar", #seen_ar,
      "residue_survivors", residue_survivors,
      "exact_tests", exact_tests, "hits", #hits;
print "kill_counts";
for p in Sort([k : k in Keys(kill_counts)]) do
    print p, kill_counts[p];
end for;

quit;
