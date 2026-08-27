//////////////////////////////////////////////////////////////////////
//  Search for simple [3,12] by halving the marked order-6 class on
//  the contact-6 [3,6] family.
//
//  The halving cover of the marked contact-6 class is the M(2,12)
//  chart
//
//      a = (1-z^2)/(4*(r+1)).
//
//  In the standard M(12) model
//
//      y^2 + (x-r)(T+1)y = a*x^2*T*(T+1),
//      T = a*x^2 - x + r,
//
//  the completed-square sextic is
//
//      W = (x-r)^2*(T+1)^2 + 4*a*x^2*T*(T+1).
//
//  The marked class has order 12.  To get [3,12], an independent
//  rational 3-torsion direction must be present.  At every good
//  prime p != 2,3, this forces J(F_p) to contain a subgroup compatible
//  with [3,12].  We precompute those residue classes, then exact-test
//  only rational survivors which also pass a geometric simplicity
//  certificate from an irreducible L-polynomial.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned height then
    height := 30;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned prime_bound then
    prime_bound := 43;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;
if not assigned max_exact then
    max_exact := 100;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned max_hits then
    max_hits := 10;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;
if not assigned max_print then
    max_print := 30;
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
            q := Q!num/Q!den;
            key := Sprint(q);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, q);
            end if;
        end for;
    end for;
    return vals;
end function;

function Has312Invariants(invs)
    return #[n : n in invs | (Z!n) mod 3 eq 0] ge 2
       and #[n : n in invs | (Z!n) mod 12 eq 0] ge 1;
end function;

function M212PolynomialFromZR(z, r)
    if r eq -1 or z^2 eq 1 then
        return false, Qx!0, Qx!0, Qx!0, Q!0;
    end if;
    a := (1-z^2)/(4*(r+1));
    if a eq 0 then
        return false, Qx!0, Qx!0, Qx!0, Q!0;
    end if;
    T := a*x^2 - x + r;
    h := (x-r)*(T+1);
    W := h^2 + 4*a*x^2*T*(T+1);
    return true, W, T, h, a;
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

function InvariantsFinitePolynomial(f)
    C := HyperellipticCurve(f);
    A, phi := AbelianGroup(Jacobian(C));
    return Invariants(A);
end function;

function PairKey(a, b)
    return <Z!a, Z!b>;
end function;

function ResidueTables(p)
    K := GF(p);
    allowed := {};
    bad := {};
    rank_counts := AssociativeArray();
    allowed_count := 0;

    for zz in K do
        for rr in K do
            key := PairKey(Z!zz, Z!rr);
            ok, W, a := M212PolynomialFinite(K, zz, rr);
            if not ok or not GoodHyperellipticPolynomial(W) then
                Include(~bad, key);
                continue;
            end if;
            invs := InvariantsFinitePolynomial(W);
            trank := #[n : n in invs | (Z!n) mod 3 eq 0];
            if IsDefined(rank_counts, trank) then
                rank_counts[trank] +:= 1;
            else
                rank_counts[trank] := 1;
            end if;
            if Has312Invariants(invs) then
                Include(~allowed, key);
                allowed_count +:= 1;
            end if;
        end for;
    end for;

    return allowed, bad, rank_counts, allowed_count;
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

function SimpleCertificate(f)
    C := HyperellipticCurve(f);
    for p in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97] do
        try
            fp := ChangeRing(f, GF(p));
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
    return false, 0, Qx!0;
end function;

function OddQuinticAtRoot(W, w)
    out := Qx!0;
    for i in [0..Degree(W)] do
        for j in [0..i] do
            out +:= Coefficient(W, i)*Binomial(i,j)*w^(i-j)*x^(6-j);
        end for;
    end for;
    return out;
end function;

primes := [p : p in PrimesUpTo(prime_bound) | p notin {2,3}];
allowed_by_p := AssociativeArray();
bad_by_p := AssociativeArray();

print "Contact-6 [3,6] -> halve marked 6-class -> simple [3,12] search";
print "height", height, "prime_bound", prime_bound, "primes", primes,
      "max_exact", max_exact, "max_hits", max_hits;

for p in primes do
    allowed, bad, rank_counts, allowed_count := ResidueTables(p);
    allowed_by_p[p] := allowed;
    bad_by_p[p] := bad;
    print "prime", p, "allowed_312", allowed_count, "bad", #bad,
          "rank_counts", Sort([<k, rank_counts[k]> : k in Keys(rank_counts)]);
end for;

params := RationalParametersOfHeight(height);
checked := 0;
smooth := 0;
residue_survivors := 0;
simple_survivors := 0;
exact_tests := 0;
hits := [];
printed := 0;
seen_ar := {};
kill_counts := AssociativeArray();

for z in params do
    for r in params do
        if #hits ge max_hits then
            break z;
        end if;

        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "progress", checked, "smooth", smooth,
                  "residue_survivors", residue_survivors,
                  "simple_survivors", simple_survivors,
                  "exact", exact_tests, "hits", #hits;
        end if;

        ok, W, T, h, a := M212PolynomialFromZR(z, r);
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
        WI, scaleW := IntegralModelPolynomial(W);
        simple, pcert, Lp := SimpleCertificate(WI);
        if not simple then
            if printed lt max_print then
                print "SURVIVOR_NONSIMPLE_CERT_UNKNOWN", "z", z, "r", r,
                      "a", a, "bad_primes", bad_primes;
                printed +:= 1;
            end if;
            continue;
        end if;

        simple_survivors +:= 1;
        if printed lt max_print then
            print "SURVIVOR_SIMPLE", "z", z, "r", r, "a", a,
                  "pcert", pcert, "bad_primes", bad_primes;
            printed +:= 1;
        end if;

        if exact_tests ge max_exact then
            continue;
        end if;

        C := HyperellipticCurve(WI);
        G, phi := TorsionSubgroup(Jacobian(C));
        invs := Invariants(G);
        exact_tests +:= 1;
        print "EXACT", "z", z, "r", r, "a", a,
              "torsion", invs, "order", TorsionOrder(invs),
              "pcert", pcert;

        if Has312Invariants(invs) then
            w := 2*(r+1)/(1+z);
            f5 := OddQuinticAtRoot(W, w);
            f5I, scale5 := IntegralModelPolynomial(f5);
            Append(~hits, <z,r,a,invs,WI,f5I,pcert,Lp>);
            print "HIT312_SIMPLE", "z", z, "r", r, "a", a,
                  "torsion", invs;
            print " W_integral", WI;
            print " odd_f5_integral", f5I;
            print " simplicity_prime", pcert, "L", Lp;
        end if;
    end for;
end for;

print "DONE contact6 halveD [3,12] search";
print "checked", checked, "smooth", smooth, "unique_ar", #seen_ar,
      "residue_survivors", residue_survivors,
      "simple_survivors", simple_survivors,
      "exact_tests", exact_tests, "hits", #hits;
print "kill_counts";
for p in Sort([k : k in Keys(kill_counts)]) do
    print p, kill_counts[p];
end for;
for H in hits do
    print "H", H;
end for;

quit;
