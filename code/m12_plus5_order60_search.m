//////////////////////////////////////////////////////////////////////
//  M(12) one-parameter family plus possible rational 5-torsion.
//
//  The family a=(1-r)/4 has a rational point D of order 12 on the
//  Jacobian, and an independent rational 2-torsion point.  A rational
//  5-torsion class would therefore give an element of order 60.
//
//  This script first applies the necessary good-reduction condition
//      5 | #J(F_p)
//  for all good p != 5 up to prime_bound.  Bad/boundary residues are
//  allowed through.  Survivors are exact-checked with TorsionSubgroup.
//
//  Typical runs:
//      magma -b height:=500 prime_bound:=251 max_exact:=200 \
//          code/m12_plus5_order60_search.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned height then
    height := 500;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned prime_bound then
    prime_bound := 251;
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
    progress_interval := 100000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Q := Rationals();
Z := Integers();
P<X> := PolynomialRing(Q);

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            r := Q!num / Q!den;
            key := Sprint(r);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, r);
            end if;
        end for;
    end for;
    return vals;
end function;

function OddQuinticForR(r)
    a := (1-r)/4;
    T := a*X^2 - X + r;
    h := (X-r)*(T+1);
    W := h^2 + 4*a*X^2*T*(T+1);

    // Move the root x=2 of T+1 to infinity:
    // x = 2 + 1/X, Y_old = Y_new/X^3.
    f5 := P!0;
    for i in [0..Degree(W)] do
        for j in [0..i] do
            f5 +:= Coefficient(W, i)*Binomial(i,j)*2^(i-j)*X^(6-j);
        end for;
    end for;
    return f5;
end function;

function OddQuinticForRFinite(r, F)
    PF<XX> := PolynomialRing(F);
    a := (1-r)/(F!4);
    T := a*XX^2 - XX + r;
    h := (XX-r)*(T+1);
    W := h^2 + 4*a*XX^2*T*(T+1);
    f5 := PF!0;
    for i in [0..Degree(W)] do
        for j in [0..i] do
            f5 +:= Coefficient(W, i)*Binomial(i,j)*(F!2)^(i-j)*XX^(6-j);
        end for;
    end for;
    return f5;
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function GoodHyperellipticPolynomial(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function ResidueStatus(p)
    F := GF(p);
    allowed := {};
    bad := {};
    good := 0;
    killed := 0;
    for a in [0..p-1] do
        r := F!a;
        f := OddQuinticForRFinite(r, F);
        if not GoodHyperellipticPolynomial(f) then
            Include(~bad, a);
            continue;
        end if;
        good +:= 1;
        C := HyperellipticCurve(f);
        N := Z!Evaluate(LPolynomial(C), 1);
        if (N mod 5) eq 0 then
            Include(~allowed, a);
        else
            killed +:= 1;
        end if;
    end for;
    return allowed, bad, good, killed;
end function;

function ResidueOfRational(r, p)
    num := Numerator(r);
    den := Denominator(r);
    if (den mod p) eq 0 then
        return false, 0;
    end if;
    F := GF(p);
    return true, Z!(F!num / F!den);
end function;

function PassesFiveResidues(r, allowed, badres, primes)
    bad_primes := [];
    good_checked := 0;
    for p in primes do
        ok, residue := ResidueOfRational(r, p);
        if not ok then
            Append(~bad_primes, p);
            continue;
        end if;
        if residue in badres[p] then
            Append(~bad_primes, p);
            continue;
        end if;
        good_checked +:= 1;
        if residue notin allowed[p] then
            return false, p, good_checked, bad_primes;
        end if;
    end for;
    return true, 0, good_checked, bad_primes;
end function;

function TorsionOrderFromInvariants(invs)
    if #invs eq 0 then
        return 1;
    end if;
    return &*invs;
end function;

function TorsionExponentFromInvariants(invs)
    e := 1;
    for n in invs do
        e := LCM(e, n);
    end for;
    return e;
end function;

function PointOfExactOrder(G, phi, target)
    for g in G do
        og := Order(g);
        if og mod target eq 0 then
            h := (og div target)*g;
            if Order(h) eq target then
                return true, phi(h), h;
            end if;
        end if;
    end for;
    return false, Codomain(phi)!0, G!0;
end function;

primes := [ p : p in PrimesUpTo(prime_bound) | p notin {2,5} ];
allowed := AssociativeArray();
badres := AssociativeArray();

print "M(12) one-parameter plus 5/order-60 search";
print "height", height, "prime_bound", prime_bound, "primes", primes,
      "max_exact", max_exact;
print "Precomputing 5-divisibility residue filters";
for p in primes do
    good_allowed, bad, good, killed := ResidueStatus(p);
    allowed[p] := good_allowed;
    badres[p] := bad;
    print " ", p, "good", good, "allowed5", #good_allowed,
          "bad", #bad, "killed", killed;
end for;

params := RationalParametersOfHeight(height);
checked := 0;
smooth := 0;
residue_survivors := 0;
exact_tests := 0;
hits := [];
printed := 0;
kill_counts := AssociativeArray();
bad_signature_counts := AssociativeArray();

for r in params do
    checked +:= 1;
    if r in {Q!0, Q!1, Q!2} then
        continue;
    end if;
    if progress_interval gt 0 and checked mod progress_interval eq 0 then
        print "progress", checked, "smooth", smooth,
              "residue_survivors", residue_survivors,
              "exact", exact_tests, "hits", #hits;
    end if;

    pass, killp, good_checked, bad_primes := PassesFiveResidues(r, allowed, badres, primes);
    if not pass then
        if IsDefined(kill_counts, killp) then
            kill_counts[killp] +:= 1;
        else
            kill_counts[killp] := 1;
        end if;
        continue;
    end if;

    f5 := OddQuinticForR(r);
    if not GoodHyperellipticPolynomial(f5) then
        continue;
    end if;
    smooth +:= 1;
    residue_survivors +:= 1;
    sig := Join([ IntegerToString(p) : p in bad_primes ], ",");
    if IsDefined(bad_signature_counts, sig) then
        bad_signature_counts[sig] +:= 1;
    else
        bad_signature_counts[sig] := 1;
    end if;

    if printed lt max_print then
        print "SURVIVOR", "r", r, "good_checked", good_checked,
              "bad_primes", bad_primes;
        printed +:= 1;
    end if;

    if exact_tests ge max_exact then
        continue;
    end if;

    fI, L := IntegralModelPolynomial(f5);
    if not GoodHyperellipticPolynomial(fI) then
        continue;
    end if;
    C := HyperellipticCurve(fI);
    J := Jacobian(C);
    G, phi := TorsionSubgroup(J);
    invs := Invariants(G);
    ord := TorsionOrderFromInvariants(invs);
    exp := TorsionExponentFromInvariants(invs);
    exact_tests +:= 1;

    // Sanity check for the known order-12 point, scaled to the integral model.
    Yp := (Evaluate(((X-r)*(((1-r)/4)*X^2 - X + r + 1)), 0))*(-Q!1/2)^3;
    D12 := J![X + Q!1/2, L*Yp];
    D12_order := Order(D12);

    print "EXACT", "r", r, "torsion", invs, "order", ord,
          "exponent", exp, "D12_order", D12_order,
          "bad_primes", bad_primes;

    if exp mod 60 eq 0 then
        ok60, P60, g60 := PointOfExactOrder(G, phi, 60);
        Append(~hits, <r, invs, ord, exp, fI, P60>);
        print "ORDER60_HIT", "r", r, "torsion", invs,
              "order", ord, "exponent", exp;
        print "  f =", fI;
        if ok60 then
            print "  P60 =", P60;
            print "  abstract_g60 =", g60;
        else
            print "  WARNING exponent divisible by 60 but explicit point extraction failed";
        end if;
    end if;
end for;

print "DONE M12 plus 5/order 60";
print "checked", checked, "smooth", smooth,
      "residue_survivors", residue_survivors,
      "exact_tests", exact_tests, "hits", #hits;
print "kill_counts";
for p in Sort([k : k in Keys(kill_counts)]) do
    print p, kill_counts[p];
end for;
print "bad_signature_counts";
for key in Sort([ k : k in Keys(bad_signature_counts) ]) do
    print " ", key, bad_signature_counts[key];
end for;
for H in hits do
    print H;
end for;
quit;
