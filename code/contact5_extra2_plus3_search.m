//////////////////////////////////////////////////////////////////////
//  Try to add rational 3-torsion to the contact-5/order-20 family on
//  the independent-extra-2 loci.
//
//  The two parametrized extra-2 loci are:
//
//      linear:
//          t = -(z^4 + 4*z + 4)/(z^4 + 4*z^3 + 8*z^2 + 8*z + 4)
//
//      qq:
//          t = -(r^6 - 2*r^5 + 2*r^4 - 4*r^3 + 4*r^2 - 8*r + 8)
//               /((r^2 - 2)^2*(r^2 - 2*r + 2)).
//
//  Each smooth specialization has torsion containing Z/2 x Z/20.  If a
//  rational 3-torsion point is also present, then for every good prime
//  p != 3 the reduction has 3 | #J(F_p).  This script precomputes the
//  allowed residues modulo many primes and scans rational parameters.
//
//  Typical runs:
//      magma -b height:=1000 prime_bound:=73 code/contact5_extra2_plus3_search.m
//      magma -b height:=1000 prime_bound:=251 code/contact5_extra2_plus3_search.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 1000;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned prime_bound then
    prime_bound := 73;
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
P<x> := PolynomialRing(Q);

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            z := Q!num/den;
            key := Sprint(z);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, z);
            end if;
        end for;
    end for;
    return vals;
end function;

function TLinear(z)
    den := z^4 + 4*z^3 + 8*z^2 + 8*z + 4;
    if den eq 0 then
        return false, Q!0;
    end if;
    return true, -(z^4 + 4*z + 4)/den;
end function;

function TQuadraticQuadratic(r)
    den := (r^2 - 2)^2*(r^2 - 2*r + 2);
    if den eq 0 then
        return false, Q!0;
    end if;
    num := r^6 - 2*r^5 + 2*r^4 - 4*r^3 + 4*r^2 - 8*r + 8;
    return true, -num/den;
end function;

function TLinearFinite(z, F)
    den := z^4 + 4*z^3 + 8*z^2 + 8*z + 4;
    if den eq 0 then
        return false, F!0;
    end if;
    return true, -(z^4 + 4*z + 4)/den;
end function;

function TQuadraticQuadraticFinite(r, F)
    den := (r^2 - 2)^2*(r^2 - 2*r + 2);
    if den eq 0 then
        return false, F!0;
    end if;
    num := r^6 - 2*r^5 + 2*r^4 - 4*r^3 + 4*r^2 - 8*r + 8;
    return true, -num/den;
end function;

function FamilyPolynomial(t)
    b := (t^2 - 1)/2;
    h := 1 + t*x + b*x^2;
    f := h^2 - ((t + 1)^4/4)*x^5;
    return f, b;
end function;

function FamilyPolynomialFinite(t, F)
    PF<xx> := PolynomialRing(F);
    inv2 := (F!2)^-1;
    inv4 := (F!4)^-1;
    b := (t^2 - 1)*inv2;
    h := 1 + t*xx + b*xx^2;
    return h^2 - ((t + 1)^4)*inv4*xx^5;
end function;

function IntegralModel(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function FactorTypeString(fac)
    degs := Sort([ Degree(ff[1]) : ff in fac ]);
    return Join([ IntegerToString(d) : d in degs ], "+");
end function;

function TorsionOrder(invs)
    if #invs eq 0 then
        return 1;
    end if;
    return &*invs;
end function;

function ResidueStatus(label, p)
    F := GF(p);
    good_allowed := {};
    good_killed := 0;
    bad := {};
    for a in [0..p-1] do
        z := F!a;
        if label eq "linear" then
            ok, t := TLinearFinite(z, F);
        else
            ok, t := TQuadraticQuadraticFinite(z, F);
        end if;
        if not ok or t eq -F!1 then
            Include(~bad, a);
            continue;
        end if;
        f := FamilyPolynomialFinite(t, F);
        if Degree(f) ne 5 or Discriminant(f) eq 0 then
            Include(~bad, a);
            continue;
        end if;
        C := HyperellipticCurve(f);
        N := Z!Evaluate(LPolynomial(C), 1);
        if (N mod 3) eq 0 then
            Include(~good_allowed, a);
        else
            good_killed +:= 1;
        end if;
    end for;
    return good_allowed, bad, good_killed;
end function;

function ResidueOfRational(z, p)
    num := Numerator(z);
    den := Denominator(z);
    if (den mod p) eq 0 then
        return false, 0;
    end if;
    F := GF(p);
    return true, Z!(F!num / F!den);
end function;

primes := [ p : p in PrimesUpTo(prime_bound) | p notin {2,3,5} ];

allowed := AssociativeArray();
badres := AssociativeArray();
for label in ["linear", "qq"] do
    allowed[label] := AssociativeArray();
    badres[label] := AssociativeArray();
    print "Precomputing 3-torsion residue filter", label;
    for p in primes do
        good_allowed, bad, good_killed := ResidueStatus(label, p);
        allowed[label][p] := good_allowed;
        badres[label][p] := bad;
        print " ", p, "allowed_good", #good_allowed, "bad", #bad,
              "killed_good", good_killed;
    end for;
end for;

params := RationalParametersOfHeight(height);
seen_t := {};
checked := AssociativeArray();
survivors := AssociativeArray();
bad_signature_counts := AssociativeArray();
exact_tests := 0;
printed := 0;
hits := [];

for label in ["linear", "qq"] do
    checked[label] := 0;
    survivors[label] := 0;
end for;

print "Contact-5 [2,20] plus 3-torsion search";
print "height", height, "parameters_per_family", #params,
      "prime_bound", prime_bound, "primes", #primes,
      "max_exact", max_exact;

for label in ["linear", "qq"] do
    for z in params do
        checked[label] +:= 1;
        total_checked := checked["linear"] + checked["qq"];
        if progress_interval gt 0 and total_checked mod progress_interval eq 0 then
            print "progress", total_checked, "survivors_linear", survivors["linear"],
                  "survivors_qq", survivors["qq"], "exact", exact_tests,
                  "hits", #hits;
        end if;

        if label eq "linear" then
            ok, t := TLinear(z);
        else
            ok, t := TQuadraticQuadratic(z);
        end if;
        if not ok or t eq -Q!1 then
            continue;
        end if;
        tkey := Sprint(t);
        if tkey in seen_t then
            continue;
        end if;
        Include(~seen_t, tkey);

        killed := false;
        good_checked := 0;
        bad_primes := [];
        for p in primes do
            okres, residue := ResidueOfRational(z, p);
            if not okres then
                Append(~bad_primes, p);
                continue;
            end if;
            if residue in badres[label][p] then
                Append(~bad_primes, p);
                continue;
            end if;
            good_checked +:= 1;
            if residue notin allowed[label][p] then
                killed := true;
                break;
            end if;
        end for;
        if killed then
            continue;
        end if;

        survivors[label] +:= 1;
        sig := Join([ IntegerToString(p) : p in bad_primes ], ",");
        if IsDefined(bad_signature_counts, sig) then
            bad_signature_counts[sig] +:= 1;
        else
            bad_signature_counts[sig] := 1;
        end if;

        if printed lt max_print then
            print "SURVIVOR", "label", label, "z", z, "t", t,
                  "good_checked", good_checked, "bad_primes", bad_primes;
            printed +:= 1;
        end if;

        if exact_tests ge max_exact then
            continue;
        end if;

        f, b := FamilyPolynomial(t);
        if Degree(f) ne 5 or Discriminant(f) eq 0 then
            continue;
        end if;
        q := ExactQuotient(f, x - 1);
        facq := Factorization(q);
        ftype := FactorTypeString(facq);

        fI, L := IntegralModel(f);
        C := HyperellipticCurve(fI);
        J := Jacobian(C);
        G, phi := TorsionSubgroup(J);
        invs := Invariants(G);
        exact_tests +:= 1;
        ord := TorsionOrder(invs);

        print "EXACT", "label", label, "z", z, "t", t, "torsion", invs,
              "order", ord, "factor_type", ftype, "bad_primes", bad_primes;

        if (ord mod 3) eq 0 then
            Append(~hits, <label,z,t,b,invs,ftype,fI,facq,bad_primes>);
            print "PLUS3_HIT", "label", label, "z", z, "t", t, "b", b,
                  "torsion", invs, "order", ord, "factor_type", ftype;
            print "  f =", fI;
            print "  residual_factorization =", facq;
        end if;
    end for;
end for;

print "DONE height", height;
print "checked_linear", checked["linear"], "checked_qq", checked["qq"],
      "unique_t", #seen_t;
print "survivors_linear", survivors["linear"], "survivors_qq", survivors["qq"],
      "exact_tests", exact_tests, "plus3_hits", #hits;
print "Bad-prime signatures among survivors";
for key in Sort([ k : k in Keys(bad_signature_counts) ]) do
    print " ", key, bad_signature_counts[key];
end for;

quit;
