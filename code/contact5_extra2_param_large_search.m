//////////////////////////////////////////////////////////////////////
//  Parametrized search for larger torsion on the contact-5/order-20
//  family with an extra rational 2-torsion class.
//
//  With u = t+1 and y = u*x, the residual quartic f/(x-1) is
//
//      u*(y^4 + 4*y^3 + 8*y^2 + 8*y + 4) - 4*y*(y+1)^2.
//
//  Hence the 1+3 factor locus is parametrized by
//
//      t = -(z^4 + 4*z + 4)/(z^4 + 4*z^3 + 8*z^2 + 8*z + 4).
//
//  The 2+2 factor locus is parametrized by
//
//      t = -(r^6 - 2*r^5 + 2*r^4 - 4*r^3 + 4*r^2 - 8*r + 8)
//           /((r^2 - 2)^2*(r^2 - 2*r + 2)).
//
//  This script enumerates rational parameters in these two loci, first
//  using precomputed finite-field residue data to bound the possible
//  rational torsion order, and exact-computing torsion only when the
//  modular gcd remains larger than the selected threshold.  The default
//  threshold is 40, the forced torsion order on the extra-2 locus.
//
//  Typical run:
//      magma -b height:=3000 code/contact5_extra2_param_large_search.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 1000;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned max_exact then
    max_exact := 1000;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned threshold then
    threshold := 40;
elif Type(threshold) eq MonStgElt then
    threshold := StringToInteger(threshold);
end if;
if not assigned progress_interval then
    progress_interval := 1000000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

primes := [3,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113];

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

function FamilyPolynomial(t)
    b := (t^2 - 1)/2;
    h := 1 + t*x + b*x^2;
    f := h^2 - ((t + 1)^4/4)*x^5;
    return f;
end function;

function IntegralModel(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function TorsionOrder(invs)
    if #invs eq 0 then
        return 1;
    end if;
    return &*invs;
end function;

function FactorTypeString(fac)
    degs := Sort([ Degree(ff[1]) : ff in fac ]);
    return Join([ IntegerToString(d) : d in degs ], "+");
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

function FamilyPolynomialFinite(t, F)
    PF<xx> := PolynomialRing(F);
    inv2 := (F!2)^-1;
    inv4 := (F!4)^-1;
    b := (t^2 - 1)*inv2;
    h := 1 + t*xx + b*xx^2;
    f := h^2 - ((t + 1)^4)*inv4*xx^5;
    return f;
end function;

function ResidueJacobianOrder(f)
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false, 0;
    end if;
    C := HyperellipticCurve(f);
    Lp := LPolynomial(C);
    return true, Z!Evaluate(Lp, 1);
end function;

function PrecomputeResidues(label, p)
    F := GF(p);
    data := AssociativeArray();
    valid := 0;
    for a in [0..p-1] do
        z := F!a;
        if label eq "linear" then
            ok, t := TLinearFinite(z, F);
        else
            ok, t := TQuadraticQuadraticFinite(z, F);
        end if;
        if not ok or t eq -F!1 then
            continue;
        end if;
        f := FamilyPolynomialFinite(t, F);
        good, N := ResidueJacobianOrder(f);
        if not good then
            continue;
        end if;
        data[a] := N;
        valid +:= 1;
    end for;
    return data, valid;
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

residue_data := AssociativeArray();
for label in ["linear", "qq"] do
    residue_data[label] := AssociativeArray();
    print "Precomputing residues", label;
    for p in primes do
        data, valid := PrecomputeResidues(label, p);
        residue_data[label][p] := data;
        print " ", p, "valid", valid;
    end for;
end for;

params := RationalParametersOfHeight(height);
seen_t := {};

checked := AssociativeArray();
survivors := AssociativeArray();
exact_tests := 0;
large_hits := [];
bound_counts := AssociativeArray();

for label in ["linear", "qq"] do
    checked[label] := 0;
    survivors[label] := 0;
end for;

print "Parametrized extra-2 larger torsion search";
print "height", height, "parameters_per_family", #params, "threshold", threshold,
      "max_exact", max_exact;

for label in ["linear", "qq"] do
    for z in params do
        checked[label] +:= 1;
        total_checked := checked["linear"] + checked["qq"];
        if progress_interval gt 0 and total_checked mod progress_interval eq 0 then
            print "progress", total_checked, "linear", checked["linear"],
                  "qq", checked["qq"], "survivors_linear", survivors["linear"],
                  "survivors_qq", survivors["qq"], "exact", exact_tests,
                  "large", #large_hits;
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

        g := 0;
        good_primes := 0;
        for p in primes do
            okres, residue := ResidueOfRational(z, p);
            if not okres then
                continue;
            end if;
            data := residue_data[label][p];
            if not IsDefined(data, residue) then
                continue;
            end if;
            N := data[residue];
            if g eq 0 then
                g := N;
            else
                g := GCD(g, N);
            end if;
            good_primes +:= 1;
            if g eq threshold then
                break;
            end if;
        end for;

        if good_primes eq 0 or g le threshold then
            continue;
        end if;
        survivors[label] +:= 1;

        bkey := Sprint(g);
        if IsDefined(bound_counts, bkey) then
            bound_counts[bkey] +:= 1;
        else
            bound_counts[bkey] := 1;
        end if;

        if exact_tests ge max_exact then
            continue;
        end if;

        f := FamilyPolynomial(t);
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
        print "EXACT", "label", label, "z", z, "t", t, "bound", g,
              "torsion", invs, "order", ord, "factor_type", ftype;
        if ord gt threshold then
            Append(~large_hits, <label,z,t,g,invs,ftype,fI,facq>);
            print "LARGE", "label", label, "z", z, "t", t, "bound", g,
                  "torsion", invs, "order", ord, "factor_type", ftype;
            print "  f =", fI;
            print "  quartic_factorization =", facq;
        end if;
    end for;
end for;

print "DONE height", height;
print "checked_linear", checked["linear"], "checked_qq", checked["qq"],
      "unique_t", #seen_t;
print "survivors_linear", survivors["linear"], "survivors_qq", survivors["qq"],
      "exact_tests", exact_tests, "large_hits", #large_hits;
print "Surviving modular bounds";
for key in Sort([ k : k in Keys(bound_counts) ]) do
    print " ", key, bound_counts[key];
end for;

quit;
