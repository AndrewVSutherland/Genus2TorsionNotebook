//////////////////////////////////////////////////////////////////////
//  Pure modular sieve for torsion larger than a chosen threshold in the
//  parametrized extra-2 loci of the contact-5/order-20 family.
//
//  This is a faster companion to contact5_extra2_param_large_search.m:
//  it only computes gcds of #J(F_p) over many small good primes and
//  reports parameters whose modular bound remains above threshold.
//
//  Typical run:
//      magma -b height:=1000 threshold:=80 prime_bound:=251 \
//          code/contact5_extra2_param_modsieve.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 1000;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned threshold then
    threshold := 80;
elif Type(threshold) eq MonStgElt then
    threshold := StringToInteger(threshold);
end if;
if not assigned prime_bound then
    prime_bound := 251;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
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

function FamilyPolynomialFinite(t, F)
    PF<x> := PolynomialRing(F);
    inv2 := (F!2)^-1;
    inv4 := (F!4)^-1;
    b := (t^2 - 1)*inv2;
    h := 1 + t*x + b*x^2;
    return h^2 - ((t + 1)^4)*inv4*x^5;
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
        if good then
            data[a] := N;
            valid +:= 1;
        end if;
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

primes := [ p : p in PrimesUpTo(prime_bound) | p notin {2,5} ];
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
bound_counts := AssociativeArray();
printed := 0;

for label in ["linear", "qq"] do
    checked[label] := 0;
    survivors[label] := 0;
end for;

print "Contact-5 extra-2 modular sieve";
print "height", height, "parameters_per_family", #params,
      "threshold", threshold, "prime_bound", prime_bound, "primes", #primes;

for label in ["linear", "qq"] do
    for z in params do
        checked[label] +:= 1;
        total_checked := checked["linear"] + checked["qq"];
        if progress_interval gt 0 and total_checked mod progress_interval eq 0 then
            print "progress", total_checked, "survivors_linear", survivors["linear"],
                  "survivors_qq", survivors["qq"], "printed", printed;
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
            if g le threshold then
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

        if printed lt max_print then
            print "SURVIVOR", "label", label, "z", z, "t", t,
                  "bound", g, "good_primes", good_primes;
            printed +:= 1;
        end if;
    end for;
end for;

print "DONE height", height;
print "checked_linear", checked["linear"], "checked_qq", checked["qq"],
      "unique_t", #seen_t;
print "survivors_linear", survivors["linear"], "survivors_qq", survivors["qq"];
print "Surviving modular bounds";
for key in Sort([ k : k in Keys(bound_counts) ]) do
    print " ", key, bound_counts[key];
end for;

quit;
