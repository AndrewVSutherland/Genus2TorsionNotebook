//////////////////////////////////////////////////////////////////////
//  Direct contact-6 [1,2,2] curve search for simple [6,6].
//
//  This searches the original contact-6 family
//
//      h6 = 1 + a*x + b*x^2 + x^3,
//      f  = h6^2 - (x-1)^6
//
//  but keeps only factor type [1,2,2], i.e. exactly the built-in rational
//  root x=0 and two irreducible quadratic factors.  Finite-field [6,6]
//  filters are applied before exact TorsionSubgroup computations.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned height then
    height := 20;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned prime_bound then
    prime_bound := 31;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;
if not assigned max_exact then
    max_exact := 5000;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;
if not assigned progress_interval then
    progress_interval := 100000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;
if not assigned simple_only then
    simple_only := true;
elif Type(simple_only) eq MonStgElt then
    simple_only := simple_only in {"true", "True", "1", "yes"};
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
            q := Q!num/den;
            key := Sprint(q);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, q);
            end if;
        end for;
    end for;
    return vals;
end function;

function Contact6Polynomial(a, b)
    h := 1 + a*x + b*x^2 + x^3;
    f := h^2 - (x-1)^6;
    return f, h;
end function;

function GoodPolynomial(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function FactorDegrees(f)
    return Sort([Degree(fe[1]) : fe in Factorization(f)]);
end function;

function HasFactorType122(f)
    return FactorDegrees(f) eq [1,2,2];
end function;

function Has66(invs)
    return #[n : n in invs | (Z!n) mod 6 eq 0] ge 2;
end function;

function SimpleCertificate(f)
    C := HyperellipticCurve(f);
    for p in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97] do
        try
            fp := ChangeRing(f, GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
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

function CurvePolynomialFinite(F, a, b)
    PF<X> := PolynomialRing(F);
    h := 1 + a*X + b*X^2 + X^3;
    f := h^2 - (X-1)^6;
    return f, h;
end function;

function FiniteData(F, a, b)
    f, h := CurvePolynomialFinite(F, a, b);
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false, [];
    end if;
    if Evaluate(h, F!1) eq 0 then
        return false, [];
    end if;
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    A, phi := AbelianGroup(J);
    return true, Invariants(A);
end function;

function ResidueKey(a, b, p)
    F := GF(p);
    try
        aa := F!a;
        bb := F!b;
    catch e
        return false, 0;
    end try;
    ok, invs := FiniteData(F, aa, bb);
    if not ok then
        return false, 0;
    end if;
    return true, Z!aa + p*(Z!bb);
end function;

function AllowedResidues(p)
    F := GF(p);
    allowed := {Z | };
    good := 0;
    target := 0;
    for a in F do
        for b in F do
            ok, invs := FiniteData(F, a, b);
            if not ok then
                continue;
            end if;
            good +:= 1;
            if Has66(invs) then
                Include(~allowed, Z!a + p*(Z!b));
                target +:= 1;
            end if;
        end for;
    end for;
    return allowed, good, target;
end function;

function PassesResidues(a, b, residue_data)
    for data in residue_data do
        p := data[1];
        allowed := data[2];
        ok, key := ResidueKey(a, b, p);
        if ok and key notin allowed then
            return false, p;
        end if;
    end for;
    return true, 0;
end function;

function IntegralModel(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function ExactData(a, b)
    f, h := Contact6Polynomial(a, b);
    if not GoodPolynomial(f) or Evaluate(h, Q!1) eq 0 then
        return false, [], f, false, 0, P!0;
    end if;
    if not HasFactorType122(f) then
        return false, [], f, false, 0, P!0;
    end if;
    simple, pcert, Lp := SimpleCertificate(f);
    if simple_only and not simple then
        return false, [], f, false, 0, P!0;
    end if;
    fI, scale := IntegralModel(f);
    C := HyperellipticCurve(fI);
    J := Jacobian(C);
    G, phi := TorsionSubgroup(J);
    return true, Invariants(G), fI, simple, pcert, Lp;
end function;

params := RationalParametersOfHeight(height);
primes := [p : p in PrimesUpTo(prime_bound) | p notin {2,3}];
residue_data := [];

print "Contact-6 [1,2,2] direct curve search for simple [6,6]";
print "height", height, "parameter_count", #params,
      "prime_bound", prime_bound, "max_exact", max_exact,
      "simple_only", simple_only;
print "precomputing finite filters";
for p in primes do
    allowed, good, target := AllowedResidues(p);
    Append(~residue_data, <p, allowed>);
    print " p", p, "good", good, "allowed66", #allowed, "target", target;
end for;

checked := 0;
smooth122 := 0;
survivors := 0;
exact := 0;
hits := [];
kill := AssociativeArray();

for a in params do
    for b in params do
        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "progress", checked, "smooth122", smooth122,
                  "survivors", survivors, "exact", exact, "hits", #hits;
        end if;

        f, h := Contact6Polynomial(a, b);
        if not GoodPolynomial(f) or Evaluate(h, Q!1) eq 0 then
            continue;
        end if;
        if not HasFactorType122(f) then
            continue;
        end if;
        smooth122 +:= 1;

        pass, pbad := PassesResidues(a, b, residue_data);
        if not pass then
            if IsDefined(kill, pbad) then kill[pbad] +:= 1; else kill[pbad] := 1; end if;
            continue;
        end if;
        survivors +:= 1;
        if exact ge max_exact then
            continue;
        end if;
        ok, invs, fI, simple, pcert, Lp := ExactData(a, b);
        if not ok then
            continue;
        end if;
        exact +:= 1;
        if Has66(invs) then
            Append(~hits, <a,b,invs,simple,pcert,fI>);
            print "HIT66", "a", a, "b", b, "invs", invs,
                  "simple", simple, "pcert", pcert;
            print " f", fI;
            if #hits ge max_hits then
                break a;
            end if;
        end if;
    end for;
end for;

print "Done";
print "checked", checked;
print "smooth122", smooth122;
print "survivors", survivors;
print "exact", exact;
print "hits", #hits;
print "kill", kill;
for H in hits do
    print "H", H;
end for;

quit;
