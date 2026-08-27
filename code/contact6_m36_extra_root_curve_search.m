//////////////////////////////////////////////////////////////////////
//  Direct curve search on the contact-6 extra-root family for [6,6].
//
//  This is broader than searching the cubic-contact variables.  We force
//  a second rational Weierstrass point by
//
//      h6(r) = eps*(r - 1)^3,    eps in {+1,-1},
//
//  so
//
//      a = (eps*(r-1)^3 - 1 - r^3 - b*r^2)/r.
//
//  Then we search the two-parameter family (r,b), using finite-field
//  [6,6] filters before exact TorsionSubgroup computations over Q.
//
//  Typical runs:
//      magma -b height:=12 prime_bound:=31 code/contact6_m36_extra_root_curve_search.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned height then
    height := 12;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned prime_bound then
    prime_bound := 31;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;
if not assigned max_exact then
    max_exact := 1000;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;
if not assigned progress_interval then
    progress_interval := 10000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;
if not assigned require_113 then
    require_113 := false;
elif Type(require_113) eq MonStgElt then
    require_113 := require_113 in {"true", "True", "1", "yes"};
end if;
if not assigned simple_only then
    simple_only := false;
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

function ExtraRootA(eps, r, b)
    return (eps*(r-1)^3 - 1 - r^3 - b*r^2)/r;
end function;

function Contact6Polynomial(a, b)
    h := 1 + a*x + b*x^2 + x^3;
    f := h^2 - (x-1)^6;
    return f, h;
end function;

function CurvePolynomial(eps, r, b)
    a := ExtraRootA(eps, r, b);
    f, h := Contact6Polynomial(a, b);
    return a, f, h;
end function;

function IntegralModel(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function GoodPolynomial(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function Has66(invs)
    return #[n : n in invs | (Z!n) mod 6 eq 0] ge 2;
end function;

function FactorDegrees(f)
    return Sort([Degree(fe[1]) : fe in Factorization(f)]);
end function;

function HasFactorType113(f)
    return FactorDegrees(f) eq [1,1,3];
end function;

function TwoRank(invs)
    return #[n : n in invs | (Z!n) mod 2 eq 0];
end function;

function ThreeRank(invs)
    return #[n : n in invs | (Z!n) mod 3 eq 0];
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

function CurvePolynomialFinite(F, eps, r, b)
    PF<X> := PolynomialRing(F);
    a := (eps*(r-1)^3 - 1 - r^3 - b*r^2)/r;
    h := 1 + a*X + b*X^2 + X^3;
    f := h^2 - (X-1)^6;
    return a, f, h;
end function;

function FiniteData(F, eps, r, b)
    a, f, h := CurvePolynomialFinite(F, eps, r, b);
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false, [], a;
    end if;
    if Evaluate(h, F!1) eq 0 or Evaluate(f, r) ne 0 then
        return false, [], a;
    end if;
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    A, phi := AbelianGroup(J);
    return true, Invariants(A), a;
end function;

function ResidueKey(eps, r, b, p)
    F := GF(p);
    try
        rr := F!r;
        bb := F!b;
    catch e
        return false, 0;
    end try;
    if rr eq 0 or rr eq 1 then
        return false, 0;
    end if;
    ee := eps eq 1 select F!1 else -F!1;
    epsbit := eps eq 1 select 0 else 1;
    ok, invs, a := FiniteData(F, ee, rr, bb);
    if not ok then
        return false, 0;
    end if;
    return true, epsbit + 2*(Z!rr + p*(Z!bb));
end function;

function AllowedResidues(p)
    F := GF(p);
    allowed := {Z | };
    good := 0;
    target := 0;
    rank_counts := AssociativeArray();
    for eps in [F!1, -F!1] do
        epsbit := eps eq F!1 select 0 else 1;
        for r in F do
            if r eq 0 or r eq 1 then
                continue;
            end if;
            for b in F do
                ok, invs, a := FiniteData(F, eps, r, b);
                if not ok then
                    continue;
                end if;
                good +:= 1;
                keyrank := Sprintf("%o,%o", TwoRank(invs), ThreeRank(invs));
                if IsDefined(rank_counts, keyrank) then
                    rank_counts[keyrank] +:= 1;
                else
                    rank_counts[keyrank] := 1;
                end if;
                if Has66(invs) then
                    Include(~allowed, epsbit + 2*(Z!r + p*(Z!b)));
                    target +:= 1;
                end if;
            end for;
        end for;
    end for;
    return allowed, good, target, rank_counts;
end function;

function PassesResidues(eps, r, b, residue_data)
    for data in residue_data do
        p := data[1];
        allowed := data[2];
        ok, key := ResidueKey(eps, r, b, p);
        if ok and key notin allowed then
            return false, p;
        end if;
    end for;
    return true, 0;
end function;

function ExactData(eps, r, b)
    a, f, h := CurvePolynomial(eps, r, b);
    if r eq 0 or r eq 1 or not GoodPolynomial(f) then
        return false, [], a, f, false, 0, P!0;
    end if;
    if Evaluate(h, Q!1) eq 0 or Evaluate(f, r) ne 0 then
        return false, [], a, f, false, 0, P!0;
    end if;
    fI, scale := IntegralModel(f);
    if not GoodPolynomial(fI) then
        return false, [], a, fI, false, 0, P!0;
    end if;
    C := HyperellipticCurve(fI);
    J := Jacobian(C);
    G, phi := TorsionSubgroup(J);
    return true, Invariants(G), a, fI, false, 0, P!0;
end function;

params := RationalParametersOfHeight(height);
primes := [p : p in PrimesUpTo(prime_bound) | p notin {2,3}];
residue_data := [];

print "Contact-6 extra-root direct curve search for [6,6]";
print "height", height, "parameter_count", #params,
      "prime_bound", prime_bound, "max_exact", max_exact,
      "require_113", require_113, "simple_only", simple_only;
print "precomputing finite filters";
for p in primes do
    allowed, good, target, rank_counts := AllowedResidues(p);
    Append(~residue_data, <p, allowed>);
    print " p", p, "good", good, "allowed66", #allowed, "target", target;
end for;

checked := 0;
smooth := 0;
survivors := 0;
exact := 0;
hits := [];
kill := AssociativeArray();

for eps in [1, -1] do
    for r in params do
        if r eq 0 or r eq 1 then
            continue;
        end if;
        for b in params do
            checked +:= 1;
            if progress_interval gt 0 and checked mod progress_interval eq 0 then
                print "progress", checked, "smooth", smooth,
                      "survivors", survivors, "exact", exact, "hits", #hits;
            end if;

            a, f, h := CurvePolynomial(Q!eps, r, b);
            if not GoodPolynomial(f) or Evaluate(h, Q!1) eq 0 then
                continue;
            end if;
            if require_113 and not HasFactorType113(f) then
                continue;
            end if;
            smooth +:= 1;

            pass, pbad := PassesResidues(eps, r, b, residue_data);
            if not pass then
                if IsDefined(kill, pbad) then kill[pbad] +:= 1; else kill[pbad] := 1; end if;
                continue;
            end if;
            survivors +:= 1;
            pre_simple := false;
            pre_pcert := 0;
            if simple_only then
                fIpre, scale_pre := IntegralModel(f);
                pre_simple, pre_pcert, pre_Lp := SimpleCertificate(fIpre);
                if not pre_simple then
                    continue;
                end if;
            end if;
            if exact ge max_exact then
                continue;
            end if;

            ok, invs, a0, fI, simple, pcert, Lp := ExactData(Q!eps, r, b);
            if not ok then
                continue;
            end if;
            exact +:= 1;
            if Has66(invs) then
                if simple_only then
                    simple := pre_simple;
                    pcert := pre_pcert;
                else
                    simple, pcert, Lp := SimpleCertificate(fI);
                end if;
                if simple_only and not simple then
                    continue;
                end if;
                Append(~hits, <eps,r,a0,b,invs,simple,pcert,fI>);
                print "HIT66", "eps", eps, "r", r, "a", a0, "b", b,
                      "factor_degrees", FactorDegrees(f),
                      "invs", invs, "simple", simple, "pcert", pcert;
                print " f", fI;
                if #hits ge max_hits then
                    break eps;
                end if;
            end if;
        end for;
    end for;
end for;

print "Done";
print "checked", checked;
print "smooth", smooth;
print "survivors", survivors;
print "exact", exact;
print "hits", #hits;
print "kill", kill;
for H in hits do
    print "H", H;
end for;

quit;
