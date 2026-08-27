//////////////////////////////////////////////////////////////////////
//  Contact-7 family and first extra-even-torsion search.
//
//  On an odd genus-2 model y^2 = f(x), the function x*y - h(x)
//  has pole order 7 at infinity when deg(h) <= 3.  If
//
//      h(x)^2 - x^2 f(x) = -(x - 1)^7,
//
//  then div(x*y - h) = 7*P - 7*infinity for
//  P = (1, h(1)), provided h(1) != 0 and the curve is smooth.
//
//  The divisibility by x^2 is forced by taking
//
//      h = 1 - (7/2)*x + a*x^2 + b*x^3.
//
//  This script searches the resulting two-parameter family for
//  extra rational even torsion by looking at factorization of f.
//
//  Typical runs:
//      magma -b mode:="verify" code/contact7_family_search.m
//      magma -b height:=20 mode:="integer" code/contact7_family_search.m
//      magma -b height:=8 mode:="rational" code/contact7_family_search.m
//////////////////////////////////////////////////////////////////////

if not assigned mode then
    mode := "verify";
end if;
if not assigned height then
    height := 12;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned max_exact then
    max_exact := 5000;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned progress_interval then
    progress_interval := 10000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Q := Rationals();
P<x> := PolynomialRing(Q);

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

function IntegralModel(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function Contact7Polynomial(a, b)
    h := 1 - (Q!7/2)*x + a*x^2 + b*x^3;
    num := h^2 + (x - 1)^7;
    if Coefficient(num, 0) ne 0 or Coefficient(num, 1) ne 0 then
        return false, P!0, h;
    end if;
    f := ExactQuotient(num, x^2);
    return true, f, h;
end function;

function GoodHyperellipticPolynomial(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
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

function HasEvenTorsionFromFactorization(fac)
    return #fac ge 2;
end function;

function SimpleCertificate(f)
    C := HyperellipticCurve(f);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97] do
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

function VerifyContact7(a, b)
    ok, f, h := Contact7Polynomial(a, b);
    if not ok or not GoodHyperellipticPolynomial(f) then
        return false, 0, f, h, "bad polynomial";
    end if;
    yP := Evaluate(h, Q!1);
    if yP eq 0 or Evaluate(f, Q!1) ne yP^2 then
        return false, 0, f, h, "bad marked point";
    end if;
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    D := J![x - 1, yP];
    ord := Order(D);
    return ord eq 7, ord, f, h, "";
end function;

procedure VerifySamples()
    samples := [ <Q!0,Q!0>, <Q!1,Q!0>, <Q!0,Q!1>, <Q!2,Q!-1>, <Q!-3,Q!2> ];
    print "CONTACT7 sample verification";
    for ab in samples do
        ok, ord, f, h, msg := VerifyContact7(ab[1], ab[2]);
        print "a", ab[1], "b", ab[2], "ok", ok, "order", ord, "msg", msg;
        if ok then
            print "  h =", h;
            print "  f =", f;
            print "  factor_type", FactorTypeString(Factorization(f));
        end if;
    end for;
end procedure;

procedure SearchContact7(params, label)
    checked := 0;
    smooth := 0;
    verified7 := 0;
    reducible := 0;
    exact_tests := 0;
    hits := [];
    type_counts := AssociativeArray();
    torsion_counts := AssociativeArray();

    print "CONTACT7 search", label;
    print "height", height, "parameters", #params, "max_exact", max_exact;

    for a in params do
        for b in params do
            checked +:= 1;
            if progress_interval gt 0 and checked mod progress_interval eq 0 then
                print "progress", "checked", checked, "smooth", smooth,
                      "verified7", verified7, "reducible", reducible,
                      "exact_tests", exact_tests, "hits", #hits;
            end if;

            ok, f, h := Contact7Polynomial(a, b);
            if not ok or not GoodHyperellipticPolynomial(f) then
                continue;
            end if;
            smooth +:= 1;

            yP := Evaluate(h, Q!1);
            if yP eq 0 then
                continue;
            end if;

            C := HyperellipticCurve(f);
            J := Jacobian(C);
            D := J![x - 1, yP];
            if Order(D) ne 7 then
                continue;
            end if;
            verified7 +:= 1;

            fac := Factorization(f);
            ftype := FactorTypeString(fac);
            if IsDefined(type_counts, ftype) then
                type_counts[ftype] +:= 1;
            else
                type_counts[ftype] := 1;
            end if;

            if not HasEvenTorsionFromFactorization(fac) then
                continue;
            end if;
            reducible +:= 1;
            if exact_tests ge max_exact then
                continue;
            end if;

            fI, L := IntegralModel(f);
            CI := HyperellipticCurve(fI);
            JI := Jacobian(CI);
            G, phi := TorsionSubgroup(JI);
            invs := Invariants(G);
            exact_tests +:= 1;

            invkey := Sprint(invs);
            if IsDefined(torsion_counts, invkey) then
                torsion_counts[invkey] +:= 1;
            else
                torsion_counts[invkey] := 1;
            end if;

            ord := TorsionOrder(invs);
            if ord gt 7 then
                simple, pcert, Lp := SimpleCertificate(fI);
                Append(~hits, <a,b,invs,ftype,simple,pcert,fI,fac>);
                print "HIT", "a", a, "b", b, "torsion", invs,
                      "order", ord, "factor_type", ftype,
                      "simple", simple, "pcert", pcert;
                print "  f =", fI;
                print "  factorization =", fac;
            end if;
        end for;
    end for;

    print "DONE contact7", label;
    print "checked", checked, "smooth", smooth, "verified7", verified7,
          "reducible", reducible, "exact_tests", exact_tests,
          "hits", #hits;
    print "Factor type counts";
    for key in Sort([ k : k in Keys(type_counts) ]) do
        print " ", key, type_counts[key];
    end for;
    print "Torsion counts among exact reducible cases";
    for key in Sort([ k : k in Keys(torsion_counts) ]) do
        print " ", key, torsion_counts[key];
    end for;
end procedure;

if mode eq "verify" then
    VerifySamples();
elif mode eq "integer" then
    params := [ Q!n : n in [-height..height] ];
    SearchContact7(params, "integer");
elif mode eq "rational" then
    params := RationalParametersOfHeight(height);
    SearchContact7(params, "rational");
else
    print "unknown mode", mode;
end if;

quit;
