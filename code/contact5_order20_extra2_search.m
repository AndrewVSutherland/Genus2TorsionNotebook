//////////////////////////////////////////////////////////////////////
//  Search for an independent rational 2-torsion point in the
//  contact-5/order-20 family.
//
//  In the family
//
//      h = 1 + t*x + ((t^2 - 1)/2)*x^2,
//      f = h^2 - ((t + 1)^4/4)*x^5,
//
//  the factor x-1 gives the rational 2-torsion class halved by the
//  order-4 class H, while the contact identity gives rational 5-torsion.
//  An independent rational 2-torsion point appears exactly when the
//  residual quartic f/(x-1) is reducible over Q.  Such examples should
//  have torsion containing Z/2 x Z/20.
//
//  Typical run:
//      magma -b height:=100 code/contact5_order20_extra2_search.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 100;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned max_hits then
    max_hits := 30;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
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
            t := Q!num/den;
            key := Sprint(t);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, t);
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

function FamilyPolynomial(t)
    b := (t^2 - 1)/2;
    h := 1 + t*x + b*x^2;
    f := h^2 - ((t + 1)^4/4)*x^5;
    return f, h, b;
end function;

function SimpleCertificate(f)
    C := HyperellipticCurve(f);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79] do
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

params := RationalParametersOfHeight(height);
checked := 0;
smooth := 0;
reducible := 0;
torsion_tests := 0;
hits := [];

print "Contact-5/order-20 search for independent rational 2-torsion";
print "height", height, "parameters", #params;

for t in params do
    if t in {Q!-1, Q!-3} then
        continue;
    end if;

    f, h, b := FamilyPolynomial(t);
    checked +:= 1;
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        continue;
    end if;
    smooth +:= 1;

    q := ExactQuotient(f, x - 1);
    fac := Factorization(q);
    if #fac eq 1 then
        continue;
    end if;
    reducible +:= 1;

    fI, L := IntegralModel(f);
    C := HyperellipticCurve(fI);
    J := Jacobian(C);
    torsion_tests +:= 1;
    G, phi := TorsionSubgroup(J);
    invs := Invariants(G);
    simple, pcert, Lp := SimpleCertificate(fI);

    Append(~hits, <t,b,invs,simple,pcert,fI,fac>);
    print "HIT", "t", t, "b", b, "torsion", invs,
          "factor_degrees", [ Degree(ff[1]) : ff in fac ],
          "simple_cert", simple, "pcert", pcert;
    print "  f =", fI;
    print "  quartic_factorization =", fac;

    if #hits ge max_hits then
        break;
    end if;
end for;

print "DONE height", height;
print "checked", checked, "smooth", smooth, "reducible_quartic", reducible,
      "torsion_tests", torsion_tests, "hits", #hits;

quit;
