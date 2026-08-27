//////////////////////////////////////////////////////////////////////
//  Contact-5 family with an explicit order-4 lift, and search for
//  order 40 specializations.
//
//  Start from
//      h = 1 + t*x + ((t^2 - 1)/2)*x^2,
//      f = h^2 - ((t + 1)^4/4)*x^5.
//
//  Then x=1 is a rational Weierstrass point and the contact identity
//  gives a rational 5-torsion class.  For nonsingular specializations,
//  the divisor
//
//      H = [ x^2 + 2*x/(t+1), (t+2)*x + 1 ]
//
//  halves the Weierstrass class [x-1,0], so H has order 4.  This gives
//  a rational point of order 20.  The script searches for t where H is
//  itself divisible by 2, giving order 40.
//
//  Typical run:
//      magma -b height:=35 code/contact5_order40_family_search.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 35;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned max_hits then
    max_hits := 20;
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

function FrobeniusIrreducibleCertificate(f)
    C := HyperellipticCurve(f);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71] do
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
order20 := 0;
order40 := 0;
hits := [];

print "Contact-5 order-20 family; searching for order 40";
print "height", height, "parameters", #params;

for t in params do
    if t eq -1 then
        continue;
    end if;

    f, h, b := FamilyPolynomial(t);
    checked +:= 1;
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        continue;
    end if;
    smooth +:= 1;

    fI, L := IntegralModel(f);
    C := HyperellipticCurve(fI);
    J := Jacobian(C);

    // Under yI = L*y, the contact point (0,1) becomes (0,L).
    D5 := J![x, Q!L];
    assert 5*D5 eq J!0;
    assert D5 ne J!0;

    H := J![x^2 + (2/(t+1))*x, L*((t+2)*x + 1)];
    assert 2*H eq J![x - 1, Q!0];
    assert Order(H) eq 4;
    order20 +:= 1;

    ok, Q8 := IsDivisibleBy(H, 2);
    if ok then
        G, phi := TorsionSubgroup(J);
        invs := Invariants(G);
        simple, pcert, Lp := FrobeniusIrreducibleCertificate(fI);
        order40 +:= 1;
        Append(~hits, <t,b,invs,simple,pcert,Lp,fI,Q8>);
        print "HIT", "t", t, "b", b, "torsion", invs,
              "Q8_order", Order(Q8), "simple_cert", simple, "pcert", pcert;
        print "  f =", fI;
        print "  Q8 =", Q8;
        if #hits ge max_hits then
            break;
        end if;
    end if;
end for;

print "DONE height", height;
print "checked", checked, "smooth", smooth, "order20", order20,
      "order40_hits", order40;

quit;
