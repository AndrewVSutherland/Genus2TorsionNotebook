//////////////////////////////////////////////////////////////////////
//  Search the known rational curve on M_1(8,2,2,2) for extra halving.
//
//  NotesAndTodo.tex gives a one-parameter family with torsion containing
//  (2,2,2,8):
//
//      y^2 = (x-1)(x-u)(x-v)(x+u+v+1) q(x).
//
//  This script specializes the parameter s and asks Magma for the full
//  rational torsion subgroup.  A hit with invariants [2,2,4,8] would be
//  an example of the desired type (8,4,2,2).
//
//  Typical run from torsion_jac:
//      magma -b height:=80 code/m18_m2228_curve_extra_halving_search.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 80;
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
            r := Q!num/den;
            key := Sprint(r);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, r);
            end if;
        end for;
    end for;
    return vals;
end function;

function IsSquareQ(q)
    q := Q!q;
    if q lt 0 then
        return false;
    end if;
    return IsSquare(Numerator(q)) and IsSquare(Denominator(q));
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function FamilyPolynomial(s)
    u := ((2*s^2 - 2*s + 1)*(s^3 - s^2 + 2*s - 1))
         /((s^2 - 1)*(s^3 - 2*s^2 + s - 1));
    v := s*(s^2 - 2*s + 2)/(s^2 - 1);
    q := -(u^2 + u*v + v^2 + u + v + 1)*x^2
         + (v + 1)*(u + 1)*(u + v)*x
         - u*v*(u + v + 1);
    f := (x - 1)*(x - u)*(x - v)*(x + u + v + 1)*q;
    return f, q, u, v;
end function;

function IrreducibleFrobeniusCertificate(fI)
    C := HyperellipticCurve(fI);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61] do
        try
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
qsplit := 0;
base_torsion := 0;
hits := [];

for s in params do
    if s in {Q!-1, Q!0, Q!1} then
        continue;
    end if;
    try
        f, q, u, v := FamilyPolynomial(s);
    catch e
        continue;
    end try;
    checked +:= 1;
    if Degree(f) ne 6 or Discriminant(f) eq 0 then
        continue;
    end if;
    smooth +:= 1;

    discq := Discriminant(q);
    if discq ne 0 and IsSquareQ(discq) then
        qsplit +:= 1;
    end if;

    fI, L := IntegralModelPolynomial(f);
    C := HyperellipticCurve(fI);
    J := Jacobian(C);
    G, phi := TorsionSubgroup(J);
    invs := Invariants(G);
    if invs eq [ 2, 2, 2, 8 ] then
        base_torsion +:= 1;
    end if;

    if #invs ge 4 and invs[#invs] ge 8 and invs[#invs - 1] ge 4 then
        simple, pcert, Lp := IrreducibleFrobeniusCertificate(fI);
        Append(~hits, <s,u,v,invs,simple,pcert,Lp,f>);
        print "HIT", "s", s, "u", u, "v", v,
              "torsion", invs, "simple", simple, "pcert", pcert;
        print "  f =", f;
        if #hits ge max_hits then
            break;
        end if;
    end if;
end for;

print "DONE height", height;
print "checked", checked, "smooth", smooth, "qsplit", qsplit,
      "base_torsion", base_torsion, "hits", #hits;

quit;
