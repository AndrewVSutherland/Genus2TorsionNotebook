//////////////////////////////////////////////////////////////////////
//  Exact search for halving the independent rational 2-torsion point
//  in the one-parameter M(12) family
//
//      a = (1-r)/4.
//
//  The family has an order-12 point and an independent rational
//  2-torsion class, hence at least Z/12 x Z/2.  This script tests
//  whether the independent 2-torsion class is divisible by 2, which
//  would give Z/12 x Z/4.
//
//  Magma's IsDivisibleBy for hyperelliptic Jacobians requires an
//  integral model, so for y^2=f(x) we use the isomorphic model
//  Y^2=L^2 f(x), with L clearing denominators.
//
//  Typical run from torsion_jac:
//      magma -b height:=40 code/m12_z12x2_halving_search.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 40;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned progress_interval then
    progress_interval := 500;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Q := Rationals();
P<X> := PolynomialRing(Q);

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

function OddQuinticForR(r)
    a := (1-r)/4;
    T := a*X^2 - X + r;
    h := (X-r)*(T+1);
    W := h^2 + 4*a*X^2*T*(T+1);

    // Move the root x=2 of T+1 to infinity:
    // x = 2 + 1/X,  Y_old = Y_new/X^3.
    f5 := P!0;
    for i in [0..Degree(W)] do
        for j in [0..i] do
            f5 +:= Coefficient(W, i)*Binomial(i,j)*2^(i-j)*X^(6-j);
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

function IrreducibleFrobeniusCertificate(f)
    C := HyperellipticCurve(f);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47] do
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

print "M(12) Z/12 x Z/4 halving search on a=(1-r)/4";
print "height", height;

hits := [];
checked := 0;
good := 0;
known_ok := 0;

for r in RationalParametersOfHeight(height) do
    if r in {Q!0, Q!1, Q!2} then
        continue;
    end if;

    f5 := OddQuinticForR(r);
    if Degree(f5) ne 5 or Discriminant(f5) eq 0 then
        continue;
    end if;

    beta_ind := (2-r)/(4*(r-1));
    beta_div := (1-r)/(4*r);

    fI, L := IntegralModelPolynomial(f5);
    if Discriminant(fI) eq 0 then
        continue;
    end if;

    checked +:= 1;
    C := HyperellipticCurve(fI);
    J := Jacobian(C);
    Tind := J![X-beta_ind, Q!0];
    Tdiv := J![X-beta_div, Q!0];

    if Tind eq J!0 or Tind eq Tdiv then
        continue;
    end if;
    good +:= 1;

    // Sanity check: beta_div is the 2-torsion class 6D, hence divisible.
    if IsDivisibleBy(Tdiv, 2) then
        known_ok +:= 1;
    else
        print "WARNING known divisible class failed at r", r;
    end if;

    if IsDivisibleBy(Tind, 2) then
        simple, pcert, Lp := IrreducibleFrobeniusCertificate(f5);
        Append(~hits, <r, f5, simple, pcert, Lp>);
        print "HIT";
        print "  r", r;
        print "  f5", f5;
        print "  simple", simple, "prime", pcert, "L", Lp;
    end if;

    if progress_interval gt 0 and checked mod progress_interval eq 0 then
        print "checked", checked, "good", good, "known_ok", known_ok, "hits", #hits;
    end if;
end for;

print "Done";
print "checked", checked;
print "good", good;
print "known_ok", known_ok;
print "hits", #hits;
for H in hits do
    print H;
end for;
quit;
