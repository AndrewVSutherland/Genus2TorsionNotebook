//////////////////////////////////////////////////////////////////////
//  Parametrized direct search for M(2,2,4,8) in cover coordinates.
//
//  This improves code/m2248_rst_direct_search.m by parametrizing the
//  conditions that the reconstructed branch values C,D are rational
//  squares.  For fixed sigma and slope m, the conic
//      rho^2 = sigma^2(1-c^2) + c^2
//  is parametrized by the line through (rho,c)=(sigma,0):
//      rho = sigma*((sigma^2-1)m^2 - 1)/(1 + (sigma^2-1)m^2),
//      c   = -2*sigma*m/(1 + (sigma^2-1)m^2).
//  For fixed rho and slope n, the same formula gives tau,d.
//
//  Typical run from torsion_jac:
//      magma -b height:=8 max_hits:=20 code/m2248_rst_param_search.m
//////////////////////////////////////////////////////////////////////

load "code/m2248_sieve.m";

if not assigned height then
    height := 8;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;

if not assigned progress_interval then
    progress_interval := 500000;
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

function IsNonzeroSquareQ(q)
    if q eq 0 then
        return false, Q!0;
    end if;
    num := Z!Numerator(q);
    den := Z!Denominator(q);
    if num lt 0 then
        return false, Q!0;
    end if;
    okN, rtN := IsSquare(num);
    if not okN then
        return false, Q!0;
    end if;
    okD, rtD := IsSquare(den);
    if not okD then
        return false, Q!0;
    end if;
    return true, Q!rtN / Q!rtD;
end function;

function ConicParam(s, m)
    den := 1 + (s^2 - 1)*m^2;
    if den eq 0 then
        return false, Q!0, Q!0;
    end if;
    r := s*((s^2 - 1)*m^2 - 1)/den;
    c := -2*s*m/den;
    return true, r, c;
end function;

function IrreducibleFrobeniusCertificateFromTuple(tup)
    a := Q!tup[1]; b := Q!tup[2]; c := Q!tup[3]; d := Q!tup[4];
    f := x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2);
    C := HyperellipticCurve(f);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71] do
        try
            fp := ChangeRing(f, GF(p));
            if Discriminant(fp) eq 0 then
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

function TupleKey(tup)
    roots := Sort([ Q!tup[i]^2 : i in [1..4] ]);
    return Sprint(roots);
end function;

params := RationalParametersOfHeight(height);
print "Parametrized (rho,sigma,tau) M2248 search";
print "height", height, "parameters", #params;

checked := 0;
valid_param := 0;
cover_square := 0;
exact_sources := 0;
simple_hits := 0;
hits := [];
seen := {};
stop := false;

for sigma in params do
    if stop then
        break;
    end if;
    if sigma eq 0 or sigma eq -1 or sigma eq 1 then
        continue;
    end if;
    for m in params do
        if stop then
            break;
        end if;
        if m eq 0 then
            continue;
        end if;
        ok1, rho, c := ConicParam(sigma, m);
        if not ok1 then
            continue;
        end if;
        if rho eq 0 or rho eq -1 or rho eq -sigma or c eq 0 then
            continue;
        end if;

        for n in params do
            if n eq 0 then
                continue;
            end if;
            checked +:= 1;
            ok2, tau, d := ConicParam(rho, n);
            if not ok2 then
                continue;
            end if;
            if tau eq 0 or tau eq -1 or tau eq 1 or tau eq -rho or tau eq -sigma or d eq 0 then
                continue;
            end if;
            tup := [rho, Q!1, c, d];
            if #(Set([ t^2 : t in tup ])) ne 4 then
                continue;
            end if;
            valid_param +:= 1;

            F1 := (1 + rho)*(1 + sigma)*(1 + tau);
            F2 := rho*(1 + rho)*(rho + sigma)*(rho + tau);
            F3 := sigma*(1 + sigma)*(rho + sigma)*(sigma + tau);
            F4 := tau*(1 + tau)*(rho + tau)*(sigma + tau);
            okF1, q1 := IsNonzeroSquareQ(F1);
            if not okF1 then
                continue;
            end if;
            okF2, q2 := IsNonzeroSquareQ(F2);
            if not okF2 then
                continue;
            end if;
            okF3, q3 := IsNonzeroSquareQ(F3);
            if not okF3 then
                continue;
            end if;
            okF4, q4 := IsNonzeroSquareQ(F4);
            if not okF4 then
                continue;
            end if;
            cover_square +:= 1;

            key := TupleKey(tup);
            if key in seen then
                continue;
            end if;
            Include(~seen, key);

            witnesses := M2248WitnessesForTupleAllPermutations(tup, true);
            if #witnesses eq 0 then
                continue;
            end if;
            exact_sources +:= 1;
            simple, pcert, Lp := IrreducibleFrobeniusCertificateFromTuple(tup);
            if simple then
                simple_hits +:= 1;
            end if;
            Append(~hits, <rho,sigma,tau,m,n,tup,simple,pcert,Lp,#witnesses>);
            print "HIT", "rho", rho, "sigma", sigma, "tau", tau,
                  "m", m, "n", n, "tuple", tup,
                  "witnesses", #witnesses,
                  "simple", simple, "pcert", pcert, "L", Lp;
            if #hits ge max_hits then
                stop := true;
                break;
            end if;

            if progress_interval gt 0 and checked mod progress_interval eq 0 then
                print "checked", checked, "valid_param", valid_param,
                      "cover_square", cover_square, "hits", #hits;
            end if;
        end for;
    end for;
end for;

print "DONE";
print "checked", checked,
      "valid_param", valid_param,
      "cover_square", cover_square,
      "exact_sources", exact_sources,
      "simple_hits", simple_hits,
      "hits", #hits;
