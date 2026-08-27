//////////////////////////////////////////////////////////////////////
//  Direct search in the M(2,2,4,8) cover coordinates (rho,sigma,tau).
//
//  In code/m2248_sieve.m, a full witness for an ordered tuple [a,b,c,d]
//  is encoded by
//      rho = +/- a/b,
//      sigma = +/- (A-C)/sqrt((A-C)(B-C)),
//      tau = +/- (A-D)/sqrt((A-D)(B-D)),
//  where A=a^2, B=b^2, C=c^2, D=d^2.
//
//  Conversely, set B=1, A=rho^2.  Then the sigma and tau equations force
//      C = (sigma^2 - A)/(sigma^2 - 1),
//      D = (tau^2   - A)/(tau^2   - 1).
//  If C and D are rational squares and the four full-cover expressions
//  are squares, [rho,1,sqrt(C),sqrt(D)] is an ordered M(2,2,4,8)
//  candidate.  We then verify with the existing exact sieve and certify
//  Q-simplicity by irreducible Frobenius polynomial when possible.
//
//  Typical run from torsion_jac:
//      magma -b height:=12 max_hits:=20 code/m2248_rst_direct_search.m
//////////////////////////////////////////////////////////////////////

load "code/m2248_sieve.m";

if not assigned height then
    height := 10;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;

if not assigned progress_interval then
    progress_interval := 200000;
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
print "Direct (rho,sigma,tau) M2248 search";
print "height", height, "parameters", #params;

checked := 0;
cover_square := 0;
cd_square := 0;
exact_sources := 0;
full_sources := 0;
simple_hits := 0;
hits := [];
seen := {};
stop := false;

for rho in params do
    if stop then
        break;
    end if;
    if rho eq 0 or rho eq -1 then
        continue;
    end if;
    A := rho^2;

    for sigma in params do
        if stop then
            break;
        end if;
        if sigma eq 0 or sigma eq -1 or sigma eq 1 or sigma eq -rho then
            continue;
        end if;
        sig2 := sigma^2;
        Cval := (sig2 - A)/(sig2 - 1);
        okC, c := IsNonzeroSquareQ(Cval);
        if not okC then
            continue;
        end if;

        for tau in params do
            if tau eq 0 or tau eq -1 or tau eq 1 or tau eq -rho or tau eq -sigma then
                continue;
            end if;
            checked +:= 1;

            F1 := (1 + rho)*(1 + sigma)*(1 + tau);
            F2 := rho*(1 + rho)*(rho + sigma)*(rho + tau);
            F3 := sigma*(1 + sigma)*(rho + sigma)*(sigma + tau);
            F4 := tau*(1 + tau)*(rho + tau)*(sigma + tau);
            ok1, q1 := IsNonzeroSquareQ(F1);
            if not ok1 then
                continue;
            end if;
            ok2, q2 := IsNonzeroSquareQ(F2);
            if not ok2 then
                continue;
            end if;
            ok3, q3 := IsNonzeroSquareQ(F3);
            if not ok3 then
                continue;
            end if;
            ok4, q4 := IsNonzeroSquareQ(F4);
            if not ok4 then
                continue;
            end if;
            cover_square +:= 1;

            tau2 := tau^2;
            Dval := (tau2 - A)/(tau2 - 1);
            okD, d := IsNonzeroSquareQ(Dval);
            if not okD then
                continue;
            end if;
            tup := [rho, Q!1, c, d];
            if #(Set([ t^2 : t in tup ])) ne 4 then
                continue;
            end if;
            cd_square +:= 1;

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
            full_sources +:= 1;

            simple, pcert, Lp := IrreducibleFrobeniusCertificateFromTuple(tup);
            if simple then
                simple_hits +:= 1;
            end if;
            Append(~hits, <rho,sigma,tau,tup,simple,pcert,Lp,#witnesses>);
            print "HIT", "rho", rho, "sigma", sigma, "tau", tau,
                  "tuple", tup, "witnesses", #witnesses,
                  "simple", simple, "pcert", pcert, "L", Lp;
            if #hits ge max_hits then
                stop := true;
                break;
            end if;

            if progress_interval gt 0 and checked mod progress_interval eq 0 then
                print "checked", checked, "cover_square", cover_square,
                      "cd_square", cd_square, "hits", #hits;
            end if;
        end for;
    end for;
end for;

print "DONE";
print "checked", checked,
      "cover_square", cover_square,
      "cd_square", cd_square,
      "exact_sources", exact_sources,
      "full_sources", full_sources,
      "simple_hits", simple_hits,
      "hits", #hits;
