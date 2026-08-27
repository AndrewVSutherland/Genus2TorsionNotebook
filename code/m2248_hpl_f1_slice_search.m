//////////////////////////////////////////////////////////////////////
//  HPL-based F1-slice search for M(2,2,4,8).
//
//  Start from a verified HPL full-cover witness (rho0,sigma0,tau0).
//  Write F1=(1+rho)(1+sigma)(1+tau).  Along the slice
//      sigma = sigma0*(1 + ds*u),
//      tau   = tau0  *(1 + dt*u),
//      q1    = q10   *(1 + dq*u),
//      rho   = q1^2/((1+sigma)(1+tau)) - 1,
//  the first full-cover square condition F1=q1^2 is automatic.
//
//  The script searches small integer directions (ds,dt,dq) and rational
//  u, then tests the remaining M2248 full-cover conditions, the branch
//  reconstruction C=c^2 and D=d^2, the exact m2248_sieve witness, and
//  a Frobenius irreducibility certificate for Q-simplicity.
//
//  Typical run from torsion_jac:
//      magma -b dir_bound:=3 height:=20 max_hits:=20 \
//          code/m2248_hpl_f1_slice_search.m
//////////////////////////////////////////////////////////////////////

load "code/m2248_sieve.m";

if not assigned dir_bound then
    dir_bound := 3;
elif Type(dir_bound) eq MonStgElt then
    dir_bound := StringToInteger(dir_bound);
end if;

if not assigned height then
    height := 20;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;

if not assigned progress_interval then
    progress_interval := 50000;
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
            if num eq 0 or GCD(num, den) ne 1 then
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

function CoverValues(rho, sigma, tau)
    F0 := rho*sigma*tau;
    F1 := (1 + rho)*(1 + sigma)*(1 + tau);
    F2 := rho*(1 + rho)*(rho + sigma)*(rho + tau);
    F3 := sigma*(1 + sigma)*(rho + sigma)*(sigma + tau);
    F4 := tau*(1 + tau)*(rho + tau)*(sigma + tau);
    return F0, F1, F2, F3, F4;
end function;

// One full HPL witness from code/m2248_hpl_normalization_check.m.
rho0 := Q!58466134224 / Q!53109477625;
sigma0 := Q!719363573659505664 / Q!749082246897952705;
tau0 := Q!307598400 / Q!352612321;
F00, F10, F20, F30, F40 := CoverValues(rho0, sigma0, tau0);
ok10, q10 := IsNonzeroSquareQ(F10);
assert ok10;
assert #M2248WitnessesForTupleAllPermutations([rho0, Q!1, Q!1211794463/Q!2124379105, Q!406014821/Q!506569821], true) gt 0;

params := RationalParametersOfHeight(height);
print "HPL F1-slice search";
print "dir_bound", dir_bound, "height", height, "parameters", #params;
print "base", "rho", rho0, "sigma", sigma0, "tau", tau0, "q1", q10;

checked := 0;
valid_slice := 0;
branch_square := 0;
cover_square := 0;
exact_sources := 0;
simple_hits := 0;
hits := [];
seen := {};
stop := false;

for ds in [-dir_bound..dir_bound] do
    if stop then break; end if;
    for dt in [-dir_bound..dir_bound] do
        if stop then break; end if;
        for dq in [-dir_bound..dir_bound] do
            if stop then break; end if;
            if ds eq 0 and dt eq 0 and dq eq 0 then
                continue;
            end if;

            for u in params do
                checked +:= 1;
                sf := 1 + Q!ds*u;
                tf := 1 + Q!dt*u;
                qf := 1 + Q!dq*u;
                if sf eq 0 or tf eq 0 or qf eq 0 then
                    continue;
                end if;

                sigma := sigma0*sf;
                tau := tau0*tf;
                q1 := q10*qf;
                den := (1 + sigma)*(1 + tau);
                if den eq 0 then
                    continue;
                end if;
                rho := q1^2/den - 1;
                if rho eq 0 or sigma eq 0 or tau eq 0 then
                    continue;
                end if;
                if rho eq rho0 and sigma eq sigma0 and tau eq tau0 then
                    continue;
                end if;
                if rho eq -1 or sigma eq -1 or tau eq -1 then
                    continue;
                end if;
                if sigma^2 eq 1 or tau^2 eq 1 then
                    continue;
                end if;
                valid_slice +:= 1;

                Cval := (sigma^2 - rho^2)/(sigma^2 - 1);
                Dval := (tau^2 - rho^2)/(tau^2 - 1);
                okC, c := IsNonzeroSquareQ(Cval);
                if not okC then
                    continue;
                end if;
                okD, d := IsNonzeroSquareQ(Dval);
                if not okD then
                    continue;
                end if;
                tup := [rho, Q!1, c, d];
                if #(Set([ t^2 : t in tup ])) ne 4 then
                    continue;
                end if;
                branch_square +:= 1;

                F0, F1, F2, F3, F4 := CoverValues(rho, sigma, tau);
                ok0, q0 := IsNonzeroSquareQ(F0);
                if not ok0 then continue; end if;
                ok2, q2 := IsNonzeroSquareQ(F2);
                if not ok2 then continue; end if;
                ok3, q3 := IsNonzeroSquareQ(F3);
                if not ok3 then continue; end if;
                ok4, q4 := IsNonzeroSquareQ(F4);
                if not ok4 then continue; end if;
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
                Append(~hits, <ds,dt,dq,u,rho,sigma,tau,tup,simple,pcert,Lp,#witnesses>);
                print "HIT", "dir", [ds,dt,dq], "u", u,
                      "rho", rho, "sigma", sigma, "tau", tau,
                      "tuple", tup, "witnesses", #witnesses,
                      "simple", simple, "pcert", pcert, "L", Lp;
                if #hits ge max_hits then
                    stop := true;
                    break;
                end if;
            end for;

            if progress_interval gt 0 and checked mod progress_interval eq 0 then
                print "checked", checked, "valid_slice", valid_slice,
                      "branch_square", branch_square,
                      "cover_square", cover_square,
                      "hits", #hits;
            end if;
        end for;
    end for;
end for;

print "DONE";
print "checked", checked,
      "valid_slice", valid_slice,
      "branch_square", branch_square,
      "cover_square", cover_square,
      "exact_sources", exact_sources,
      "simple_hits", simple_hits,
      "hits", #hits;
