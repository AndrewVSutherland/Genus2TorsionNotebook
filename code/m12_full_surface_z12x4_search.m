//////////////////////////////////////////////////////////////////////
//  Exact rational search for Z/12 x Z/4 on the full two-dimensional
//  extra-Weierstrass surface inside M(12).
//
//  Split T+1 by a=(1-z^2)/(4*(r+1)).  Each rational root u of
//  Q4 = W/(T+1) gives an independent rational 2-torsion class.  For
//  each root w of T+1, move w to infinity and test whether the class
//  [u-w] is divisible by 2 in the Jacobian.
//
//  Typical run from torsion_jac:
//      magma -b height:=15 code/m12_full_surface_z12x4_search.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 15;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
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

Q := Rationals();
P<x> := PolynomialRing(Q);
PX<X> := PolynomialRing(Q);

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

function M12Data(a, r)
    T := a*x^2 - x + r;
    h := (x-r)*(T+1);
    W := h^2 + 4*a*x^2*T*(T+1);
    Q4 := ExactQuotient(W, T+1);
    return W, T, h, Q4;
end function;

function OddQuinticAtRoot(W, w)
    out := PX!0;
    for i in [0..Degree(W)] do
        for j in [0..i] do
            out +:= Coefficient(W, i)*Binomial(i,j)*w^(i-j)*X^(6-j);
        end for;
    end for;
    return out;
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return PX!(L^2*f), L;
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
    return false, 0, PX!0;
end function;

print "M(12) full extra-Weierstrass Z/12 x Z/4 search";
print "height", height;

params := RationalParametersOfHeight(height);
hits := [];
checked_split := 0;
extra_roots := 0;
exact_tests := 0;

for r in params do
    for z in params do
        if #hits ge max_hits then
            break r;
        end if;
        if r eq -1 or z^2 eq 1 then
            continue;
        end if;
        a := (1-z^2)/(4*(r+1));
        if a eq 0 then
            continue;
        end if;

        W, T, h, Q4 := M12Data(a, r);
        if Degree(W) ne 6 or Discriminant(W) eq 0 then
            continue;
        end if;

        rootsT := Roots(T+1);
        if #rootsT lt 2 then
            continue;
        end if;
        checked_split +:= 1;

        rootsQ := [ rt[1] : rt in Roots(Q4) | rt[2] eq 1 ];
        if #rootsQ eq 0 then
            continue;
        end if;
        extra_roots +:= #rootsQ;

        for wd in rootsT do
            w := wd[1];
            if w eq 0 then
                continue;
            end if;

            f5 := OddQuinticAtRoot(W, w);
            if Degree(f5) ne 5 or Discriminant(f5) eq 0 then
                continue;
            end if;

            Y0 := Evaluate(h, 0);
            Xp := -1/w;
            Yp := Y0*Xp^3;
            fI, L := IntegralModelPolynomial(f5);
            if Discriminant(fI) eq 0 then
                continue;
            end if;

            C := HyperellipticCurve(fI);
            J := Jacobian(C);
            D := J![X-Xp, L*Yp];
            try
                ord := Order(D);
            catch e
                continue;
            end try;
            if ord ne 12 then
                continue;
            end if;
            Tdiv := 6*D;

            for u in rootsQ do
                if u eq w then
                    continue;
                end if;
                beta := 1/(u-w);
                Tbeta := J![X-beta, Q!0];
                if Tbeta eq J!0 or Tbeta eq Tdiv then
                    continue;
                end if;

                exact_tests +:= 1;
                if IsDivisibleBy(Tbeta, 2) then
                    simple, pcert, Lp := IrreducibleFrobeniusCertificate(f5);
                    Append(~hits, <r,z,a,w,u,beta,f5,simple,pcert,Lp>);
                    print "HIT";
                    print "  r,z,a", r, z, a;
                    print "  w,u,beta", w, u, beta;
                    print "  simple", simple, "prime", pcert, "L", Lp;
                    print "  f5", f5;
                end if;
            end for;
        end for;

        if progress_interval gt 0 and checked_split mod progress_interval eq 0 then
            print "checked_split", checked_split, "extra_roots", extra_roots,
                  "exact_tests", exact_tests, "hits", #hits;
        end if;
    end for;
end for;

print "Done";
print "checked_split", checked_split;
print "extra_roots", extra_roots;
print "exact_tests", exact_tests;
print "hits", #hits;
for H in hits do
    print H;
end for;
quit;
