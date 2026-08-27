//////////////////////////////////////////////////////////////////////
//  Search M(12) for strictly simple Jacobians.
//
//  Start with the model
//      y^2 + (x-r)(T+1)y = a x^2 T(T+1),  T=a*x^2-x+r.
//  Completing the square gives Y^2 = W(x).  We impose that the quadratic
//  factor T+1 has a rational root by writing
//      1 - 4*a*(r+1) = z^2,
//      a = (1-z^2)/(4*(r+1)).
//
//  For each rational root w of T+1, move w to infinity to get an odd
//  quintic.  The point P=(0,0) maps to P' and should give P'-infinity
//  of order 12.  We then look for a prime where the L-polynomial is
//  irreducible, which certifies that the Jacobian is Q-simple.
//
//  Typical run from torsion_jac:
//      magma -b height:=12 max_hits:=20 code/m12_simple_search.m
//////////////////////////////////////////////////////////////////////

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
if not assigned require_full_split then
    require_full_split := false;
elif Type(require_full_split) eq MonStgElt then
    require_full_split := require_full_split eq "true";
end if;
if not assigned min_rational_roots then
    min_rational_roots := 0;
elif Type(min_rational_roots) eq MonStgElt then
    min_rational_roots := StringToInteger(min_rational_roots);
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

function M12Polynomial(a, r)
    T := a*x^2 - x + r;
    h := (x-r)*(T+1);
    f := a*x^2*T*(T+1);
    W := h^2 + 4*f;
    return W, T, h;
end function;

function OddQuinticAtRoot(W, w)
    // x = w + 1/X, Y_old = Y_new/X^3.
    coeffs := [ Coefficient(W, i) : i in [0..Degree(W)] ];
    out := PX!0;
    for i in [0..#coeffs-1] do
        for j in [0..i] do
            out +:= coeffs[i+1] * Binomial(i,j) * w^(i-j) * X^(6-j);
        end for;
    end for;
    // W has a root at w, so the constant term vanishes and the degree is <=5.
    return out;
end function;

function IrreducibleFrobeniusCertificate(f5)
    C := HyperellipticCurve(f5);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47] do
        try
            Cp := ChangeRing(C, GF(p));
            Lp := LPolynomial(Cp);
            fac := Factorization(Lp);
            if #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4 then
                return true, p, Lp;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false, 0, Parent(X)!0;
end function;

function RationalRootCount(W)
    roots := Roots(W);
    return &+[ rt[2] : rt in roots ];
end function;

function FullSplit(W)
    return RationalRootCount(W) eq 6;
end function;

function HasIndependentRationalTwoTorsion(J, f5, D)
    for rt in Roots(f5) do
        if rt[2] ne 1 then
            continue;
        end if;
        beta := rt[1];
        Tbeta := J![X-beta, Q!0];
        if Tbeta ne J!0 and Tbeta ne 6*D then
            return true, beta;
        end if;
    end for;
    return false, Q!0;
end function;

print "M(12) simple search";
print "height", height;
print "max_hits", max_hits;
print "require_full_split", require_full_split;
print "min_rational_roots", min_rational_roots;

params := RationalParametersOfHeight(height);
hits := [];
seen := {};
checked := 0;

for z in params do
    for r in params do
        if #hits ge max_hits then
            break z;
        end if;
        if r eq -1 or z^2 eq 1 then
            continue;
        end if;
        a := (1-z^2)/(4*(r+1));
        if a eq 0 then
            continue;
        end if;
        W, T, h := M12Polynomial(a, r);
        if Degree(W) ne 6 or Discriminant(W) eq 0 then
            continue;
        end if;
        root_count := RationalRootCount(W);
        if min_rational_roots gt 0 and root_count lt min_rational_roots then
            continue;
        end if;
        if require_full_split and root_count lt 6 then
            continue;
        end if;

        roots_quad := Roots(T+1);
        if #roots_quad eq 0 then
            continue;
        end if;

        checked +:= 1;
        for root_data in roots_quad do
            if root_data[2] ne 1 then
                continue;
            end if;
            w := root_data[1];
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
            if Evaluate(f5, Xp) ne Yp^2 then
                print "internal point mismatch", a, r, z, w;
                continue;
            end if;

            C5 := HyperellipticCurve(f5);
            J := Jacobian(C5);
            D := J![X-Xp, Yp];
            ord := Order(D);
            if ord ne 12 then
                continue;
            end if;

            simple, pcert, Lp := IrreducibleFrobeniusCertificate(f5);
            if not simple then
                continue;
            end if;
            independent2, independent_beta := HasIndependentRationalTwoTorsion(J, f5, D);

            key := Sprint([a,r,w]);
            if key in seen then
                continue;
            end if;
            Include(~seen, key);
            Append(~hits, <a, r, z, w, f5, root_count, independent2, independent_beta, FullSplit(W), pcert, Lp>);
            print "HIT";
            print "  a,r,z,w", a, r, z, w;
            print "  rational_root_count", root_count, "independent2", independent2, "independent_beta", independent_beta, "full_split", FullSplit(W);
            print "  f5", f5;
            print "  simplicity_prime", pcert, "L", Lp;
        end for;
    end for;
end for;

print "Done";
print "checked parameter pairs with split T+1", checked;
print "hits", #hits;
for H in hits do
    print H;
end for;
quit;
