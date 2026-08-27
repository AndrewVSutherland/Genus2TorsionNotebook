//////////////////////////////////////////////////////////////////////
//  Finite-field scan for the full-rational-Weierstrass condition inside
//  M(12).  This tests the completed-square polynomial
//      W = ((x-r)(T+1))^2 + 4*a*x^2*T*(T+1), T=a*x^2-x+r,
//  and records good parameter points where W has six distinct roots in F_p.
//
//  Typical run from torsion_jac:
//      magma code/m12_full_split_finite_field_sieve.m
//////////////////////////////////////////////////////////////////////

prime_list := [3,5,7,11,13,17,19,23,29,31,37,41,43,47];

function ProjectiveOrAffineParams(F)
    // M(12) is affine A^2_{a,r}; scan affine points mod p.
    return [ [a,r] : a in F, r in F ];
end function;

function WPolynomial(F, vals)
    P<x> := PolynomialRing(F);
    a := vals[1];
    r := vals[2];
    T := a*x^2 - x + r;
    h := (x-r)*(T+1);
    f := a*x^2*T*(T+1);
    W := h^2 + 4*f;
    return W;
end function;

function IsGoodFullSplit(W)
    if Degree(W) ne 6 or Discriminant(W) eq 0 then
        return false;
    end if;
    roots := Roots(W);
    return #roots eq 6 and &and [ rt[2] eq 1 : rt in roots ];
end function;

for p in prime_list do
    F := GF(p);
    good := 0;
    split := [];
    for vals in ProjectiveOrAffineParams(F) do
        W := WPolynomial(F, vals);
        if Degree(W) eq 6 and Discriminant(W) ne 0 then
            good +:= 1;
            if IsGoodFullSplit(W) then
                Append(~split, vals);
            end if;
        end if;
    end for;
    print "p", p, "good", good, "full_split", #split;
    if #split gt 0 then
        print "  examples", split[1..Minimum(#split, 10)];
    end if;
    if good gt 0 and #split eq 0 then
        print "  NO GOOD FULL-SPLIT POINTS MOD", p;
        break;
    end if;
end for;
quit;
