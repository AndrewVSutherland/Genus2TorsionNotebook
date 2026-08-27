//////////////////////////////////////////////////////////////////////
//  Search on the derived M_1(8,4) cover of M_1(8,2^w).
//
//  Set n=1 and R=m/n.  The first tangent-equation condition for
//  halving W_0 - infinity in
//
//      y^2 = x*A(x)*B(x)
//
//  is
//
//      (2*m^2 + m*n + 4*t)/(m*n + 2*n^2 + 4*t) = w^2.
//
//  Solving gives t in terms of R,w.  The remaining tangent equation
//  factors into two quadratics in the tangent coefficient U.  Their
//  discriminants are:
//
//  plus:
//      -4*(w-R)*(R-1)^2*(R+1)*(R+w)*(w+1)
//        *(R*w - 3*R + 3*w - 1)
//
//  minus:
//       4*(w-1)*(R+1)*(R*w + 3*R + 3*w + 1)
//        *(R^4 - 2*R^3 + R^2*w^2 - R^2 + 2*R*w^2 - w^2).
//
//  This script searches rational R,w for which one discriminant is a
//  square, then verifies exact divisibility in Magma.
//
//  Typical run from torsion_jac:
//      magma -b height:=20 code/m18_m14_cover_search.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 20;
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

function FamilyPolynomial(R, w)
    m := R;
    n := Q!1;
    t := (2*m^2 + (1-w^2)*m*n - 2*w^2*n^2)/(4*(w^2-1));
    A := n^4*x^2
         + (m^3*n + 4*m^2*t + m*n^3 - 8*m*n*t + 4*n^2*t)*x
         + m^4;
    B := (m*n + 2*n^2 + 4*t)*x^2
         + (m^2 + 4*m*n + n^2 + 8*t)*x
         + (2*m^2 + m*n + 4*t);
    return x*A*B, t, A, B;
end function;

function PlusDisc(R, w)
    return -4*(w-R)*(R-1)^2*(R+1)*(R+w)*(w+1)
            *(R*w - 3*R + 3*w - 1);
end function;

function MinusDisc(R, w)
    return 4*(w-1)*(R+1)*(R*w + 3*R + 3*w + 1)
           *(R^4 - 2*R^3 + R^2*w^2 - R^2 + 2*R*w^2 - w^2);
end function;

params := RationalParametersOfHeight(height);
checked := 0;
cover_points := 0;
verified := 0;
hits := [];

for R in params do
    for w in params do
        if w in {Q!-1, Q!0, Q!1} then
            continue;
        end if;
        if R eq 0 then
            continue;
        end if;
        checked +:= 1;

        plus_ok := IsSquareQ(PlusDisc(R,w));
        minus_ok := IsSquareQ(MinusDisc(R,w));
        if not (plus_ok or minus_ok) then
            continue;
        end if;
        cover_points +:= 1;

        f, t, A, B := FamilyPolynomial(R,w);
        if Degree(f) ne 5 or Discriminant(f) eq 0 then
            continue;
        end if;
        fI, L := IntegralModelPolynomial(f);
        C := HyperellipticCurve(fI);
        J := Jacobian(C);
        D := J![x, Q!0];
        divisible, half := IsDivisibleBy(D, 2);
        if not divisible then
            print "WARNING cover point did not verify", R, w,
                  "plus", plus_ok, "minus", minus_ok, "t", t;
            continue;
        end if;

        verified +:= 1;
        G, phi := TorsionSubgroup(J);
        invs := Invariants(G);
        Append(~hits, <R,w,t,plus_ok,minus_ok,invs,f>);
        print "HIT", "R", R, "w", w, "t", t,
              "plus", plus_ok, "minus", minus_ok,
              "half_order", Order(half), "torsion", invs;
        if #hits ge max_hits then
            break R;
        end if;
    end for;
end for;

print "DONE height", height;
print "checked", checked, "cover_points", cover_points,
      "verified", verified, "hits", #hits;

quit;
