//////////////////////////////////////////////////////////////////////
//  B-split cover search inside the M_1(8,4) tangent locus.
//
//  In the n=1 chart, the B quadratic has normalized form
//
//      x^2 + S*x + w^2.
//
//  Requiring it to split is parameterized by S = w*(z + 1/z), giving
//
//      R = (-3*w^2*z + 2*w*z^2 + 2*w - z)
//          /(w^2*z - 2*w*z^2 - 2*w + 3*z).
//
//  This script searches small rational (w,z), then imposes A split and
//  the exact tangent equation for halving W_0 - infinity.
//
//  Typical run from torsion_jac:
//      magma -b height:=25 code/m18_m14_bsplit_cover_search.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 25;
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

function BSplitR(w, z)
    den := w^2*z - 2*w*z^2 - 2*w + 3*z;
    if den eq 0 then
        return false, Q!0;
    end if;
    R := (-3*w^2*z + 2*w*z^2 + 2*w - z)/den;
    return true, R;
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

function TangentCandidates(f, R, w)
    h := ExactQuotient(f, x);
    c0 := Coefficient(h, 0);
    c1 := Coefficient(h, 1);
    c2 := Coefficient(h, 2);
    c3 := Coefficient(h, 3);
    c4 := Coefficient(h, 4);

    out := [];
    PR<U> := PolynomialRing(Q);

    for V in [R^2*w, -R^2*w] do
        F := 4*(c3 - 2*c4*U)*(c1 - 2*c4*U*V)
             - (c2 - c4*(U^2 + 2*V))^2;
        for rt in Roots(F) do
            U0 := rt[1];
            M2 := c3 - 2*c4*U0;
            N2 := c1 - 2*c4*U0*V;
            if IsSquareQ(M2) and IsSquareQ(N2) then
                Append(~out, <U0, V, M2, N2>);
            end if;
        end for;
    end for;
    return out;
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
valid_chart := 0;
a_split := 0;
tangent_points := 0;
verified := 0;
hits := [];

for w in params do
    if w in {Q!-1, Q!0, Q!1} then
        continue;
    end if;
    for z in params do
        if z eq 0 then
            continue;
        end if;
        checked +:= 1;
        okR, R := BSplitR(w, z);
        if not okR or R eq 0 then
            continue;
        end if;
        valid_chart +:= 1;

        f, t, A, B := FamilyPolynomial(R,w);
        if Degree(f) ne 5 or Discriminant(f) eq 0 then
            continue;
        end if;

        discA := Discriminant(A);
        discB := Discriminant(B);
        if discA eq 0 or discB eq 0 then
            continue;
        end if;
        if not (IsSquareQ(discA) and IsSquareQ(discB)) then
            continue;
        end if;
        a_split +:= 1;

        candidates := TangentCandidates(f, R, w);
        if #candidates eq 0 then
            continue;
        end if;
        tangent_points +:= 1;

        fI, L := IntegralModelPolynomial(f);
        C := HyperellipticCurve(fI);
        J := Jacobian(C);
        D := J![x, Q!0];
        divisible, half := IsDivisibleBy(D, 2);
        if not divisible then
            print "WARNING tangent equations did not verify", R, w, z, t, candidates;
            continue;
        end if;

        verified +:= 1;
        G, phi := TorsionSubgroup(J);
        invs := Invariants(G);
        simple, pcert, Lp := IrreducibleFrobeniusCertificate(fI);
        Append(~hits, <R,w,z,t,invs,simple,pcert,Lp,f,discA,discB>);
        print "HIT", "R", R, "w", w, "z", z, "t", t,
              "discA", discA, "discB", discB,
              "candidates", #candidates, "half_order", Order(half),
              "torsion", invs, "simple", simple, "pcert", pcert;
        print "  f =", f;
        if #hits ge max_hits then
            break w;
        end if;
    end for;
end for;

print "DONE height", height;
print "checked", checked, "valid_chart", valid_chart,
      "a_split", a_split, "tangent_points", tangent_points,
      "verified", verified, "hits", #hits;

quit;
