//////////////////////////////////////////////////////////////////////
//  Halve the extra rational 2-torsion in the [2,4,8] one-split family.
//
//  Base model:
//      C: y^2 = f(x) = x*A(x)*B(x)
//  on the M_1(8,4) tangent locus.  If exactly one of A,B splits over Q,
//  the torsion is typically [2,4,8].
//
//  Algebraic halving condition for a split branch point alpha:
//      f(x) = (x-alpha)*h_alpha(x).
//  The class (alpha,0)-infinity is divisible by 2 iff there are rational
//  U,V,M,N such that, after shifting alpha to 0,
//
//      h_alpha(z) - z*(M*z+N)^2 = lc(h_alpha)*(z^2+U*z+V)^2.
//
//  This script searches [2,4,8] specializations and imposes this tangent
//  system before exact Jacobian verification.
//
//  Typical run:
//      magma -b height:=50 max_hits:=20 \
//          code/m18_m14_248_halve_split_root_search.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned height then
    height := 40;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;

if not assigned progress_interval then
    progress_interval := 1000000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

if not assigned require_exact_one_split then
    require_exact_one_split := true;
elif Type(require_exact_one_split) eq MonStgElt then
    require_exact_one_split := require_exact_one_split in {"true", "True", "1", "yes"};
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

function SquareRootsQ(q)
    q := Q!q;
    if q eq 0 then
        return [Q!0];
    end if;
    if not IsSquareQ(q) then
        return [];
    end if;
    okn, rn := IsSquare(Numerator(q));
    okd, rd := IsSquare(Denominator(q));
    assert okn and okd;
    r := Q!rn/Q!rd;
    return [r, -r];
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function GoodHyperellipticPolynomial(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function FamilyPolynomial(R, w)
    t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
    A := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
    B := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x
         + (2*R^2 + R + 4*t);
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

function M14CoverPossible(R, w)
    return IsSquareQ(PlusDisc(R,w)) or IsSquareQ(MinusDisc(R,w));
end function;

function TangentCandidates(f, R, w)
    h := ExactQuotient(f, x);
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

function RootHalvingCandidates(f, alpha)
    // Shift x=alpha+z, so alpha is the root z=0.
    g := Evaluate(f, x + alpha);
    if Coefficient(g, 0) ne 0 then
        return [];
    end if;
    h := ExactQuotient(g, x);
    c0 := Coefficient(h, 0);
    c1 := Coefficient(h, 1);
    c2 := Coefficient(h, 2);
    c3 := Coefficient(h, 3);
    c4 := Coefficient(h, 4);
    if c4 eq 0 then
        return [];
    end if;

    out := [];
    PR<U> := PolynomialRing(Q);
    for V in SquareRootsQ(c0/c4) do
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

function RationalRoots(poly)
    return [rt[1] : rt in Roots(poly)];
end function;

function TorsionOrder(invs)
    if #invs eq 0 then
        return 1;
    end if;
    return &*invs;
end function;

function Exponent(invs)
    if #invs eq 0 then
        return 1;
    end if;
    return invs[#invs];
end function;

function IrreducibleFrobeniusCertificate(fI)
    C := HyperellipticCurve(fI);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79] do
        try
            fp := ChangeRing(fI, GF(p));
            if not GoodHyperellipticPolynomial(fp) then
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

params := RationalParametersOfHeight(height);
checked := 0;
cover := 0;
smooth := 0;
one_split := 0;
first_tangent := 0;
root_tests := 0;
root_candidates := 0;
exact_verified := 0;
hits := [];
torsion_counts := AssociativeArray();

print "Halve split-root 2-torsion in the [2,4,8] one-split family";
print "height", height, "params", #params, "pairs", #params^2,
      "max_hits", max_hits, "progress_interval", progress_interval,
      "require_exact_one_split", require_exact_one_split;

for R in params do
    for w in params do
        if #hits ge max_hits then
            break R;
        end if;
        if R eq 0 or w in {Q!-1, Q!0, Q!1} then
            continue;
        end if;
        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "progress", checked, "cover", cover, "smooth", smooth,
                  "one_split", one_split, "first_tangent", first_tangent,
                  "root_tests", root_tests,
                  "root_candidates", root_candidates,
                  "exact_verified", exact_verified,
                  "hits", #hits;
        end if;

        if not M14CoverPossible(R,w) then
            continue;
        end if;
        cover +:= 1;

        f, t, A, B := FamilyPolynomial(R,w);
        if not GoodHyperellipticPolynomial(f) then
            continue;
        end if;
        smooth +:= 1;

        discA := Discriminant(A);
        discB := Discriminant(B);
        if discA eq 0 or discB eq 0 then
            continue;
        end if;
        splitA := IsSquareQ(discA);
        splitB := IsSquareQ(discB);
        if not (splitA or splitB) then
            continue;
        end if;
        if require_exact_one_split and splitA and splitB then
            continue;
        end if;
        one_split +:= 1;

        first := TangentCandidates(f, R, w);
        if #first eq 0 then
            continue;
        end if;
        first_tangent +:= 1;

        split_roots := [];
        if splitA then
            for alpha in RationalRoots(A) do
                Append(~split_roots, <"A", alpha>);
            end for;
        end if;
        if splitB then
            for alpha in RationalRoots(B) do
                Append(~split_roots, <"B", alpha>);
            end for;
        end if;

        for item in split_roots do
            block := item[1];
            alpha := item[2];
            root_tests +:= 1;
            rhalves := RootHalvingCandidates(f, alpha);
            if #rhalves eq 0 then
                continue;
            end if;
            root_candidates +:= 1;

            fI, L := IntegralModelPolynomial(f);
            C := HyperellipticCurve(fI);
            J := Jacobian(C);
            Dalpha := J![x-alpha, Q!0];
            ok, half := IsDivisibleBy(Dalpha, 2);
            if not ok then
                print "WARNING tangent root-halving equations did not verify",
                      "R", R, "w", w, "t", t,
                      "block", block, "alpha", alpha,
                      "candidates", rhalves;
                continue;
            end if;
            exact_verified +:= 1;

            G, phi := TorsionSubgroup(J);
            invs := Invariants(G);
            key := Sprint(invs);
            if IsDefined(torsion_counts, key) then
                torsion_counts[key] +:= 1;
            else
                torsion_counts[key] := 1;
            end if;
            simple, pcert, Lp := IrreducibleFrobeniusCertificate(fI);
            Append(~hits, <R,w,t,block,alpha,invs,simple,pcert,Lp,fI,rhalves[1]>);
            print "HIT", "R", R, "w", w, "t", t,
                  "block", block, "alpha", alpha,
                  "root_halves", #rhalves,
                  "half_order", Order(half),
                  "torsion", invs,
                  "order", TorsionOrder(invs),
                  "exponent", Exponent(invs),
                  "simple", simple, "pcert", pcert;
            print "  first_cover_candidates", #first;
            print "  root_halving_candidate", rhalves[1];
            print "  fI", fI;
            if #hits ge max_hits then
                break R;
            end if;
        end for;
    end for;
end for;

print "DONE height", height;
print "checked", checked, "cover", cover, "smooth", smooth,
      "one_split", one_split, "first_tangent", first_tangent,
      "root_tests", root_tests, "root_candidates", root_candidates,
      "exact_verified", exact_verified,
      "hits", #hits;
print "TORSION_COUNTS";
for key in Sort([k : k in Keys(torsion_counts)]) do
    print " ", key, torsion_counts[key];
end for;

quit;
