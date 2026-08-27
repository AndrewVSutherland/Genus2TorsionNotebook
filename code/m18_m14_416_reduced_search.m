//////////////////////////////////////////////////////////////////////
//  Targeted search for simple [4,16] examples in the M_1(8,4)
//  tangent-cover family.
//
//  This uses the reduced [4,16] necessary condition
//
//      ss^2 = -c4*R,        c4 = lc(A*B),
//
//  before doing exact Jacobian arithmetic.  A surviving specialization is
//  then checked by exact divisibility of both T_x and P_R, and finally by
//  exact TorsionSubgroup plus a fast irreducible Frobenius certificate.
//
//  Typical run:
//      magma -b height:=80 max_pr_tests:=200 \
//          code/m18_m14_416_reduced_search.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned height then
    height := 60;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

if not assigned max_pr_tests then
    max_pr_tests := 200;
elif Type(max_pr_tests) eq MonStgElt then
    max_pr_tests := StringToInteger(max_pr_tests);
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

function M14CoverPossible(R, w)
    return IsSquareQ(PlusDisc(R,w)) or IsSquareQ(MinusDisc(R,w));
end function;

function C4Value(R, w)
    return (2*R^2 - 2)/(w^2 - 1);
end function;

function Square416Possible(R, w)
    return IsSquareQ(-C4Value(R,w)*R);
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

function Monic(g)
    return g/LeadingCoefficient(g);
end function;

function YRValue(R, w)
    Qfac := R^2 - (Q!1/2)*R*w^2 + (Q!1/2)*R - w^2;
    return -2*R*(R-1)^2*Qfac/(w^2-1);
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

function TwoAdicExponents(invs)
    vals := Sort([Valuation(n, 2) : n in invs]);
    return Reverse(vals);
end function;

function Has416(invs)
    vals := TwoAdicExponents(invs);
    return #vals ge 2 and vals[1] ge 4 and vals[2] ge 2;
end function;

function IrreducibleFrobeniusCertificate(fI)
    C := HyperellipticCurve(fI);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97] do
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
    return false, 0, P!0;
end function;

params := RationalParametersOfHeight(height);
checked := 0;
square416 := 0;
cover := 0;
smooth := 0;
tangent_points := 0;
tx_halves := 0;
pr_tests := 0;
pr_halves := 0;
torsion_tests := 0;
hits := [];
simple_hits := [];

print "M_1(8,4) targeted [4,16] search";
print "height", height, "parameters", #params,
      "max_pr_tests", max_pr_tests, "max_hits", max_hits,
      "progress_interval", progress_interval;

for R in params do
    for w in params do
        if R eq 0 or w in {Q!-1, Q!0, Q!1} then
            continue;
        end if;
        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "progress", checked,
                  "square416", square416,
                  "cover", cover,
                  "smooth", smooth,
                  "tangent", tangent_points,
                  "tx_halves", tx_halves,
                  "pr_tests", pr_tests,
                  "pr_halves", pr_halves,
                  "torsion_tests", torsion_tests,
                  "hits", #hits,
                  "simple_hits", #simple_hits;
        end if;

        if not Square416Possible(R,w) then
            continue;
        end if;
        square416 +:= 1;

        if not M14CoverPossible(R,w) then
            continue;
        end if;
        cover +:= 1;

        f, t, A, B := FamilyPolynomial(R,w);
        if Degree(f) ne 5 or Discriminant(f) eq 0 then
            continue;
        end if;
        smooth +:= 1;

        candidates := TangentCandidates(f, R, w);
        if #candidates eq 0 then
            continue;
        end if;
        tangent_points +:= 1;

        fI, L := IntegralModelPolynomial(f);
        C := HyperellipticCurve(fI);
        J := Jacobian(C);
        Tx := J![x, Q!0];
        okTx, Hx := IsDivisibleBy(Tx, 2);
        if not okTx then
            print "WARNING tangent equations did not verify", R, w, t, candidates;
            continue;
        end if;
        tx_halves +:= 1;

        YR := YRValue(R,w);
        if YR^2 ne Evaluate(f, -R) then
            print "WARNING P_R formula failed", R, w;
            continue;
        end if;
        PR := J![x + R, Q!(L*YR)];

        pr_tests +:= 1;
        okPR, HPR := IsDivisibleBy(PR, 2);
        if not okPR then
            if pr_tests ge max_pr_tests then
                break R;
            end if;
            continue;
        end if;
        pr_halves +:= 1;

        torsion_tests +:= 1;
        G, phi := TorsionSubgroup(J);
        invs := Invariants(G);
        simple, pcert, Lp := IrreducibleFrobeniusCertificate(fI);
        is416 := Has416(invs);

        print "PR_HALF", "R", R, "w", w, "t", t,
              "torsion", invs, "has416", is416,
              "simple", simple, "pcert", pcert,
              "Lp", Lp;
        print "  fI", fI;
        print "  half_order_Tx", Order(Hx), "half_order_PR", Order(HPR);
        Append(~hits, <R,w,t,invs,simple,pcert,Lp,fI>);
        if is416 and simple then
            Append(~simple_hits, <R,w,t,invs,pcert,Lp,fI>);
            print "SIMPLE_416_HIT", "R", R, "w", w, "t", t,
                  "torsion", invs, "pcert", pcert;
            if #simple_hits ge max_hits then
                break R;
            end if;
        end if;

        if pr_tests ge max_pr_tests then
            break R;
        end if;
    end for;
end for;

print "DONE targeted [4,16] search";
print "height", height;
print "checked", checked;
print "square416", square416;
print "cover", cover;
print "smooth", smooth;
print "tangent_points", tangent_points;
print "tx_halves", tx_halves;
print "pr_tests", pr_tests;
print "pr_halves", pr_halves;
print "torsion_tests", torsion_tests;
print "hits", #hits;
print "simple_hits", #simple_hits;

quit;
