//////////////////////////////////////////////////////////////////////
//  Exact rational search for [8,8] on the reduced square subcover.
//
//  The reduced [8,8] conditions force V to be a square on the first
//  [4,8] cover.  Write w=s^2 and use the first-cover sign
//      V = R^2*s^2.
//
//  For each rational R,s of bounded height this script:
//    1. builds the M_1(8,4) curve with w=s^2,
//    2. solves the first-cover tangent equations with V=R^2*s^2,
//    3. solves the reduced second-halving equations by resultants,
//    4. verifies H_x is divisible by 2 in the exact Jacobian,
//    5. computes exact rational torsion.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

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

if not assigned use_filter then
    use_filter := true;
elif Type(use_filter) eq MonStgElt then
    use_filter := use_filter in {"true", "True", "1", "yes"};
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

function BoundaryBad(R, w)
    if R eq 0 or w eq 0 or w eq 1 or w eq -1 then
        return true;
    end if;
    if R eq 1 or R eq -1 or R eq w or R eq -w then
        return true;
    end if;
    if R*w - 3*R + 3*w - 1 eq 0 then
        return true;
    end if;
    if R*w + 3*R + 3*w + 1 eq 0 then
        return true;
    end if;
    if 2*R^2 - R*w^2 + R - 2*w^2 eq 0 then
        return true;
    end if;
    if R^4 - 2*R^3 + R^2*w^2 - R^2 + 2*R*w^2 - w^2 eq 0 then
        return true;
    end if;
    return false;
end function;

function ResidueInt(q, p)
    den := Z!Denominator(q);
    if den mod p eq 0 then
        return false, 0;
    end if;
    num := Z!Numerator(q);
    return true, ((num mod p) * InverseMod(den mod p, p)) mod p;
end function;

function PassesReducedResidueFilter(R, s)
    allowed7 := {<3,3>, <3,4>, <5,2>, <5,5>};
    allowed11 := {<2,2>, <2,9>, <5,5>, <5,6>, <6,5>, <6,6>,
                  <7,4>, <7,7>, <8,3>, <8,8>, <9,2>, <9,9>};
    allowed13 := {<3,2>, <3,11>, <6,6>, <6,7>,
                  <9,6>, <9,7>, <11,2>, <11,11>};
    filters := [<7, allowed7>, <11, allowed11>, <13, allowed13>];

    for filt in filters do
        p := filt[1];
        allowed := filt[2];
        okR, rp := ResidueInt(R, p);
        okS, sp := ResidueInt(s, p);
        if not (okR and okS) then
            continue;
        end if;
        if <rp,sp> notin allowed then
            return false;
        end if;
    end for;
    return true;
end function;

function FirstCoverSolutions(f, R, s)
    w := s^2;
    V := R^2*s^2;
    h := ExactQuotient(f, x);
    c1 := Coefficient(h, 1);
    c2 := Coefficient(h, 2);
    c3 := Coefficient(h, 3);
    c4 := Coefficient(h, 4);

    PR<Uvar> := PolynomialRing(Q);
    F := 4*(c3 - 2*c4*Uvar)*(c1 - 2*c4*Uvar*V)
         - (c2 - c4*(Uvar^2 + 2*V))^2;

    sols := [];
    for rt in Roots(F) do
        U0 := rt[1];
        M2 := c3 - 2*c4*U0;
        N2 := c1 - 2*c4*U0*V;
        for M0 in SquareRootsQ(M2) do
            for N0 in SquareRootsQ(N2) do
                a0 := x^2 + U0*x + V;
                identity := h - x*(M0*x+N0)^2 - c4*a0^2;
                if identity eq 0 then
                    Append(~sols, <U0, V, M0, N0, c4>);
                end if;
            end for;
        end for;
    end for;
    return sols;
end function;

function SpecializeTauPoly(poly, z0)
    QT<T> := PolynomialRing(Q);
    out := QT!0;
    for i in [0..Degree(poly)] do
        out +:= Q!Evaluate(Coefficient(poly, i), z0)*T^i;
    end for;
    return out;
end function;

function ReducedSecondSolutions(R, s, U0, V, M0, N0, c4)
    out := [];
    S<z> := PolynomialRing(Q);
    T<tau> := PolynomialRing(S);
    PX<X> := PolynomialRing(T);

    for eta in [Q!1, Q!-1] do
        a0 := X^2 + (T!U0)*X + T!V;
        ell0 := X*((T!M0)*X + T!N0);
        a := T!(U0/2) + (T!M0)*(T!z) + tau - T!(c4/2)*(T!z)^2;
        b := T!(eta*R*s)*tau;
        q := X^2 + a*X + b;
        second := (T!c4)*(T!z)^2*X*a0
                  - 2*(T!z)*ell0*(X+tau)
                  - a0*(X+tau)^2 + q^2;
        E1 := Coefficient(second, 1);
        E2 := Coefficient(second, 2);
        res := Resultant(E1, E2);
        if res eq 0 then
            print "WARNING zero resultant", "R", R, "s", s,
                  "U", U0, "M", M0, "N", N0, "eta", eta;
            continue;
        end if;
        res := S!res;
        for zrt in Roots(res) do
            z0 := zrt[1];
            if z0 eq 0 then
                continue;
            end if;
            E1z := SpecializeTauPoly(E1, z0);
            E2z := SpecializeTauPoly(E2, z0);
            tau_roots := [];
            if E1z eq 0 and E2z eq 0 then
                print "WARNING line of tau solutions", "R", R, "s", s,
                      "U", U0, "M", M0, "N", N0, "eta", eta, "z", z0;
                continue;
            elif E1z eq 0 then
                tau_roots := Roots(E2z);
            elif E2z eq 0 then
                tau_roots := Roots(E1z);
            else
                tau_roots := Roots(GCD(E1z, E2z));
            end if;
            for trt in tau_roots do
                tau0 := trt[1];
                a_val := U0/2 + M0*z0 + tau0 - (c4/2)*z0^2;
                b_val := eta*R*s*tau0;
                Append(~out, <z0, tau0, eta, a_val, b_val>);
            end for;
        end for;
    end for;
    return out;
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

function Has88(invs)
    vals := Sort([Valuation(n, 2) : n in invs]);
    vals := Reverse(vals);
    return #vals ge 2 and vals[1] ge 3 and vals[2] ge 3;
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
residue_filtered := 0;
open_good := 0;
first_bases := 0;
first_solutions := 0;
second_solutions := 0;
exact_verified := 0;
torsion_tests := 0;
hits := [];
torsion_counts := AssociativeArray();

print "Reduced [8,8] exact search";
print "height", height, "params", #params, "pairs", #params^2,
      "max_hits", max_hits, "progress_interval", progress_interval,
      "use_filter", use_filter;

for R in params do
    for s in params do
        if #hits ge max_hits then
            break R;
        end if;
        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "progress", checked, "residue_filtered", residue_filtered,
                  "open_good", open_good,
                  "first_bases", first_bases,
                  "first_solutions", first_solutions,
                  "second_solutions", second_solutions,
                  "verified", exact_verified,
                  "hits", #hits;
        end if;

        w := s^2;
        if use_filter and not PassesReducedResidueFilter(R, s) then
            residue_filtered +:= 1;
            continue;
        end if;
        if BoundaryBad(R, w) then
            continue;
        end if;
        f, t, A, B := FamilyPolynomial(R, w);
        if not GoodHyperellipticPolynomial(f) then
            continue;
        end if;
        open_good +:= 1;

        first := FirstCoverSolutions(f, R, s);
        if #first eq 0 then
            continue;
        end if;
        first_bases +:= 1;
        first_solutions +:= #first;

        fI, L := IntegralModelPolynomial(f);
        C := HyperellipticCurve(fI);
        J := Jacobian(C);
        Tx := J![x, Q!0];

        for fs in first do
            U0 := fs[1];
            V := fs[2];
            M0 := fs[3];
            N0 := fs[4];
            c4 := fs[5];
            seconds := ReducedSecondSolutions(R, s, U0, V, M0, N0, c4);
            second_solutions +:= #seconds;
            if #seconds eq 0 then
                continue;
            end if;

            a0 := x^2 + U0*x + V;
            v0 := (M0*U0 - N0)*x + M0*V;
            Hx := J![a0, Q!(L*v0)];
            ok_half, Q16 := IsDivisibleBy(Hx, 2);
            if not ok_half then
                print "WARNING reduced equations did not verify exact halving",
                      "R", R, "s", s, "U", U0, "M", M0, "N", N0,
                      "seconds", seconds;
                continue;
            end if;
            exact_verified +:= 1;
            if 2*Hx ne Tx then
                print "WARNING Hx is not half of Tx", "R", R, "s", s;
                continue;
            end if;

            torsion_tests +:= 1;
            G, phi := TorsionSubgroup(J);
            invs := Invariants(G);
            key := Sprint(invs);
            if IsDefined(torsion_counts, key) then
                torsion_counts[key] +:= 1;
            else
                torsion_counts[key] := 1;
            end if;

            if Has88(invs) then
                simple, pcert, Lp := IrreducibleFrobeniusCertificate(fI);
                Append(~hits, <R,s,w,t,invs,simple,pcert,Lp,fI,fs,seconds[1]>);
                print "HIT", "R", R, "s", s, "w", w, "t", t,
                      "torsion", invs, "order", TorsionOrder(invs),
                      "exponent", Exponent(invs),
                      "simple", simple, "pcert", pcert;
                print "  first", fs;
                print "  second", seconds[1];
                print "  fI", fI;
                if #hits ge max_hits then
                    break R;
                end if;
            else
                print "NON_88_VERIFIED", "R", R, "s", s, "w", w,
                      "torsion", invs;
            end if;
        end for;
    end for;
end for;

print "DONE height", height;
print "checked", checked, "residue_filtered", residue_filtered,
      "open_good", open_good,
      "first_bases", first_bases, "first_solutions", first_solutions,
      "second_solutions", second_solutions,
      "exact_verified", exact_verified,
      "torsion_tests", torsion_tests,
      "hits", #hits;
print "TORSION_COUNTS";
for key in Sort([k : k in Keys(torsion_counts)]) do
    print " ", key, torsion_counts[key];
end for;

quit;
