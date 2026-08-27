//////////////////////////////////////////////////////////////////////
//  Sample search on the A/M(2,4,4) elliptic-fiber-product model.
//
//  E2: y^2 = X(X+t^2)(X+4s+t^2)
//  E3: y^2 = X(X^2 - 4(t^2+2s)X + 16s^2)
//
//  A point R in E3 and P=(x1,y1) in E2 gives Q=phi(R)+P in E2.
//  With Q=(x2,y2), set u=y1/(2*x1), v=y2/(2*x2).  The genus-2 curve is
//
//      y^2 = x(x+u^2)(x+v^2)(x^2 + (t^2+2s)x + s^2).
//
//  This script samples small rational s,t and small rational points on E2,E3,
//  then computes exact rational torsion.  Hits containing [2,4,8] are printed.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned height then
    height := 6;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

if not assigned point_bound then
    point_bound := 80;
elif Type(point_bound) eq MonStgElt then
    point_bound := StringToInteger(point_bound);
end if;

if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;

if not assigned max_curves then
    max_curves := 0;
elif Type(max_curves) eq MonStgElt then
    max_curves := StringToInteger(max_curves);
end if;

if not assigned targeted_only then
    targeted_only := false;
elif Type(targeted_only) eq MonStgElt then
    targeted_only := targeted_only in {"true", "True", "1", "yes"};
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

function GoodHyperellipticPolynomial(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function Curves(s,t)
    A := t^2 + 2*s;
    E2 := EllipticCurve([Q!0, Q!(2*A), Q!0, Q!(A^2 - 4*s^2), Q!0]);
    E3 := EllipticCurve([Q!0, Q!(-4*A), Q!0, Q!(16*s^2), Q!0]);
    return E2, E3;
end function;

function DualPhiE3ToE2(R, s, E2)
    E3 := Curve(R);
    if R eq E3!0 then
        return E2!0;
    end if;
    X := R[1];
    Y := R[2];
    if X eq 0 then
        return E2!0;
    end if;
    return E2![Y^2/(4*X^2), Y*(16*s^2 - X^2)/(8*X^2)];
end function;

function Has248(invs)
    vals := Sort([Valuation(n, 2) : n in invs]);
    vals := Reverse(vals);
    return #vals ge 3 and vals[1] ge 3 and vals[2] ge 2 and vals[3] ge 1;
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

function CurvePolynomial(s,t,P1,R)
    E2 := Curve(P1);
    Q2 := P1 + DualPhiE3ToE2(R, s, E2);
    if P1 eq E2!0 or Q2 eq E2!0 then
        return false, P!0, _, _, _;
    end if;
    if P1[1] eq 0 or Q2[1] eq 0 then
        return false, P!0, _, _, _;
    end if;
    u := P1[2]/(2*P1[1]);
    v := Q2[2]/(2*Q2[1]);
    if u eq 0 or v eq 0 or u^2 eq v^2 then
        return false, P!0, _, _, _;
    end if;
    f := x*(x+u^2)*(x+v^2)*(x^2 + (t^2+2*s)*x + s^2);
    if not GoodHyperellipticPolynomial(f) then
        return false, P!0, _, _, _;
    end if;
    return true, f, u, v, Q2;
end function;

params := RationalParametersOfHeight(height);
hits := [];
tested_curves := 0;
torsion_counts := AssociativeArray();
seen_curves := {};

print "M(2,4,4) -> [2,4,8] sample search";
print "height", height, "parameter_count", #params,
      "point_bound", point_bound, "max_hits", max_hits,
      "max_curves", max_curves,
      "targeted_only", targeted_only;

for s in params do
    if #hits ge max_hits then
        break;
    end if;
    for t in params do
        if #hits ge max_hits then
            break;
        end if;
        if s eq 0 or t eq 0 or t^2 + 4*s eq 0 then
            continue;
        end if;
        E2, E3 := Curves(s,t);
        if Discriminant(E2) eq 0 or Discriminant(E3) eq 0 then
            continue;
        end if;
        pts2 := Points(E2 : Bound := point_bound);
        pts3 := Points(E3 : Bound := point_bound);
        if #pts2 le 4 or #pts3 le 2 then
            continue;
        end if;
        print "fiber", "s", s, "t", t, "#E2pts", #pts2, "#E3pts", #pts3;

        for P1 in pts2 do
            if #hits ge max_hits then
                break;
            end if;
            for R in pts3 do
                if #hits ge max_hits then
                    break;
                end if;
                ok, f, u, v, Q2 := CurvePolynomial(s,t,P1,R);
                if not ok then
                    continue;
                end if;
                fI, L := IntegralModelPolynomial(f);
                key_curve := Sprint(fI);
                if key_curve in seen_curves then
                    continue;
                end if;
                Include(~seen_curves, key_curve);
                tested_curves +:= 1;
                if max_curves gt 0 and tested_curves gt max_curves then
                    break s;
                end if;

                C := HyperellipticCurve(fI);
                J := Jacobian(C);

                // Target the order-4 class above the irreducible quadratic
                // q=x^2+(t^2+2s)x+s^2.  The M(2,4,4) construction should make
                // Tq=[q,0] divisible by 2; [2,4,8] on this branch means that
                // a chosen half Hq is again divisible by 2 over Q.
                qpoly := x^2 + (t^2 + 2*s)*x + s^2;
                Tq := J![qpoly, Q!0];
                ok4, Hq := IsDivisibleBy(Tq, 2);
                ok8 := false;
                if ok4 then
                    ok8, Q8 := IsDivisibleBy(Hq, 2);
                end if;
                if targeted_only and not ok8 then
                    continue;
                end if;

                G, phi := TorsionSubgroup(J);
                invs := Invariants(G);
                key := Sprint(invs);
                if IsDefined(torsion_counts, key) then
                    torsion_counts[key] +:= 1;
                else
                    torsion_counts[key] := 1;
                end if;

                if ok8 or Has248(invs) then
                    simple, pcert, Lp := IrreducibleFrobeniusCertificate(fI);
                    Append(~hits, <s,t,u,v,P1,R,Q2,invs,simple,pcert,Lp,fI>);
                    print "HIT", "s", s, "t", t,
                          "u", u, "v", v,
                          "quad_half_divisible", ok8,
                          "torsion", invs,
                          "order", TorsionOrder(invs),
                          "exponent", Exponent(invs),
                          "simple", simple, "pcert", pcert;
                    print "  P", P1;
                    print "  R", R;
                    print "  Q", Q2;
                    print "  fI", fI;
                end if;
            end for;
        end for;
    end for;
end for;

print "DONE";
print "tested_curves", tested_curves, "hits", #hits;
print "TORSION_COUNTS";
for key in Sort([k : k in Keys(torsion_counts)]) do
    print " ", key, torsion_counts[key];
end for;

quit;
