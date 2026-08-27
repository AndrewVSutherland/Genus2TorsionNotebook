//////////////////////////////////////////////////////////////////////
//  Finite-field sieve for the route: some order-6 class g+T is
//  divisible by 2 in J.
//
//  This tests every good parameter point in P^2(F_p).  Over finite
//  fields it uses AbelianGroup(J) and solves divisibility by 2 in the
//  resulting finite abelian group.  If no good F_p point works, any
//  rational example must reduce to the boundary modulo p.
//
//  Typical run from torsion_jac:
//      magma code/m2226_order6_halving_finite_field_sieve.m
//////////////////////////////////////////////////////////////////////

prime_list := [3,5,7,11,13,17,19,23,29,31,37,41,43];

function ProjectivePointsP2(F)
    pts := [];
    for s in F do
        for m in F do
            for n in F do
                if s eq 0 and m eq 0 and n eq 0 then
                    continue;
                end if;
                if s ne 0 and s ne 1 then
                    continue;
                end if;
                if s eq 0 and m ne 0 and m ne 1 then
                    continue;
                end if;
                if s eq 0 and m eq 0 and n ne 1 then
                    continue;
                end if;
                Append(~pts, [s,m,n]);
            end for;
        end for;
    end for;
    return pts;
end function;

function CurveData(F, vals)
    P<x> := PolynomialRing(F);
    s := vals[1];
    m := vals[2];
    n := vals[3];
    A := [ F!1, F!1, F!1, F!2, F!2 ];
    B := [
        2*s^2 - s*n,
        2*s^2 + s*m - 2*s*n - m*n,
        2*s^2 + s*m - s*n - m*n,
        -m*n,
        4*s^2 - 4*s*n - m*n
    ];
    if &or [ b eq 0 : b in B ] then
        return false, _, _, _, _;
    end if;
    betas := [ -A[i]/B[i] : i in [1..5] ];
    if #(Set(betas)) ne 5 then
        return false, _, _, _, _;
    end if;
    f := &*[ A[i] + B[i]*x : i in [1..5] ];
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false, _, _, _, _;
    end if;
    return true, f, betas, P, x;
end function;

function IsDivisibleBy2Finite(J, D)
    G, phi := AbelianGroup(J);
    a := D @@ phi;
    coords := Eltseq(a);
    invs := Invariants(G);
    for i in [1..#coords] do
        if GCD(2, invs[i]) ne 1 and coords[i] mod 2 ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function TwoTorsionClasses(J, betas, x)
    F := BaseRing(Parent(x));
    out := [ <"0", J!0> ];
    for i in [1..5] do
        Append(~out, <Sprintf("W%o-inf", i), J![x-betas[i], F!0]>);
    end for;
    for i in [1..4] do
        for j in [i+1..5] do
            Append(~out, <Sprintf("W%o+W%o", i, j), J![(x-betas[i])*(x-betas[j]), F!0]>);
        end for;
    end for;
    return out;
end function;

for p in prime_list do
    F := GF(p);
    pts := ProjectivePointsP2(F);
    good := 0;
    order6 := 0;
    divs := [];
    print "prime", p, "P2 points", #pts;

    for vals in pts do
        ok, f, betas, P, x := CurveData(F, vals);
        if not ok then
            continue;
        end if;
        good +:= 1;
        C := HyperellipticCurve(f);
        J := Jacobian(C);
        g := J![x, F!2];
        for data in TwoTorsionClasses(J, betas, x) do
            label := data[1];
            D := g + data[2];
            if Order(D) ne 6 then
                continue;
            end if;
            order6 +:= 1;
            if IsDivisibleBy2Finite(J, D) then
                Append(~divs, <vals, label>);
                print "  divisible", vals, label;
            end if;
        end for;
    end for;

    print "  good", good, "order6 classes", order6, "divisible", #divs;
    if #divs eq 0 and good gt 0 then
        print "  NO GOOD DIVISIBLE POINTS MOD", p;
        break;
    end if;
end for;
quit;
