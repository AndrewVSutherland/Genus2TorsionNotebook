//////////////////////////////////////////////////////////////////////
//  Finite-field diagnostic for the contact-7 first-halving surface.
//
//  This enumerates the surface from contact7_halving_surface_search.m
//  over F_p and checks whether the explicit order-4 half H4 of the
//  rational Weierstrass class is divisible by 2 in J(F_p).  If no open
//  point over F_p has H4 divisible by 2, then rational [56] examples on
//  this open surface must reduce to the boundary at p.
//
//  Typical run:
//      magma code/contact7_halving_surface_finite.m
//////////////////////////////////////////////////////////////////////

if not assigned max_print then
    max_print := 8;
elif Type(max_print) eq MonStgElt then
    max_print := StringToInteger(max_print);
end if;

prime_list := [3,5,7,11,13,17,19,23,29,31];

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

for p in prime_list do
    F := GF(p);
    P<x> := PolynomialRing(F);
    eps_vals := [F!-1, F!1];

    checked := 0;
    root_good := 0;
    surface := 0;
    exact_h4 := 0;
    h4_div2 := 0;
    surface_samples := [];
    samples := [];

    for s in F do
        for u in F do
            for z in F do
                for eps in eps_vals do
                    checked +:= 1;
                    if s eq 0 or s^2 eq 1 or u eq 0 then
                        continue;
                    end if;
                    S := eps*s;
                    if S eq 0 or S + 1 eq 0 then
                        continue;
                    end if;

                    A0 := 3*S^4 + 9*S^3 + 11*S^2 + 9*S + 3;
                    c00 := S^7*A0/(S + 1)^3;
                    lam := 2*S^7;
                    if lam eq 0 then
                        continue;
                    end if;
                    b := (z^2 - c00)/lam;

                    r := 1 - s^2;
                    if r eq 0 then
                        continue;
                    end if;
                    a := (eps*s^7 - 1 + (F!7/F!2)*r - b*r^3)/r^2;
                    h := 1 - (F!7/F!2)*x + a*x^2 + b*x^3;
                    f := ExactQuotient(h^2 + (x - 1)^7, x^2);
                    if Degree(f) ne 5 or Discriminant(f) eq 0 then
                        continue;
                    end if;
                    if Evaluate(f, r) ne 0 then
                        continue;
                    end if;
                    root_good +:= 1;

                    g := ExactQuotient(Evaluate(f, x + r), x);
                    c3 := Coefficient(g, 3);
                    c2 := Coefficient(g, 2);
                    c1 := Coefficient(g, 1);
                    c0 := Coefficient(g, 0);
                    if c0 ne z^2 then
                        continue;
                    end if;
                    v := (u^2 - c3)/2;
                    w := (v^2 + 2*z - c2)/(2*u);
                    if c1 ne w^2 - 2*v*z then
                        continue;
                    end if;
                    surface +:= 1;

                    C := HyperellipticCurve(f);
                    J := Jacobian(C);
                    D2 := J![x - r, F!0];
                    D7 := J![x - 1, Evaluate(h, F!1)];
                    X := x - r;
                    Qpoly := X^2 - v*X + z;
                    alpha := (u*v - w)*X - u*z;
                    H4 := J![Qpoly, alpha];
                    if not (Order(D2) eq 2 and Order(H4) eq 4 and 2*H4 eq D2) then
                        continue;
                    end if;
                    exact_h4 +:= 1;
                    if #surface_samples lt max_print then
                        Append(~surface_samples, <s,u,z,eps,a,b,r,#J,Order(D7)>);
                    end if;

                    if IsDivisibleBy2Finite(J, H4) then
                        h4_div2 +:= 1;
                        if #samples lt max_print then
                            Append(~samples, <s,u,z,eps,a,b,r,#J,Order(D7)>);
                        end if;
                    end if;
                end for;
            end for;
        end for;
    end for;

    print "p", p,
          "checked", checked,
          "root_good", root_good,
          "surface", surface,
          "exact_h4", exact_h4,
          "h4_divisible_by_2", h4_div2;
    if #samples gt 0 then
        print "  samples", samples;
    else
        print "  NO OPEN H4 HALVES";
        if #surface_samples gt 0 then
            print "  surface_samples", surface_samples;
        end if;
    end if;
end for;

quit;
