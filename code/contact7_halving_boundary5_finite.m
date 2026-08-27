//////////////////////////////////////////////////////////////////////
//  F_5 boundary diagnostic for contact-7 first-halving surface.
//
//  The earlier finite script solved the first-halving equations by dividing
//  by u, so it deliberately removed the u=0 chart.  This script enumerates
//  the unsolved symmetric equations
//
//      c3 = u^2 - 2v,
//      c2 = v^2 - 2uw + 2z,
//      c1 = w^2 - 2vz,
//      c0 = z^2,
//
//  over F_5 in the open contact-7 rational-root chart for s.  It therefore
//  includes u=0 boundary points and checks whether the resulting order-4
//  class H4 is divisible by 2 in J(F_5).
//////////////////////////////////////////////////////////////////////

if not assigned p then
    p := 5;
elif Type(p) eq MonStgElt then
    p := StringToInteger(p);
end if;

if not assigned max_print then
    max_print := 30;
elif Type(max_print) eq MonStgElt then
    max_print := StringToInteger(max_print);
end if;

F := GF(p);
P<x> := PolynomialRing(F);

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

print "contact-7 first-halving unsolved finite boundary diagnostic";
print "p", p;

root_checked := 0;
root_good := 0;
surface := 0;
surface_u0 := 0;
exact_h4 := 0;
exact_h4_u0 := 0;
h4_div2 := 0;
h4_div2_u0 := 0;
samples := [];
samples_u0 := [];
div2_samples := [];
div2_samples_u0 := [];

for s in F do
    if s eq 0 or s^2 eq 1 then
        continue;
    end if;
    for eps in [F!-1, F!1] do
        r := 1 - s^2;
        if r eq 0 then
            continue;
        end if;
        for b in F do
            root_checked +:= 1;
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

            C := HyperellipticCurve(f);
            J := Jacobian(C);
            D2 := J![x - r, F!0];
            D7 := J![x - 1, Evaluate(h, F!1)];
            X := x - r;

            for u in F do
                for v in F do
                    if c3 ne u^2 - 2*v then
                        continue;
                    end if;
                    for w in F do
                        for z in F do
                            if c2 ne v^2 - 2*u*w + 2*z then
                                continue;
                            end if;
                            if c1 ne w^2 - 2*v*z then
                                continue;
                            end if;
                            if c0 ne z^2 then
                                continue;
                            end if;

                            surface +:= 1;
                            if u eq 0 then
                                surface_u0 +:= 1;
                            end if;
                            if #samples lt max_print then
                                Append(~samples, <s,u,v,w,z,eps,a,b,r,#J,Order(D7)>);
                            end if;
                            if u eq 0 and #samples_u0 lt max_print then
                                Append(~samples_u0, <s,u,v,w,z,eps,a,b,r,#J,Order(D7)>);
                            end if;

                            Qpoly := X^2 - v*X + z;
                            alpha := (u*v - w)*X - u*z;
                            H4 := J![Qpoly, alpha];
                            if not (Order(D2) eq 2 and Order(H4) eq 4 and 2*H4 eq D2) then
                                continue;
                            end if;
                            exact_h4 +:= 1;
                            if u eq 0 then
                                exact_h4_u0 +:= 1;
                            end if;

                            if IsDivisibleBy2Finite(J, H4) then
                                h4_div2 +:= 1;
                                if #div2_samples lt max_print then
                                    Append(~div2_samples, <s,u,v,w,z,eps,a,b,r,#J,Order(D7)>);
                                end if;
                                if u eq 0 then
                                    h4_div2_u0 +:= 1;
                                    if #div2_samples_u0 lt max_print then
                                        Append(~div2_samples_u0, <s,u,v,w,z,eps,a,b,r,#J,Order(D7)>);
                                    end if;
                                end if;
                            end if;
                        end for;
                    end for;
                end for;
            end for;
        end for;
    end for;
end for;

print "root_checked", root_checked;
print "root_good", root_good;
print "surface", surface;
print "surface_u0", surface_u0;
print "exact_h4", exact_h4;
print "exact_h4_u0", exact_h4_u0;
print "h4_divisible_by_2", h4_div2;
print "h4_divisible_by_2_u0", h4_div2_u0;
print "samples", samples;
print "samples_u0", samples_u0;
print "div2_samples", div2_samples;
print "div2_samples_u0", div2_samples_u0;

quit;
