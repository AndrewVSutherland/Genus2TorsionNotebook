//////////////////////////////////////////////////////////////////////
//  Search for [6,6] on the contact-6 chart by an algebraic route.
//
//  The contact-6 family is
//
//      h6 = 1 + a*x + b*x^2 + x^3,
//      f  = h6^2 - (x - 1)^6.
//
//  It already has the rational Weierstrass point x=0.  To get genuine
//  [6,6] from an independent 3-torsion class, we also need an independent
//  rational 2-torsion class.  The clean condition is to force another
//  rational root r != 0,1:
//
//      h6(r) = eps*(r - 1)^3,   eps in {+1,-1}.
//
//  This gives
//
//      a = (eps*(r-1)^3 - 1 - r^3 - b*r^2)/r.
//
//  We then impose the cubic-contact 3-torsion equations.  For fixed
//  (eps,r,L,U,v), the first contact equation is quadratic in b, so the
//  search solves for b algebraically and checks the remaining equations.
//
//  Typical runs:
//      magma -b height:=4 code/contact6_m36_extra_root_search.m
//      magma -b height:=6 progress_interval:=1000000 \
//          code/contact6_m36_extra_root_search.m > data/contact6_m36_extra_root_h6.txt
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned height then
    height := 5;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;
if not assigned progress_interval then
    progress_interval := 200000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;
if not assigned max_contact_print then
    max_contact_print := 20;
elif Type(max_contact_print) eq MonStgElt then
    max_contact_print := StringToInteger(max_contact_print);
end if;
if not assigned fixed_eps then
    fixed_eps := 0;
elif Type(fixed_eps) eq MonStgElt then
    fixed_eps := StringToInteger(fixed_eps);
end if;
if not assigned fixed_r_num then
    fixed_r_num := 0;
elif Type(fixed_r_num) eq MonStgElt then
    fixed_r_num := StringToInteger(fixed_r_num);
end if;
if not assigned fixed_r_den then
    fixed_r_den := 1;
elif Type(fixed_r_den) eq MonStgElt then
    fixed_r_den := StringToInteger(fixed_r_den);
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);
PB<z> := PolynomialRing(Q);

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            q := Q!num/den;
            key := Sprint(q);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, q);
            end if;
        end for;
    end for;
    return vals;
end function;

function Contact6Polynomial(a, b)
    h := 1 + a*x + b*x^2 + x^3;
    f := h^2 - (x-1)^6;
    return f, h;
end function;

function ExtraRootA(eps, r, b)
    return (eps*(r-1)^3 - 1 - r^3 - b*r^2)/r;
end function;

function GoodPolynomial(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function ThreeRankFromInvariants(invs)
    return #[n : n in invs | (Z!n) mod 3 eq 0];
end function;

function Has66(invs)
    return #[n : n in invs | (Z!n) mod 6 eq 0] ge 2;
end function;

function SimpleCertificate(f)
    C := HyperellipticCurve(f);
    for p in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97] do
        try
            fp := ChangeRing(f, GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
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

function CubicContactData(eps, r, b, L, U, v)
    a := ExtraRootA(eps, r, b);
    f, h6 := Contact6Polynomial(a, b);
    c1 := Coefficient(f, 1);
    c2 := Coefficient(f, 2);
    c3 := Coefficient(f, 3);
    c4 := Coefficient(f, 4);
    c5 := Coefficient(f, 5);

    B := c5*L^2 + 3*U;
    Delta := 4*c4*L^2 + 12*(U^2 + v^2) - B^2;
    F3 := B*Delta + 16*v^3 - 8*c3*L^2 - 8*U^3 - 48*U*v^2;
    F2 := Delta^2 + 64*B*v^3 - 64*c2*L^2
          - 192*(U^2*v^2 + v^4);
    F1 := Delta*v^3 - 4*c1*L^2 - 12*U*v^4;
    q := x^2 + U*x + v^2;
    h3 := (1/L)*x^3 + (B/(2*L))*x^2 + (Delta/(8*L))*x + v^3/L;
    return a, F1, F2, F3, q, h3, f, h6;
end function;

function QuadraticRootsFromValues(y0, y1, y2)
    c2 := (y2 - 2*y1 + y0)/2;
    c1 := y1 - y0 - c2;
    c0 := y0;
    poly := PB![c0, c1, c2];
    roots := [];
    if poly eq 0 then
        return roots, true;
    end if;
    for rt in Roots(poly) do
        Append(~roots, rt[1]);
    end for;
    return roots, false;
end function;

function ExactCheck(eps, r, b, L, U, v)
    a, F1, F2, F3, q, h3, f, h6 := CubicContactData(eps, r, b, L, U, v);
    if F1 ne 0 or F2 ne 0 or F3 ne 0 then
        return false, [], 0, 0, 0, false, 0, P!0, a, f;
    end if;
    if not GoodPolynomial(f) or Evaluate(h6, Q!1) eq 0 then
        return false, [], 0, 0, 0, false, 0, P!0, a, f;
    end if;
    if r eq 0 or r eq 1 or Evaluate(f, r) ne 0 then
        return false, [], 0, 0, 0, false, 0, P!0, a, f;
    end if;
    if Degree(GCD(q, f)) gt 0 or Discriminant(q) eq 0 then
        return false, [], 0, 0, 0, false, 0, P!0, a, f;
    end if;

    C := HyperellipticCurve(f);
    J := Jacobian(C);
    D := J![x-1, Evaluate(h6, Q!1)];
    E := J![q, h3 mod q];
    T0 := J![x, 0];
    Tr := J![x-r, 0];
    ordD := Order(D);
    ordE := Order(E);
    G, phi := TorsionSubgroup(J);
    simple, pcert, Lp := SimpleCertificate(f);

    // Independence of the 2-parts is equivalent to T0 != Tr.
    independent2 := T0 ne Tr and T0 ne J!0 and Tr ne J!0;
    return true, Invariants(G), ordD, ordE, independent2, simple, pcert, Lp, a, f;
end function;

params := RationalParametersOfHeight(height);
if fixed_eps eq 1 then
    eps_values := [Q!1];
elif fixed_eps eq -1 then
    eps_values := [-Q!1];
else
    eps_values := [Q!1, -Q!1];
end if;
if fixed_r_num ne 0 then
    r_values := [Q!fixed_r_num / Q!fixed_r_den];
else
    r_values := params;
end if;
checked := 0;
nonboundary := 0;
quadratic_roots := 0;
zero_f1_polys := 0;
contact_hits := 0;
exact_hits66 := 0;
hits := [];

print "Contact-6 extra-root algebraic search for [6,6]";
print "height", height, "parameter_count", #params,
      "eps_values", eps_values, "r_values", r_values,
      "max_hits", max_hits, "progress_interval", progress_interval;

for eps in eps_values do
    for r in r_values do
        if r eq 0 or r eq 1 then
            continue;
        end if;
        for L in params do
            if L eq 0 then
                continue;
            end if;
            for U in params do
                for v in params do
                    checked +:= 1;
                    if progress_interval gt 0 and checked mod progress_interval eq 0 then
                        print "progress", checked,
                              "nonboundary", nonboundary,
                              "quadratic_roots", quadratic_roots,
                              "contact_hits", contact_hits,
                              "hits66", exact_hits66;
                    end if;
                    if v eq 0 or U^2 - 4*v^2 eq 0 then
                        continue;
                    end if;
                    nonboundary +:= 1;

                    a0, y0, y20, y30, q0, h30, f0, hh0 :=
                        CubicContactData(eps, r, Q!0, L, U, v);
                    a1, y1, y21, y31, q1, h31, f1, hh1 :=
                        CubicContactData(eps, r, Q!1, L, U, v);
                    a2, y2, y22, y32, q2, h32, f2, hh2 :=
                        CubicContactData(eps, r, Q!2, L, U, v);

                    broots, zero_poly := QuadraticRootsFromValues(y0, y1, y2);
                    if zero_poly then
                        zero_f1_polys +:= 1;
                        continue;
                    end if;
                    for b in broots do
                        quadratic_roots +:= 1;
                        a, F1, F2, F3, q, h3, f, h6 :=
                            CubicContactData(eps, r, b, L, U, v);
                        if F1 ne 0 or F2 ne 0 or F3 ne 0 then
                            continue;
                        end if;
                        contact_hits +:= 1;
                        if contact_hits le max_contact_print then
                            print "CONTACT", "eps", eps, "r", r,
                                  "a", a, "b", b, "L", L, "U", U, "v", v;
                            print " f", f;
                            print " q", q;
                            print " h3modq", h3 mod q;
                        end if;

                        ok, invs, ordD, ordE, independent2, simple, pcert, Lp, a0, f0 :=
                            ExactCheck(eps, r, b, L, U, v);
                        if not ok then
                            continue;
                        end if;
                        if ordD eq 6 and ordE eq 3 and independent2 and Has66(invs) then
                            exact_hits66 +:= 1;
                            Append(~hits, <eps,r,a0,b,L,U,v,invs,simple,pcert,f0>);
                            print "HIT66", "eps", eps, "r", r,
                                  "a", a0, "b", b, "L", L, "U", U, "v", v,
                                  "invs", invs, "simple", simple, "pcert", pcert;
                            print " f", f0;
                            if #hits ge max_hits then
                                break eps;
                            end if;
                        end if;
                    end for;
                end for;
            end for;
        end for;
    end for;
end for;

print "Done";
print "checked", checked;
print "nonboundary", nonboundary;
print "quadratic_roots", quadratic_roots;
print "zero_f1_polys", zero_f1_polys;
print "contact_hits", contact_hits;
print "exact_hits66", exact_hits66;
print "hits", #hits;
for H in hits do
    print "H", H;
end for;

quit;
