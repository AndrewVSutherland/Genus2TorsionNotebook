//////////////////////////////////////////////////////////////////////
//  Targeted search for the [3,6] cover inside the contact-6 family.
//
//  This uses the eliminated cubic-contact equations from
//  code/contact6_m36_symbolic.m.  Instead of scanning (a,b), it scans
//  (b,L,U,v), solves F1=0 for a, and then checks F2=F3=0 exactly.
//
//  Typical run:
//      magma -b height:=4 code/contact6_m36_three_torsion_search.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned height then
    height := 4;
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
if not assigned max_contact_print then
    max_contact_print := 20;
elif Type(max_contact_print) eq MonStgElt then
    max_contact_print := StringToInteger(max_contact_print);
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

function GoodPolynomial(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function ThreeRankFromInvariants(invs)
    return #[n : n in invs | (Z!n) mod 3 eq 0];
end function;

function AnyDivisible(invs, d)
    for n in invs do
        if (Z!n) mod d eq 0 then
            return true;
        end if;
    end for;
    return false;
end function;

function Has36(invs)
    return ThreeRankFromInvariants(invs) ge 2 and AnyDivisible(invs, 6);
end function;

function Has66(invs)
    return #[n : n in invs | (Z!n) mod 6 eq 0] ge 2;
end function;

function CubicContactData(a, b, L, U, v)
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
    return F1, F2, F3, q, h3, f, h6;
end function;

function ExactCheck(a, b, L, U, v)
    F1, F2, F3, q, h3, f, h6 := CubicContactData(a, b, L, U, v);
    if F1 ne 0 or F2 ne 0 or F3 ne 0 then
        return false, [], 0, 0, 0;
    end if;
    if not GoodPolynomial(f) or Evaluate(h6, Q!1) eq 0 then
        return false, [], 0, 0, 0;
    end if;
    if Degree(GCD(q, f)) gt 0 or Discriminant(q) eq 0 then
        return false, [], 0, 0, 0;
    end if;

    C := HyperellipticCurve(f);
    J := Jacobian(C);
    D := J![x-1, Evaluate(h6, Q!1)];
    E := J![q, h3 mod q];
    ordD := Order(D);
    ordE := Order(E);
    G, phi := TorsionSubgroup(J);
    return true, Invariants(G), ordD, ordE, f;
end function;

params := RationalParametersOfHeight(height);
PA<z> := PolynomialRing(Q);

checked := 0;
nonboundary := 0;
linear_roots := 0;
contact_hits := 0;
exact_hits36 := 0;
exact_hits66 := 0;
special_f1 := 0;
hits := [];

print "Contact-6 targeted [3,6] cubic-contact search";
print "height", height, "parameter_count", #params, "max_hits", max_hits;

for b in params do
    for L in params do
        if L eq 0 then
            continue;
        end if;
        for U in params do
            for v in params do
                checked +:= 1;
                if progress_interval gt 0 and checked mod progress_interval eq 0 then
                    print "progress", checked, "nonboundary", nonboundary,
                          "linear_roots", linear_roots, "contact_hits", contact_hits,
                          "exact36", exact_hits36, "exact66", exact_hits66;
                end if;
                if v eq 0 or U^2 - 4*v^2 eq 0 then
                    continue;
                end if;
                nonboundary +:= 1;

                F1_0, F2_0, F3_0, q0, h30, f0, h60 :=
                    CubicContactData(Q!0, b, L, U, v);
                F1_1, F2_1, F3_1, q1, h31, f1, h61 :=
                    CubicContactData(Q!1, b, L, U, v);
                coeff_a := F1_1 - F1_0;
                candidate_as := { Q | };
                if coeff_a eq 0 then
                    if F1_0 ne 0 then
                        continue;
                    end if;
                    special_f1 +:= 1;

                    // Special branch v^3=1: F1 no longer determines a.
                    f3_coeff := F3_1 - F3_0;
                    if f3_coeff ne 0 then
                        Include(~candidate_as, -F3_0/f3_coeff);
                    elif F3_0 eq 0 then
                        F1_2, F2_2, F3_2, q2, h32, f2, h62 :=
                            CubicContactData(Q!2, b, L, U, v);
                        qa2 := (F2_2 - 2*F2_1 + F2_0)/2;
                        qa1 := F2_1 - F2_0 - qa2;
                        qa0 := F2_0;
                        poly := PA![qa0, qa1, qa2];
                        if poly ne 0 then
                            for rt in Roots(poly) do
                                Include(~candidate_as, rt[1]);
                            end for;
                        end if;
                    end if;
                else
                    Include(~candidate_as, -F1_0/coeff_a);
                end if;

                for a in Sort(Setseq(candidate_as)) do
                    F1, F2, F3, q, h3, f, h6 := CubicContactData(a, b, L, U, v);
                    if F1 ne 0 then
                        continue;
                    end if;
                    linear_roots +:= 1;
                    if F2 ne 0 or F3 ne 0 then
                        continue;
                    end if;
                    contact_hits +:= 1;
                    if contact_hits le max_contact_print then
                        print "CONTACT", "a", a, "b", b, "L", L, "U", U, "v", v;
                        print " f", f;
                        print " q", q;
                        print " h3modq", h3 mod q;
                    end if;
                    ok, invs, ordD, ordE, f0 := ExactCheck(a, b, L, U, v);
                    if not ok then
                        continue;
                    end if;
                    is36 := ordD eq 6 and ordE eq 3 and Has36(invs);
                    is66 := is36 and Has66(invs);
                    if is36 then
                        exact_hits36 +:= 1;
                    end if;
                    if is66 then
                        exact_hits66 +:= 1;
                    end if;
                    Append(~hits, <a,b,L,U,v,invs,ordD,ordE,f0>);
                    print "HIT", "a", a, "b", b, "L", L, "U", U, "v", v,
                          "invs", invs, "ordD", ordD, "ordE", ordE;
                    print " f", f0;
                    if #hits ge max_hits then
                        break b;
                    end if;
                end for;
            end for;
        end for;
    end for;
end for;

print "Done";
print "checked", checked;
print "nonboundary", nonboundary;
print "linear_roots", linear_roots;
print "special_f1", special_f1;
print "contact_hits", contact_hits;
print "exact_hits36", exact_hits36;
print "exact_hits66", exact_hits66;
print "hits", #hits;
for H in hits do
    print "H", H;
end for;

quit;
