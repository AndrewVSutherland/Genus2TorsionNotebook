//////////////////////////////////////////////////////////////////////
//  Contact-6 [1,2,2] core-cover search for simple [6,6].
//
//  The contact-6 family factors as
//
//      f = h6^2 - (x-1)^6
//        = x*((b+3)*x^2 + (a-3)*x + 2)
//           *(2*x^2 + (b-3)*x + (a+3)).
//
//  We keep the simple-friendly factor type [1,2,2], so the two quadratics
//  are irreducible over Q.  We impose an independent 3-torsion class by
//  cubic contact, but write M=L^2 in the contact equations.  The search
//  scans rational (b,L,U,v), solves the first contact equation for a
//  algebraically, and then checks the remaining two equations.
//
//  Typical runs:
//      magma -b height:=5 code/contact6_m36_core_search.m
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
if not assigned simple_only then
    simple_only := true;
elif Type(simple_only) eq MonStgElt then
    simple_only := simple_only in {"true", "True", "1", "yes"};
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);
PA<z> := PolynomialRing(Q);

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

function FactorDegrees(f)
    return Sort([Degree(fe[1]) : fe in Factorization(f)]);
end function;

function HasFactorType122(f)
    return FactorDegrees(f) eq [1,2,2];
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

function CubicContactDataM(a, b, M, U, v)
    f, h6 := Contact6Polynomial(a, b);
    c1 := Coefficient(f, 1);
    c2 := Coefficient(f, 2);
    c3 := Coefficient(f, 3);
    c4 := Coefficient(f, 4);
    c5 := Coefficient(f, 5);

    B := c5*M + 3*U;
    Delta := 4*c4*M + 12*(U^2 + v^2) - B^2;
    F3 := B*Delta + 16*v^3 - 8*c3*M - 8*U^3 - 48*U*v^2;
    F2 := Delta^2 + 64*B*v^3 - 64*c2*M
          - 192*(U^2*v^2 + v^4);
    F1 := Delta*v^3 - 4*c1*M - 12*U*v^4;
    q := x^2 + U*x + v^2;
    return F1, F2, F3, q, B, Delta, f, h6;
end function;

function CandidateAsFromContact(b, M, U, v)
    F1_0,F2_0,F3_0,q0,B0,D0,f0,h0 := CubicContactDataM(Q!0,b,M,U,v);
    F1_1,F2_1,F3_1,q1,B1,D1,f1,h1 := CubicContactDataM(Q!1,b,M,U,v);
    coeff_a := F1_1 - F1_0;
    candidates := {Q | };

    if coeff_a ne 0 then
        Include(~candidates, -F1_0/coeff_a);
        return candidates, false;
    end if;

    if F1_0 ne 0 then
        return candidates, false;
    end if;

    // Special branch, rationally v=1.  Use F3 if it determines a;
    // otherwise solve the quadratic F2 in a.
    f3_coeff := F3_1 - F3_0;
    if f3_coeff ne 0 then
        Include(~candidates, -F3_0/f3_coeff);
        return candidates, true;
    end if;

    if F3_0 ne 0 then
        return candidates, true;
    end if;

    F1_2,F2_2,F3_2,q2,B2,D2,f2,h2 := CubicContactDataM(Q!2,b,M,U,v);
    c2 := (F2_2 - 2*F2_1 + F2_0)/2;
    c1 := F2_1 - F2_0 - c2;
    c0 := F2_0;
    poly := PA![c0,c1,c2];
    if poly ne 0 then
        for rt in Roots(poly) do
            Include(~candidates, rt[1]);
        end for;
    end if;
    return candidates, true;
end function;

function ExactCheck(a, b, L, U, v)
    M := L^2;
    F1,F2,F3,q,B,Delta,f,h6 := CubicContactDataM(a,b,M,U,v);
    if F1 ne 0 or F2 ne 0 or F3 ne 0 then
        return false, [], 0, 0, false, 0, P!0, f;
    end if;
    if not GoodPolynomial(f) or Evaluate(h6, Q!1) eq 0 then
        return false, [], 0, 0, false, 0, P!0, f;
    end if;
    if not HasFactorType122(f) then
        return false, [], 0, 0, false, 0, P!0, f;
    end if;
    if Discriminant(q) eq 0 or Degree(GCD(q, f)) gt 0 then
        return false, [], 0, 0, false, 0, P!0, f;
    end if;

    simple, pcert, Lp := SimpleCertificate(f);
    if simple_only and not simple then
        return false, [], 0, 0, false, 0, P!0, f;
    end if;

    h3 := (1/L)*x^3 + (B/(2*L))*x^2 + (Delta/(8*L))*x + v^3/L;
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    D := J![x-1, Evaluate(h6, Q!1)];
    E := J![q, h3 mod q];
    ordD := Order(D);
    ordE := Order(E);
    G, phi := TorsionSubgroup(J);
    return true, Invariants(G), ordD, ordE, simple, pcert, Lp, f;
end function;

params := RationalParametersOfHeight(height);
checked := 0;
nonboundary := 0;
candidate_as := 0;
special_branch := 0;
contact_hits := 0;
factor122 := 0;
exact_tests := 0;
hits := [];

print "Contact-6 [1,2,2] core search for simple [6,6]";
print "height", height, "parameter_count", #params,
      "simple_only", simple_only, "max_hits", max_hits;

for b in params do
    for L in params do
        if L eq 0 then
            continue;
        end if;
        M := L^2;
        for U in params do
            for v in params do
                checked +:= 1;
                if progress_interval gt 0 and checked mod progress_interval eq 0 then
                    print "progress", checked,
                          "nonboundary", nonboundary,
                          "candidate_as", candidate_as,
                          "contact_hits", contact_hits,
                          "factor122", factor122,
                          "exact_tests", exact_tests,
                          "hits", #hits;
                end if;
                if v eq 0 or U^2 - 4*v^2 eq 0 then
                    continue;
                end if;
                nonboundary +:= 1;

                as, special := CandidateAsFromContact(b,M,U,v);
                if special then
                    special_branch +:= 1;
                end if;
                for a in Sort(Setseq(as)) do
                    candidate_as +:= 1;
                    F1,F2,F3,q,B,Delta,f,h6 := CubicContactDataM(a,b,M,U,v);
                    if F1 ne 0 or F2 ne 0 or F3 ne 0 then
                        continue;
                    end if;
                    contact_hits +:= 1;
                    if contact_hits le max_contact_print then
                        print "CONTACT", "a", a, "b", b, "L", L, "U", U, "v", v,
                              "factor_degrees", FactorDegrees(f);
                        print " f", f;
                    end if;
                    if not GoodPolynomial(f) or not HasFactorType122(f) then
                        continue;
                    end if;
                    factor122 +:= 1;
                    exact_tests +:= 1;
                    ok, invs, ordD, ordE, simple, pcert, Lp, f0 :=
                        ExactCheck(a,b,L,U,v);
                    if not ok then
                        continue;
                    end if;
                    if ordD eq 6 and ordE eq 3 and Has66(invs) then
                        Append(~hits, <a,b,L,U,v,invs,ordD,ordE,simple,pcert,f0>);
                        print "HIT66", "a", a, "b", b, "L", L, "U", U, "v", v,
                              "invs", invs, "ordD", ordD, "ordE", ordE,
                              "simple", simple, "pcert", pcert;
                        print " f", f0;
                        if #hits ge max_hits then
                            break b;
                        end if;
                    end if;
                end for;
            end for;
        end for;
    end for;
end for;

print "Done";
print "checked", checked;
print "nonboundary", nonboundary;
print "candidate_as", candidate_as;
print "special_branch", special_branch;
print "contact_hits", contact_hits;
print "factor122", factor122;
print "exact_tests", exact_tests;
print "hits", #hits;
for H in hits do
    print "H", H;
end for;

quit;
