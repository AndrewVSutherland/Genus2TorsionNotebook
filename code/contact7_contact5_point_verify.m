//////////////////////////////////////////////////////////////////////
//  Verify candidates from contact7_contact5_point_enum.py.
//
//  Input rows:
//      a b r c0 c1 c2 d e
//
//  The script checks:
//      contact-7 marked divisor has order 7,
//      contact-5 marked divisor at x=r has order 5,
//      exact torsion contains an invariant divisible by 35,
//      optional Q-simplicity certificate.
//////////////////////////////////////////////////////////////////////

if not assigned input_file then
    input_file := "data/contact7_contact5_point_candidates_h8.txt";
end if;
if not assigned max_exact then
    max_exact := 200;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;

Q := Rationals();
P<x> := PolynomialRing(Q);

function ParseRational(s)
    if "/" in s then
        parts := Split(s, "/");
        return Q!StringToInteger(parts[1]) / Q!StringToInteger(parts[2]);
    end if;
    return Q!StringToInteger(s);
end function;

function IntegralModel(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function TorsionOrder(invs)
    if #invs eq 0 then
        return 1;
    end if;
    return &*invs;
end function;

function HasInvariantDivisibleBy(invs, n)
    return &or [ m mod n eq 0 : m in invs ];
end function;

function SimpleCertificate(f)
    C := HyperellipticCurve(f);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149] do
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

lines := Split(Read(input_file), "\n");
rows := [];
for line in lines do
    s := StripWhiteSpace(line);
    if #s eq 0 or s[1] eq "#" then
        continue;
    end if;
    Append(~rows, Split(s, " "));
end for;

print "VERIFY simultaneous contact7/contact5 point candidates";
print "input_file", input_file, "rows", #rows, "max_exact", max_exact;

checked := 0;
smooth := 0;
verified_contacts := 0;
exact_tests := 0;
hits := 0;
torsion_counts := AssociativeArray();

for row in rows do
    checked +:= 1;
    if exact_tests ge max_exact then
        continue;
    end if;

    a := ParseRational(row[1]);
    b := ParseRational(row[2]);
    r := ParseRational(row[3]);
    c0 := ParseRational(row[4]);
    c1 := ParseRational(row[5]);
    c2 := ParseRational(row[6]);

    h7 := 1 - (Q!7/2)*x + a*x^2 + b*x^3;
    f := ExactQuotient(h7^2 + (x - 1)^7, x^2);
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        continue;
    end if;
    smooth +:= 1;

    q := c0 + c1*x + c2*x^2;
    contact5 := q^2 - f;
    if contact5 ne LeadingCoefficient(contact5)*(x-r)^5 then
        print "BAD_CONTACT_POLY", row;
        continue;
    end if;

    C := HyperellipticCurve(f);
    J := Jacobian(C);
    D7 := J![x - 1, Evaluate(h7, Q!1)];
    D5 := J![x - r, Evaluate(q, r)];
    if Order(D7) ne 7 or Order(D5) ne 5 then
        print "BAD_ORDER", "a", a, "b", b, "r", r,
              "ord7", Order(D7), "ord5", Order(D5);
        continue;
    end if;
    verified_contacts +:= 1;

    fI, L := IntegralModel(f);
    CI := HyperellipticCurve(fI);
    JI := Jacobian(CI);
    G, phi := TorsionSubgroup(JI);
    invs := Invariants(G);
    exact_tests +:= 1;
    key := Sprint(invs);
    if IsDefined(torsion_counts, key) then
        torsion_counts[key] +:= 1;
    else
        torsion_counts[key] := 1;
    end if;

    if HasInvariantDivisibleBy(invs, 35) then
        hits +:= 1;
        simple, pcert, Lp := SimpleCertificate(fI);
        print "HIT35", "a", a, "b", b, "r", r,
              "torsion", invs, "order", TorsionOrder(invs),
              "simple", simple, "pcert", pcert;
        print "  f =", fI;
        print "  q =", q;
    else
        print "CONTACTS_NOT_35", "a", a, "b", b, "r", r, "torsion", invs;
    end if;
end for;

print "DONE verify simultaneous contacts";
print "checked", checked, "smooth", smooth,
      "verified_contacts", verified_contacts,
      "exact_tests", exact_tests, "hits", hits;
print "Torsion counts";
for key in Sort([ k : k in Keys(torsion_counts) ]) do
    print " ", key, torsion_counts[key];
end for;

quit;
