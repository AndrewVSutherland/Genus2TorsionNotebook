//////////////////////////////////////////////////////////////////////
//  Verify rational fibers of the simultaneous contact-5/contact-6 cover.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned input_file then
    input_file := "data/contact5_contact6_order30_fibers_h500.txt";
end if;
if not assigned output_file then
    output_file := "data/contact5_contact6_order30_verify_h500.txt";
end if;
if not assigned max_exact then
    max_exact := 20;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

function ParseRational(s)
    if "/" in s then
        parts := Split(s, "/");
        return Q!StringToInteger(parts[1]) / Q!StringToInteger(parts[2]);
    end if;
    return Q!StringToInteger(s);
end function;

function IntegralModel(f)
    L := Z!1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function SimpleCertificate(f)
    C := HyperellipticCurve(f);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199] do
        try
            fp := ChangeRing(f, GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
                continue;
            end if;
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

function TorsionOrder(invs)
    if #invs eq 0 then
        return 1;
    end if;
    return &*invs;
end function;

function Reconstruct(u, s)
    Cc := (u^2 + 1)/(2*u);
    c := (u^2 - 1)/(2*u);
    den := u^6 + 6*u^4*s - 2*u^4 + 15*u^3*s - u*s^3 + u^2;
    num := 15*u^5 + 90*u^4 + 20*u^3*s - 6*u^2*s^2 + 231*u^3
        + 2*u^2*s - 15*u*s^2 + 90*u^2 - 20*u*s + 15*u - 2*s;
    if den eq 0 or c eq 0 then
        return false, 0, 0, 0, 0, 0, 0, P!0, P!0, P!0;
    end if;
    q := num/den;
    A := (s + q)/2;
    e := (s - q)/2;
    B := (15 - s*q)/2;
    d := (B*Cc + 3)/c;
    h6 := x^3 + A*x^2 + B*x + Cc;
    h5 := e*x^2 + d*x + c;
    f := h6^2 - (x-1)^6;
    return true, A, B, Cc, e, d, c, h6, h5, f;
end function;

rows := [];
for raw_line in Split(Read(input_file), "\n") do
    line := StripWhiteSpace(raw_line);
    if #line eq 0 or line[1] eq "#" then
        continue;
    end if;
    parts := Split(line, " ");
    if #parts lt 2 then
        continue;
    end if;
    Append(~rows, <ParseRational(parts[1]), ParseRational(parts[2])>);
end for;

if #rows eq 0 then
    rows := [
        <Q!125, Q!5415>,
        <Q!125, Q!2715>,
        <Q!1/8, Q!21/32>,
        <Q!1/125, Q!831/3125>,
        <Q!1/125, -Q!69/3125>
    ];
end if;

out := Open(output_file, "w");
fprintf out, "# verify simultaneous contact5/contact6 rational fibers\n";
fprintf out, "# input_file %o rows %o max_exact %o\n", input_file, #rows, max_exact;

exact_tests := 0;
for row in rows do
    u := row[1];
    s := row[2];
    ok, A, B, Cc, e, d, c, h6, h5, f := Reconstruct(u, s);
    fprintf out, "============================================================\n";
    fprintf out, "fiber u=%o s=%o\n", u, s;
    if not ok then
        fprintf out, "SKIP reconstruction denominator zero or c=0\n";
        continue;
    end if;
    fprintf out, "A=%o B=%o C=%o e=%o d=%o c=%o\n", A, B, Cc, e, d, c;
    fprintf out, "h6=%o\n", h6;
    fprintf out, "h5=%o\n", h5;
    fprintf out, "contact5 check h5^2-f over x^5 = %o\n", ExactQuotient(h5^2 - f, x^5);
    fprintf out, "contact6 h6(1)=%o h5(0)=%o\n", Evaluate(h6,1), Evaluate(h5,0);
    fprintf out, "f=%o\n", f;

    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        fprintf out, "SINGULAR_OR_BAD_DEGREE degree=%o disc=%o\n", Degree(f), Discriminant(f);
        continue;
    end if;

    Ccurve := HyperellipticCurve(f);
    J := Jacobian(Ccurve);
    D5 := J![x, c];
    D6 := J![x-1, Evaluate(h6,1)];
    ord5 := Order(D5);
    ord6 := Order(D6);
    ordsum := Order(D5 + D6);
    fprintf out, "orders D5=%o D6=%o D5+D6=%o\n", ord5, ord6, ordsum;

    if exact_tests lt max_exact then
        fI, scale := IntegralModel(f);
        CI := HyperellipticCurve(fI);
        JI := Jacobian(CI);
        G, phi := TorsionSubgroup(JI);
        invs := Invariants(G);
        simple, pcert, Lp := SimpleCertificate(fI);
        exact_tests +:= 1;
        fprintf out, "integral_scale=%o\n", scale;
        fprintf out, "integral_f=%o\n", fI;
        fprintf out, "torsion_invariants=%o torsion_order=%o\n", invs, TorsionOrder(invs);
        fprintf out, "simple_certificate=%o p=%o Lp=%o\n", simple, pcert, Lp;
    end if;
end for;

fprintf out, "DONE exact_tests %o\n", exact_tests;
delete out;

print "Wrote", output_file;
quit;
