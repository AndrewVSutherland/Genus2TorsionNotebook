//////////////////////////////////////////////////////////////////////
//  Component analysis for the [3,6] cover on the contact-6 M(6) chart.
//
//  We start from the cubic-contact equations
//
//      F3 = 0, F2 = 0, F1 = 0
//
//  in (a,b,L,U,v).  The equation F1 is linear in a with coefficient
//
//      8*L^2*(v^3 - 1).
//
//  Thus the generic branch v^3 != 1 is analyzed by eliminating a
//  explicitly.  The rational special branch is v=1 and is handled
//  separately by direct finite-field counts.
//
//  Typical runs:
//      magma -b mode:=summary code/contact6_m36_component_analysis.m
//      magma -b mode:=finite prime_bound:=31 code/contact6_m36_component_analysis.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned mode then
    mode := "summary";
end if;
if not assigned output_file then
    output_file := "data/contact6_m36_component_analysis.txt";
end if;
if not assigned prime_bound then
    prime_bound := 31;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;
if not assigned decompose_prime then
    decompose_prime := 5;
elif Type(decompose_prime) eq MonStgElt then
    decompose_prime := StringToInteger(decompose_prime);
end if;
if not assigned do_primary then
    do_primary := true;
elif Type(do_primary) eq MonStgElt then
    do_primary := do_primary in {"true", "True", "1", "yes"};
end if;

Q := Rationals();
Z := Integers();

function PrimitivePolynomial(f)
    if f eq 0 then
        return f;
    end if;
    coeffs := Coefficients(f);
    den := LCM([Denominator(c) : c in coeffs]);
    g := Parent(f)!(den*f);
    nums := [Z!c : c in Coefficients(g)];
    cont := GCD([Abs(n) : n in nums | n ne 0]);
    if cont gt 1 then
        g := Parent(f)!(g/cont);
    end if;
    return g;
end function;

function RepeatedBoundaryRemoval(f, boundary)
    core := f;
    removed := Parent(f)!1;
    while core ne 0 do
        g := GCD(core, boundary);
        if TotalDegree(g) eq 0 then
            break;
        end if;
        removed *:= g;
        core := ExactQuotient(core, g);
    end while;
    return core, removed;
end function;

function FactorSummary(f)
    if f eq 0 then
        return "0";
    end if;
    fac := Factorization(f);
    parts := [];
    for fe in fac do
        Append(~parts, Sprintf("(deg %o, terms %o)^%o",
                               TotalDegree(fe[1]), #Terms(fe[1]), fe[2]));
    end for;
    return Sprint(parts);
end function;

procedure GenericBranchSummary(file)
    S<b,L,U,v> := PolynomialRing(Q, 4);
    K := FieldOfFractions(S);
    PA<a> := PolynomialRing(K);

    c1 := 2*a + 6;
    c2 := a^2 + 2*b - 15;
    c3 := 2*a*b + 22;
    c4 := 2*a + b^2 - 15;
    c5 := 2*b + 6;

    B := c5*L^2 + 3*U;
    Delta := 4*c4*L^2 + 12*(U^2 + v^2) - B^2;
    F3 := B*Delta + 16*v^3 - 8*c3*L^2 - 8*U^3 - 48*U*v^2;
    F2 := Delta^2 + 64*B*v^3 - 64*c2*L^2
          - 192*(U^2*v^2 + v^4);
    F1 := Delta*v^3 - 4*c1*L^2 - 12*U*v^4;

    D := Coefficient(F1, 1);
    N := -Coefficient(F1, 0);
    assert D eq 8*L^2*(v^3 - 1);

    function SubstituteANumerator(P)
        d := Degree(P);
        value := K!0;
        for i in [0..d] do
            value +:= Coefficient(P, i)*N^i*D^(d-i);
        end for;
        return PrimitivePolynomial(S!Numerator(value));
    end function;

    G3 := SubstituteANumerator(F3);
    G2 := SubstituteANumerator(F2);
    common := GCD(G2, G3);
    G2red := ExactQuotient(G2, common);
    G3red := ExactQuotient(G3, common);

    h1_num := PrimitivePolynomial(S!(N + (b+2)*D));  // numerator of a+b+2
    boundary := L*v*(v^3 - 1)*(U^2 - 4*v^2)*(b+3)*h1_num;
    G2core, G2bd := RepeatedBoundaryRemoval(G2red, boundary);
    G3core, G3bd := RepeatedBoundaryRemoval(G3red, boundary);

    out := Open(file, "w");
    fprintf out, "Generic branch v^3 != 1\n";
    fprintf out, "Variables: b,L,U,v\n";
    fprintf out, "F1 coefficient in a D = %o\n", D;
    fprintf out, "a numerator N = %o\n", PrimitivePolynomial(S!N);
    fprintf out, "a = N/D\n";
    fprintf out, "h6(1) numerator after substitution = %o\n\n", h1_num;

    fprintf out, "G2 = numerator(D^2 * F2(a=N/D))\n";
    fprintf out, "  total_degree %o, degrees %o, terms %o\n",
            TotalDegree(G2), [Degree(G2,i) : i in [1..4]], #Terms(G2);
    fprintf out, "  factor_summary %o\n", FactorSummary(G2);
    fprintf out, "G3 = numerator(D * F3(a=N/D))\n";
    fprintf out, "  total_degree %o, degrees %o, terms %o\n",
            TotalDegree(G3), [Degree(G3,i) : i in [1..4]], #Terms(G3);
    fprintf out, "  factor_summary %o\n\n", FactorSummary(G3);

    fprintf out, "Common gcd(G2,G3)\n";
    fprintf out, "  %o\n", common;
    fprintf out, "  total_degree %o, factor_summary %o\n\n",
            TotalDegree(common), FactorSummary(common);

    fprintf out, "After removing gcd and repeated nonboundary factors\n";
    fprintf out, "Boundary product = L*v*(v^3-1)*(U^2-4*v^2)*(b+3)*h1_num\n";
    fprintf out, "G2 boundary removed factor = %o\n", G2bd;
    fprintf out, "G2 core total_degree %o, degrees %o, terms %o\n",
            TotalDegree(G2core), [Degree(G2core,i) : i in [1..4]], #Terms(G2core);
    fprintf out, "G2 core factor_summary %o\n", FactorSummary(G2core);
    fprintf out, "G2 core = %o\n\n", G2core;
    fprintf out, "G3 boundary removed factor = %o\n", G3bd;
    fprintf out, "G3 core total_degree %o, degrees %o, terms %o\n",
            TotalDegree(G3core), [Degree(G3core,i) : i in [1..4]], #Terms(G3core);
    fprintf out, "G3 core factor_summary %o\n", FactorSummary(G3core);
    fprintf out, "G3 core = %o\n", G3core;
    delete out;

    print "Generic branch v^3 != 1";
    print "D =", D;
    print "a numerator degree", TotalDegree(PrimitivePolynomial(S!N)),
          "terms", #Terms(PrimitivePolynomial(S!N));
    print "h6(1) numerator degree", TotalDegree(h1_num), "terms", #Terms(h1_num);
    print "G2 total_degree", TotalDegree(G2), "terms", #Terms(G2),
          "factor_summary", FactorSummary(G2);
    print "G3 total_degree", TotalDegree(G3), "terms", #Terms(G3),
          "factor_summary", FactorSummary(G3);
    print "common total_degree", TotalDegree(common),
          "factor_summary", FactorSummary(common);
    print "G2core total_degree", TotalDegree(G2core), "terms", #Terms(G2core),
          "factor_summary", FactorSummary(G2core);
    print "G3core total_degree", TotalDegree(G3core), "terms", #Terms(G3core),
          "factor_summary", FactorSummary(G3core);
    print "wrote", file;
end procedure;

function ContactDataFinite(F, a, b, L, U, v)
    PF<x> := PolynomialRing(F);
    h6 := 1 + a*x + b*x^2 + x^3;
    f := h6^2 - (x-1)^6;
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
    return F1, F2, F3, f, h6, q;
end function;

function IsGoodFinite(F, a, b, L, U, v)
    F1,F2,F3,f,h6,q := ContactDataFinite(F,a,b,L,U,v);
    if F1 ne 0 or F2 ne 0 or F3 ne 0 then
        return false;
    end if;
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false;
    end if;
    if Evaluate(h6, F!1) eq 0 then
        return false;
    end if;
    if Discriminant(q) eq 0 or Degree(GCD(q, f)) gt 0 then
        return false;
    end if;
    return true;
end function;


procedure SpecialBranchSummary(file)
    S<b,L,U> := PolynomialRing(Q, 3);
    K := FieldOfFractions(S);
    PA<a> := PolynomialRing(K);
    v := K!1;

    c1 := 2*a + 6;
    c2 := a^2 + 2*b - 15;
    c3 := 2*a*b + 22;
    c4 := 2*a + b^2 - 15;
    c5 := 2*b + 6;

    B := c5*L^2 + 3*U;
    Delta := 4*c4*L^2 + 12*(U^2 + v^2) - B^2;
    F3 := B*Delta + 16*v^3 - 8*c3*L^2 - 8*U^3 - 48*U*v^2;
    F2 := Delta^2 + 64*B*v^3 - 64*c2*L^2
          - 192*(U^2*v^2 + v^4);
    F1 := Delta*v^3 - 4*c1*L^2 - 12*U*v^4;

    assert Degree(F1) eq 0;
    H := PrimitivePolynomial(S!Coefficient(F1, 0));
    D3 := Coefficient(F3, 1);
    N3 := -Coefficient(F3, 0);

    function SubstituteANumerator(P)
        d := Degree(P);
        value := K!0;
        for i in [0..d] do
            value +:= Coefficient(P, i)*N3^i*D3^(d-i);
        end for;
        return PrimitivePolynomial(S!Numerator(value));
    end function;

    J := SubstituteANumerator(F2);
    h1_num := PrimitivePolynomial(S!(N3 + (b+2)*D3));
    boundary := L*(U^2-4)*(b+3)*h1_num*PrimitivePolynomial(S!D3);
    Hcore, Hbd := RepeatedBoundaryRemoval(H, boundary);
    Jcore, Jbd := RepeatedBoundaryRemoval(J, boundary);

    out := Open(file, "a");
    fprintf out, "\n\nSpecial rational branch v=1\n";
    fprintf out, "Variables: b,L,U\n";
    fprintf out, "F1(v=1) = H = %o\n", H;
    fprintf out, "H factor_summary %o\n", FactorSummary(H);
    fprintf out, "F3 coefficient in a D3 = %o\n", PrimitivePolynomial(S!D3);
    fprintf out, "F3 a numerator N3 = %o\n", PrimitivePolynomial(S!N3);
    fprintf out, "h6(1) numerator after F3 substitution = %o\n", h1_num;
    fprintf out, "J = numerator(D3^2 * F2(a=N3/D3, v=1))\n";
    fprintf out, "J total_degree %o, degrees %o, terms %o\n",
            TotalDegree(J), [Degree(J,i) : i in [1..3]], #Terms(J);
    fprintf out, "J factor_summary %o\n", FactorSummary(J);
    fprintf out, "H boundary removed factor = %o\n", Hbd;
    fprintf out, "H core = %o\n", Hcore;
    fprintf out, "H core factor_summary %o\n", FactorSummary(Hcore);
    fprintf out, "J boundary removed factor = %o\n", Jbd;
    fprintf out, "J core total_degree %o, degrees %o, terms %o\n",
            TotalDegree(Jcore), [Degree(Jcore,i) : i in [1..3]], #Terms(Jcore);
    fprintf out, "J core factor_summary %o\n", FactorSummary(Jcore);
    fprintf out, "J core = %o\n", Jcore;

    fprintf out, "\nSpecial branch components by H factor\n";
    for fe in Factorization(Hcore) do
        hfac := fe[1];
        gj := GCD(hfac, Jcore);
        fprintf out, "H factor degree %o terms %o exponent %o; gcd_with_Jcore_degree %o\n",
                TotalDegree(hfac), #Terms(hfac), fe[2], TotalDegree(gj);
        fprintf out, "  Hfac = %o\n", hfac;
    end for;
    delete out;

    print "Special branch v=1";
    print "H total_degree", TotalDegree(H), "terms", #Terms(H),
          "factor_summary", FactorSummary(H);
    print "D3 total_degree", TotalDegree(PrimitivePolynomial(S!D3)),
          "factor_summary", FactorSummary(PrimitivePolynomial(S!D3));
    print "J total_degree", TotalDegree(J), "terms", #Terms(J),
          "factor_summary", FactorSummary(J);
    print "Hcore total_degree", TotalDegree(Hcore), "terms", #Terms(Hcore),
          "factor_summary", FactorSummary(Hcore);
    print "Jcore total_degree", TotalDegree(Jcore), "terms", #Terms(Jcore),
          "factor_summary", FactorSummary(Jcore);
end procedure;

procedure FiniteCounts()
    print "Finite-field nonboundary cubic-contact counts";
    print "prime_bound", prime_bound;
    for p in [p : p in PrimesUpTo(prime_bound) | p notin {2,3}] do
        F := GF(p);
        generic_contact := 0;
        generic_good := 0;
        special_contact := 0;
        special_good := 0;
        boundary_bminus3 := 0;
        sample_generic := [];
        sample_special := [];

        // Generic branch: v^3 != 1, so F1 determines a.
        for b in F do
            for L in F do
                if L eq 0 then
                    continue;
                end if;
                for U in F do
                    for v in F do
                        if v eq 0 or v^3 eq 1 or U^2 - 4*v^2 eq 0 then
                            continue;
                        end if;
                        F1_0,F2_0,F3_0,f0,h0,q0 := ContactDataFinite(F,F!0,b,L,U,v);
                        F1_1,F2_1,F3_1,f1,h1,q1 := ContactDataFinite(F,F!1,b,L,U,v);
                        coeff := F1_1 - F1_0;
                        if coeff eq 0 then
                            continue;
                        end if;
                        a := -F1_0/coeff;
                        F1,F2,F3,f,h6,q := ContactDataFinite(F,a,b,L,U,v);
                        if F1 eq 0 and F2 eq 0 and F3 eq 0 then
                            generic_contact +:= 1;
                            if b eq -F!3 then
                                boundary_bminus3 +:= 1;
                            end if;
                            good := IsGoodFinite(F,a,b,L,U,v);
                            if good then
                                generic_good +:= 1;
                                if #sample_generic lt 5 then
                                    Append(~sample_generic,
                                           <Z!a,Z!b,Z!L,Z!U,Z!v>);
                                end if;
                            end if;
                        end if;
                    end for;
                end for;
            end for;
        end for;

        // Rational special branch: v=1.
        v := F!1;
        for a in F do
            for b in F do
                for L in F do
                    if L eq 0 then
                        continue;
                    end if;
                    for U in F do
                        if U^2 - 4 eq 0 then
                            continue;
                        end if;
                        F1,F2,F3,f,h6,q := ContactDataFinite(F,a,b,L,U,v);
                        if F1 eq 0 and F2 eq 0 and F3 eq 0 then
                            special_contact +:= 1;
                            good := IsGoodFinite(F,a,b,L,U,v);
                            if good then
                                special_good +:= 1;
                                if #sample_special lt 5 then
                                    Append(~sample_special,
                                           <Z!a,Z!b,Z!L,Z!U,Z!v>);
                                end if;
                            end if;
                        end if;
                    end for;
                end for;
            end for;
        end for;

        print "p", p,
              "generic_contact", generic_contact,
              "generic_good", generic_good,
              "generic_b=-3", boundary_bminus3,
              "special_v=1_contact", special_contact,
              "special_v=1_good", special_good;
        print " sample_generic", sample_generic;
        print " sample_special", sample_special;
    end for;
end procedure;


procedure PrintIdealComponents(label, I, primary)
    print label;
    Isat := I;
    dim, component_degrees := Dimension(Isat);
    print " dimension", dim, "component_degrees", component_degrees;
    if not primary then
        return;
    end if;
    try
        print " degree", Degree(Isat);
    catch e
        print " degree unavailable";
    end try;
    comps := PrimaryDecomposition(Isat);
    print " components", #comps;
    for i in [1..#comps] do
        C := comps[i];
        cdim, cdegseq := Dimension(C);
        print "  component", i,
              "dimension", cdim,
              "component_degrees", cdegseq,
              "degree", Degree(C),
              "prime", IsPrime(C),
              "basis_degrees", [TotalDegree(g) : g in Basis(C)];
    end for;
end procedure;

procedure DecomposeFinite(p, primary)
    F := GF(p);
    print "Finite saturated component decomposition";
    print "p", p;

    // Generic branch v^3 != 1.
    S<b,L,U,v> := PolynomialRing(F, 4);
    K := FieldOfFractions(S);
    PA<a> := PolynomialRing(K);
    c1 := 2*a + 6;
    c2 := a^2 + 2*b - 15;
    c3 := 2*a*b + 22;
    c4 := 2*a + b^2 - 15;
    c5 := 2*b + 6;
    B := c5*L^2 + 3*U;
    Delta := 4*c4*L^2 + 12*(U^2 + v^2) - B^2;
    F3 := B*Delta + 16*v^3 - 8*c3*L^2 - 8*U^3 - 48*U*v^2;
    F2 := Delta^2 + 64*B*v^3 - 64*c2*L^2
          - 192*(U^2*v^2 + v^4);
    F1 := Delta*v^3 - 4*c1*L^2 - 12*U*v^4;
    D := Coefficient(F1, 1);
    N := -Coefficient(F1, 0);

    function SubNum(P)
        d := Degree(P);
        value := K!0;
        for i in [0..d] do
            value +:= Coefficient(P, i)*N^i*D^(d-i);
        end for;
        return S!Numerator(value);
    end function;

    G3 := SubNum(F3);
    G2 := SubNum(F2);
    common := GCD(G2, G3);
    G2red := ExactQuotient(G2, common);
    G3red := ExactQuotient(G3, common);
    h1_num := S!(N + (b+2)*D);
    boundary := L*v*(v^3 - 1)*(U^2 - 4*v^2)*(b+3)*h1_num;
    Igeneric := ideal<S | G2red, G3red>;
    Igeneric_sat := Saturation(Igeneric, ideal<S | boundary>);
    PrintIdealComponents("generic_branch", Igeneric_sat, primary);

    // Special branch v=1.
    S2<bb,LL,UU> := PolynomialRing(F, 3);
    K2 := FieldOfFractions(S2);
    PA2<aa> := PolynomialRing(K2);
    vv := K2!1;
    c1s := 2*aa + 6;
    c2s := aa^2 + 2*bb - 15;
    c3s := 2*aa*bb + 22;
    c4s := 2*aa + bb^2 - 15;
    c5s := 2*bb + 6;
    Bs := c5s*LL^2 + 3*UU;
    Deltas := 4*c4s*LL^2 + 12*(UU^2 + vv^2) - Bs^2;
    F3s := Bs*Deltas + 16*vv^3 - 8*c3s*LL^2 - 8*UU^3 - 48*UU*vv^2;
    F2s := Deltas^2 + 64*Bs*vv^3 - 64*c2s*LL^2
           - 192*(UU^2*vv^2 + vv^4);
    F1s := Deltas*vv^3 - 4*c1s*LL^2 - 12*UU*vv^4;
    H := S2!Coefficient(F1s, 0);
    D3 := Coefficient(F3s, 1);
    N3 := -Coefficient(F3s, 0);

    d := Degree(F2s);
    value := K2!0;
    for i in [0..d] do
        value +:= Coefficient(F2s, i)*N3^i*D3^(d-i);
    end for;
    J := S2!Numerator(value);
    h1s := S2!(N3 + (bb+2)*D3);
    boundary_s := LL*(UU^2 - 4)*(bb+3)*h1s*(S2!D3);
    Ispecial := ideal<S2 | H, J>;
    Ispecial_sat := Saturation(Ispecial, ideal<S2 | boundary_s>);
    PrintIdealComponents("special_v1_branch", Ispecial_sat, primary);
end procedure;

if mode eq "summary" then
    GenericBranchSummary(output_file);
    SpecialBranchSummary(output_file);
elif mode eq "finite" then
    FiniteCounts();
elif mode eq "decompose" then
    DecomposeFinite(decompose_prime, do_primary);
elif mode eq "dims" then
    DecomposeFinite(decompose_prime, false);
else
    error "unknown mode";
end if;

quit;
