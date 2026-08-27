//////////////////////////////////////////////////////////////////////
//  Try to extract a rational one-parameter [6,6] family from the
//  contact-6 [1,2,2] core cover.
//
//  We use the generic branch v^3 != 1, where F1 is linear in a.  After
//  substituting a=N/D, the remaining equations are two polynomials
//  G2(b,L,U,v)=G3(b,L,U,v)=0.  The known simple exact [6,6] point is
//
//      b=-7/13, L=29/16, U=-9/4, v=5/2.
//
//  A coordinate hyperplane through this point gives a curve.  This script
//  probes those curves by resultants and, when feasible, by curve invariants.
//
//  Typical run:
//      magma code/contact6_m36_core_family_slices.m \
//          > data/contact6_m36_core_family_slices.txt
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
Z := Integers();

if not assigned genus_bound then
    genus_bound := 6;
elif Type(genus_bound) eq MonStgElt then
    genus_bound := StringToInteger(genus_bound);
end if;

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

function BuildGenericCore()
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

    G2 := SubstituteANumerator(F2);
    G3 := SubstituteANumerator(F3);
    common := GCD(G2, G3);
    G2red := ExactQuotient(G2, common);
    G3red := ExactQuotient(G3, common);

    h1_num := PrimitivePolynomial(S!(N + (b+2)*D));
    boundary := L*v*(v^3 - 1)*(U^2 - 4*v^2)*(b+3)*h1_num;
    G2core, G2bd := RepeatedBoundaryRemoval(G2red, boundary);
    G3core, G3bd := RepeatedBoundaryRemoval(G3red, boundary);

    return S, G2core, G3core, N, D, h1_num;
end function;

function EvalRat(r, pt)
    num := Numerator(r);
    den := Denominator(r);
    return Evaluate(num, pt)/Evaluate(den, pt);
end function;

procedure PrintEval(label, f, pt)
    print label, Evaluate(f, pt);
end procedure;

procedure AnalyzePlaneCurve(label, R, F, G, known2)
    print "";
    print "Plane projection", label;
    print "  R rank", Rank(R);
    print "  F degree", TotalDegree(F), "terms", #Terms(F);
    print "  G degree", TotalDegree(G), "terms", #Terms(G);
    print "  gcd degree", TotalDegree(GCD(F,G)),
          "summary", FactorSummary(GCD(F,G));
    print "  F(known)", Evaluate(F, known2), "G(known)", Evaluate(G, known2);
end procedure;

function HomogenizePlanePolynomial(f)
    d := TotalDegree(f);
    P2<X,Y,Z> := ProjectiveSpace(Q, 2);
    R3 := CoordinateRing(P2);
    Fh := R3!0;
    mons := Monomials(f);
    coeffs := Coefficients(f);
    for i in [1..#mons] do
        exps := Exponents(mons[i]);
        Fh +:= coeffs[i]*R3.1^exps[1]*R3.2^exps[2]*R3.3^(d-exps[1]-exps[2]);
    end for;
    return P2, Fh;
end function;

procedure AnalyzePlaneFactor(label, f, known2)
    print "  factor polynomial", f;
    P2, Fh := HomogenizePlanePolynomial(f);
    print "  projective known value", Evaluate(Fh, [known2[1], known2[2], Q!1]);
    try
        C := Curve(P2, Fh);
        print "  plane factor genus", Genus(C);
        try
            sing := SingularPoints(C);
            print "  singular_points", #sing;
            if #sing le 8 then
                print "  singular_points_list", sing;
            end if;
        catch e
            print "  singular_points unavailable", e`Object;
        end try;
        if Genus(C) eq 0 then
            try
                par := Parametrization(C);
                print "  parametrization", par;
            catch e
                print "  parametrization unavailable", e`Object;
            end try;
        end if;
    catch e
        print "  plane factor curve unavailable", e`Object;
    end try;
end procedure;

procedure AnalyzeSlice(label, fixed_index, fixed_value, G2, G3)
    names := ["b", "L", "U", "v"];
    print "";
    print "============================================================";
    print "Slice", label, names[fixed_index], "=", fixed_value;

    // Make a 3-variable coordinate ring for the slice.
    R3<X1,X2,X3> := PolynomialRing(Q, 3);
    vals := [* R3!0, R3!0, R3!0, R3!0 *];
    j := 1;
    for i in [1..4] do
        if i eq fixed_index then
            vals[i] := Q!fixed_value;
        else
            vals[i] := R3.j;
            j +:= 1;
        end if;
    end for;

    F := PrimitivePolynomial(Evaluate(G2, <vals[1], vals[2], vals[3], vals[4]>));
    G := PrimitivePolynomial(Evaluate(G3, <vals[1], vals[2], vals[3], vals[4]>));
    H := GCD(F,G);
    if TotalDegree(H) gt 0 then
        F := ExactQuotient(F,H);
        G := ExactQuotient(G,H);
    end if;

    print "slice variables", [names[i] : i in [1..4] | i ne fixed_index];
    print "F degree", TotalDegree(F), "degrees", [Degree(F,i) : i in [1..3]],
          "terms", #Terms(F), "factor", FactorSummary(F);
    print "G degree", TotalDegree(G), "degrees", [Degree(G,i) : i in [1..3]],
          "terms", #Terms(G), "factor", FactorSummary(G);
    print "common removed degree", TotalDegree(H), "summary", FactorSummary(H);

    b0 := Q!(-7/13); L0 := Q!(29/16); U0 := Q!(-9/4); v0 := Q!(5/2);
    known_all := [b0, L0, U0, v0];
    known3 := <known_all[i] : i in [1..4] | i ne fixed_index>;
    print "known3", known3;
    print "F(known)", Evaluate(F, known3), "G(known)", Evaluate(G, known3);

    I := ideal<R3 | F, G>;
    try
        dim, degs := Dimension(I);
        print "ideal dimension", dim, "component_degrees", degs;
        print "ideal degree", Degree(I);
    catch e
        print "ideal dimension/degree unavailable", e`Object;
    end try;

    // Projection resultants.  For a low-genus rational family, a low-degree
    // factor through the known projected point is the main thing to look for.
    for elim in [1..3] do
        keep := [i : i in [1..3] | i ne elim];
        R2<Y1,Y2> := PolynomialRing(Q, 2);
        R1<T> := PolynomialRing(R2);

        imgs := [];
        for ii in [1..3] do
            if ii eq elim then
                Append(~imgs, T);
            elif ii eq keep[1] then
                Append(~imgs, R1!(R2.1));
            else
                Append(~imgs, R1!(R2.2));
            end if;
        end for;
        embF := hom<R3 -> R1 | imgs[1], imgs[2], imgs[3]>;
        RF := embF(F);
        RG := embF(G);
        res := Resultant(RF, RG);
        res2 := PrimitivePolynomial(R2!res);
        if res2 ne 0 then
            print "resultant eliminating", elim, "degree", TotalDegree(res2),
                  "terms", #Terms(res2), "factor", FactorSummary(res2);
            known2 := <known3[keep[1]], known3[keep[2]]>;
            fac := Factorization(res2);
            for k in [1..#fac] do
                facval := Evaluate(fac[k][1], known2);
                if facval eq 0 then
                    print "  factor through known projection", k,
                          "degree", TotalDegree(fac[k][1]), "terms", #Terms(fac[k][1]),
                          "exponent", fac[k][2];
                    if TotalDegree(fac[k][1]) le genus_bound then
                        AnalyzePlaneFactor(Sprintf("%o elim %o factor %o", label, elim, k),
                                           fac[k][1], known2);
                    end if;
                end if;
            end for;
        else
            print "resultant eliminating", elim, "is zero";
        end if;
    end for;
end procedure;

S, G2, G3, N, D, h1 := BuildGenericCore();
b0 := Q!(-7/13); L0 := Q!(29/16); U0 := Q!(-9/4); v0 := Q!(5/2);
pt := <b0,L0,U0,v0>;

print "Contact-6 core family-slice analysis";
print "G2 degree", TotalDegree(G2), "degrees", [Degree(G2,i) : i in [1..4]],
      "terms", #Terms(G2), "factor", FactorSummary(G2);
print "G3 degree", TotalDegree(G3), "degrees", [Degree(G3,i) : i in [1..4]],
      "terms", #Terms(G3), "factor", FactorSummary(G3);
PrintEval("G2 known", G2, pt);
PrintEval("G3 known", G3, pt);
print "a recovered", EvalRat(N, pt)/EvalRat(D, pt);
print "h6(1) numerator known", Evaluate(h1, pt);

AnalyzeSlice("fix v at known value", 4, v0, G2, G3);
AnalyzeSlice("fix b at known value", 1, b0, G2, G3);
AnalyzeSlice("fix U at known value", 3, U0, G2, G3);
AnalyzeSlice("fix L at known value", 2, L0, G2, G3);

quit;
