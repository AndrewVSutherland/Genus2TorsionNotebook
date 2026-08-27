//////////////////////////////////////////////////////////////////////
//  Symbolic finite-field solve for projective conics through the known
//  simple [6,6] point on the contact-6 core surface.
//
//  For each tangent direction d over GF(p), introduce unknowns
//
//      e1,e2,e3,e4,w1,w2
//
//  in
//
//      X_i(t) = P_i + (d_i*t + e_i*t^2)/(1 + w1*t + w2*t^2).
//
//  Clearing denominators in G2=G3 gives coefficient equations in these six
//  unknowns.  This avoids brute-force enumeration over p^6 choices.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned p0 then
    p0 := 7;
elif Type(p0) eq MonStgElt then
    p0 := StringToInteger(p0);
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
    while core ne 0 do
        g := GCD(core, boundary);
        if TotalDegree(g) eq 0 then
            break;
        end if;
        core := ExactQuotient(core, g);
    end while;
    return core;
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
    G2core := RepeatedBoundaryRemoval(G2red, boundary);
    G3core := RepeatedBoundaryRemoval(G3red, boundary);
    return G2core, G3core;
end function;

function IsNormalizedNonzero(vec)
    if vec eq [Parent(vec[1])!0 : i in [1..#vec]] then
        return false;
    end if;
    first := Min([i : i in [1..#vec] | vec[i] ne 0]);
    return vec[first] eq 1;
end function;

function Dot(a, b)
    F := Parent(a[1]);
    return &+[F!(a[i]*b[i]) : i in [1..#a]];
end function;

function TangentLinearForms(G2p, G3p, pt)
    grads2 := [Evaluate(Derivative(G2p, i), <pt[1],pt[2],pt[3],pt[4]>)
               : i in [1..4]];
    grads3 := [Evaluate(Derivative(G3p, i), <pt[1],pt[2],pt[3],pt[4]>)
               : i in [1..4]];
    return grads2, grads3;
end function;

function HomogeneousEvaluate(poly, X, W)
    deg := TotalDegree(poly);
    mons := Monomials(poly);
    coeffs := Coefficients(poly);
    R := Parent(W);
    value := R!0;
    for k in [1..#mons] do
        exps := Exponents(mons[k]);
        term := R!coeffs[k];
        edeg := 0;
        for i in [1..4] do
            if exps[i] ne 0 then
                term *:= X[i]^exps[i];
                edeg +:= exps[i];
            end if;
        end for;
        value +:= term*W^(deg-edeg);
    end for;
    return value;
end function;

G2, G3 := BuildGenericCore();
F := GF(p0);
S4<b,L,U,v> := PolynomialRing(F, 4);
phi := hom<Parent(G2) -> S4 | b,L,U,v>;
G2p := phi(G2);
G3p := phi(G3);
pt := [F!(-7)/F!13, F!29/F!16, -F!9/F!4, F!5/F!2];
assert Evaluate(G2p, <pt[1],pt[2],pt[3],pt[4]>) eq 0;
assert Evaluate(G3p, <pt[1],pt[2],pt[3],pt[4]>) eq 0;

grad2, grad3 := TangentLinearForms(G2p, G3p, pt);
tangent_dirs := [];
for d1 in F do
    for d2 in F do
        for d3 in F do
            for d4 in F do
                d := [d1,d2,d3,d4];
                if not IsNormalizedNonzero(d) then
                    continue;
                end if;
                if Dot(grad2,d) eq 0 and Dot(grad3,d) eq 0 then
                    Append(~tangent_dirs, d);
                end if;
            end for;
        end for;
    end for;
end for;

print "Contact-6 core projective conic coefficient solve";
print "p", p0, "tangent_directions", #tangent_dirs;

total_solutions := 0;
positive_dim := 0;
samples := [];

for idx in [1..#tangent_dirs] do
    d := tangent_dirs[idx];
    R<e1,e2,e3,e4,w1,w2> := PolynomialRing(F, 6);
    T<t> := PolynomialRing(R);
    W := T!1 + T!w1*t + T!w2*t^2;
    es := [e1,e2,e3,e4];
    X := [T!pt[i]*W + T!d[i]*t + T!es[i]*t^2 : i in [1..4]];
    H2 := HomogeneousEvaluate(G2p, X, W);
    H3 := HomogeneousEvaluate(G3p, X, W);
    coeffs := [R!c : c in Coefficients(H2)] cat [R!c : c in Coefficients(H3)];
    coeffs := [c : c in coeffs | c ne 0];
    I := ideal<R | coeffs>;
    try
        dim := Dimension(I);
        if dim gt 0 then
            positive_dim +:= 1;
            print " direction", idx, "d", [Z!d[i] : i in [1..4]],
                  "positive_dimension", dim,
                  "basis_degrees", [TotalDegree(g) : g in Basis(I)];
            continue;
        end if;
    catch e
        print " direction", idx, "dimension unavailable", e`Object;
    end try;
    try
        V := Variety(I);
        if #V gt 0 then
            total_solutions +:= #V;
            print " direction", idx, "d", [Z!d[i] : i in [1..4]], "solutions", #V;
            for P in V do
                if #samples lt 5 then
                    Append(~samples, <[Z!d[i] : i in [1..4]],
                                      [Z!P[i] : i in [1..4]],
                                      Z!P[5], Z!P[6]>);
                end if;
            end for;
        end if;
    catch e
        print " direction", idx, "variety unavailable", e`Object,
              "basis_degrees", [TotalDegree(g) : g in Basis(I)];
    end try;
end for;

print "summary p", p0,
      "tangent_directions", #tangent_dirs,
      "positive_dim_directions", positive_dim,
      "zero_dim_solutions", total_solutions,
      "samples", samples;

quit;
