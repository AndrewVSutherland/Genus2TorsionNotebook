//////////////////////////////////////////////////////////////////////
//  Guided low-degree test: affine polynomial cubics through the known
//  simple [6,6] point and a second known surface point.
//
//  For P fixed and Q chosen from known exact [6,6] surface points, test
//
//      X(t) = P + d*t + e*t^2 + f*t^3,  f = Q-P-d-e,
//
//  so that X(0)=P and X(1)=Q.  The tangent d is restricted to the tangent
//  plane at P.  This is a guided degree-3 ansatz rather than a blind search
//  in all cubic coefficients.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned p0 then
    p0 := 7;
elif Type(p0) eq MonStgElt then
    p0 := StringToInteger(p0);
end if;

Qrat := Rationals();
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
    S<b,L,U,v> := PolynomialRing(Qrat, 4);
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

function CubicWorks(G2p, G3p, P, Q, d, e)
    F<t> := PolynomialRing(Parent(P[1]));
    f := [Q[i] - P[i] - d[i] - e[i] : i in [1..4]];
    vals := <F!P[i] + F!d[i]*t + F!e[i]*t^2 + F!f[i]*t^3 : i in [1..4]>;
    return Evaluate(G2p, vals) eq 0 and Evaluate(G3p, vals) eq 0;
end function;

function ReducePoint(F, pt)
    return [F!Numerator(pt[i])/F!Denominator(pt[i]) : i in [1..4]];
end function;

G2, G3 := BuildGenericCore();
F := GF(p0);
S4<b,L,U,v> := PolynomialRing(F, 4);
phi := hom<Parent(G2) -> S4 | b,L,U,v>;
G2p := phi(G2);
G3p := phi(G3);

Pq := [Qrat!(-7/13), Qrat!(29/16), Qrat!(-9/4), Qrat!(5/2)];
anchors := [
    <"same_simple_alt", [Qrat!(-7/13), Qrat!(13/4), Qrat!(-17/2), Qrat!(-5/4)]>,
    <"nonsimple_1a", [Qrat!(3/2), Qrat!(28/27), Qrat!(-64/27), Qrat!(-2/3)]>,
    <"nonsimple_1b", [Qrat!(3/2), Qrat!(10/27), Qrat!(-14/27), Qrat!(2/3)]>,
    <"nonsimple_2a", [Qrat!(1/8), Qrat!(91/125), Qrat!(-491/250), Qrat!(-4/5)]>,
    <"nonsimple_2b", [Qrat!(1/8), Qrat!(81/125), Qrat!(-329/250), Qrat!(4/5)]>,
    <"nonsimple_3a", [Qrat!(5/9), Qrat!(55/64), Qrat!(-199/96), Qrat!(-3/4)]>,
    <"nonsimple_3b", [Qrat!(5/9), Qrat!(35/64), Qrat!(-101/96), Qrat!(3/4)]>
];

P := ReducePoint(F, Pq);
assert Evaluate(G2p, <P[1],P[2],P[3],P[4]>) eq 0;
assert Evaluate(G3p, <P[1],P[2],P[3],P[4]>) eq 0;

grad2, grad3 := TangentLinearForms(G2p, G3p, P);
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

print "Contact-6 core guided affine cubic test";
print "p", p0, "tangent_directions", #tangent_dirs;

for A in anchors do
    label := A[1];
    Qp := ReducePoint(F, A[2]);
    if Evaluate(G2p, <Qp[1],Qp[2],Qp[3],Qp[4]>) ne 0
       or Evaluate(G3p, <Qp[1],Qp[2],Qp[3],Qp[4]>) ne 0 then
        print "anchor", label, "does_not_reduce_to_surface";
        continue;
    end if;
    cubics := 0;
    samples := [];
    for d in tangent_dirs do
        for e1 in F do
            for e2 in F do
                for e3 in F do
                    for e4 in F do
                        e := [e1,e2,e3,e4];
                        if CubicWorks(G2p, G3p, P, Qp, d, e) then
                            cubics +:= 1;
                            if #samples lt 3 then
                                Append(~samples,
                                       <[Z!d[i] : i in [1..4]],
                                        [Z!e[i] : i in [1..4]]>);
                            end if;
                        end if;
                    end for;
                end for;
            end for;
        end for;
    end for;
    print "anchor", label, "guided_cubics", cubics, "samples", samples;
end for;

quit;
