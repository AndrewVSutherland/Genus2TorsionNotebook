//////////////////////////////////////////////////////////////////////
//  Finite-prime conic ansatz on the contact-6 [6,6] core surface.
//
//  We test quadratic polynomial maps through the known simple point
//
//      X(t) = P + d*t + e*t^2,  X=(b,L,U,v),
//
//  on the generic eliminated core surface G2=G3.  The tangent vector d is
//  first constrained by the two linear tangent equations at P.  For each
//  tangent direction, all e are tried over GF(p).  A rational polynomial
//  conic through P with good reduction would survive this test modulo every
//  good prime.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned prime_bound then
    prime_bound := 23;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
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

function ConicWorks(G2p, G3p, pt, d, e)
    F<t> := PolynomialRing(Parent(pt[1]));
    vals := <F!pt[i] + F!d[i]*t + F!e[i]*t^2 : i in [1..4]>;
    return Evaluate(G2p, vals) eq 0 and Evaluate(G3p, vals) eq 0;
end function;

function TangentLinearForms(G2p, G3p, pt)
    grads2 := [Evaluate(Derivative(G2p, i), <pt[1],pt[2],pt[3],pt[4]>)
               : i in [1..4]];
    grads3 := [Evaluate(Derivative(G3p, i), <pt[1],pt[2],pt[3],pt[4]>)
               : i in [1..4]];
    return grads2, grads3;
end function;

function Dot(a, b)
    F := Parent(a[1]);
    return &+[F!(a[i]*b[i]) : i in [1..#a]];
end function;

G2, G3 := BuildGenericCore();

print "Contact-6 core polynomial conic test through known simple [6,6] point";
print "ansatz X(t)=P+d*t+e*t^2";

for p in [q : q in PrimesUpTo(prime_bound) | q notin {2,3,5,13}] do
    F := GF(p);
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

    conics := 0;
    samples := [];
    for d in tangent_dirs do
        for e1 in F do
            for e2 in F do
                for e3 in F do
                    for e4 in F do
                        e := [e1,e2,e3,e4];
                        if ConicWorks(G2p, G3p, pt, d, e) then
                            conics +:= 1;
                            if #samples lt 5 then
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

    print "p", p,
          "tangent_directions", #tangent_dirs,
          "polynomial_conics", conics,
          "samples", samples;
end for;

quit;
