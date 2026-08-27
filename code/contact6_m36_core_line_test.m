//////////////////////////////////////////////////////////////////////
//  Test for lines through the known simple [6,6] point on the generic
//  contact-6 core-cover surface.
//
//  If a rational line through the point existed, its reduction mod every
//  good prime would give a line direction over GF(p).  Thus a single good
//  prime with no direction rules out this simplest route to an infinite
//  family through the known point.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned prime_bound then
    prime_bound := 31;
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

function DirectionWorks(G2p, G3p, pt, dir)
    F<t> := PolynomialRing(Parent(pt[1]));
    vals := <F!pt[i] + F!dir[i]*t : i in [1..4]>;
    return Evaluate(G2p, vals) eq 0 and Evaluate(G3p, vals) eq 0;
end function;

G2, G3 := BuildGenericCore();

print "Contact-6 core line test through known simple [6,6] point";
for p in [q : q in PrimesUpTo(prime_bound) | q notin {2,3,13}] do
    F := GF(p);
    S4<b,L,U,v> := PolynomialRing(F, 4);
    phi := hom<Parent(G2) -> S4 | b,L,U,v>;
    G2p := phi(G2);
    G3p := phi(G3);
    pt := [F!(-7)/F!13, F!29/F!16, -F!9/F!4, F!5/F!2];
    assert Evaluate(G2p, <pt[1],pt[2],pt[3],pt[4]>) eq 0;
    assert Evaluate(G3p, <pt[1],pt[2],pt[3],pt[4]>) eq 0;

    directions := 0;
    samples := [];
    for d1 in F do
        for d2 in F do
            for d3 in F do
                for d4 in F do
                    dir := [d1,d2,d3,d4];
                    if dir eq [F!0,F!0,F!0,F!0] then
                        continue;
                    end if;
                    first := Min([i : i in [1..4] | dir[i] ne 0]);
                    if dir[first] ne 1 then
                        continue;
                    end if;
                    if DirectionWorks(G2p, G3p, pt, dir) then
                        directions +:= 1;
                        if #samples lt 5 then
                            Append(~samples, <Z!d1,Z!d2,Z!d3,Z!d4>);
                        end if;
                    end if;
                end for;
            end for;
        end for;
    end for;
    print "p", p, "line_directions", directions, "samples", samples;
end for;

quit;
