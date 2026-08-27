//////////////////////////////////////////////////////////////////////
//  Slice solver for the contact-6 [1,2,2] core cover.
//
//  Fix (b,v), eliminate a from F1, and solve the remaining two equations
//  in (M,U), where M=L^2.  Keep only rational-square M and verify the
//  resulting curve/divisors exactly.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned height then
    height := 6;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;
if not assigned progress_interval then
    progress_interval := 50;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;
if not assigned simple_only then
    simple_only := true;
elif Type(simple_only) eq MonStgElt then
    simple_only := simple_only in {"true", "True", "1", "yes"};
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

function FactorDegrees(f)
    return Sort([Degree(fe[1]) : fe in Factorization(f)]);
end function;

function HasFactorType122(f)
    return FactorDegrees(f) eq [1,2,2];
end function;

function Has66(invs)
    return #[n : n in invs | (Z!n) mod 6 eq 0] ge 2;
end function;

function IntegralModel(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
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
    return F1,F2,F3,q,B,Delta,f,h6;
end function;

function VerifyCandidate(a,b,L,U,v)
    M := L^2;
    F1,F2,F3,q,B,Delta,f,h6 := CubicContactDataM(a,b,M,U,v);
    if F1 ne 0 or F2 ne 0 or F3 ne 0 then
        return false, [], 0, 0, false, 0, f;
    end if;
    if not GoodPolynomial(f) or Evaluate(h6, Q!1) eq 0 then
        return false, [], 0, 0, false, 0, f;
    end if;
    if not HasFactorType122(f) then
        return false, [], 0, 0, false, 0, f;
    end if;
    if Discriminant(q) eq 0 or Degree(GCD(q,f)) gt 0 then
        return false, [], 0, 0, false, 0, f;
    end if;
    simple, pcert, Lp := SimpleCertificate(f);
    if simple_only and not simple then
        return false, [], 0, 0, simple, pcert, f;
    end if;
    h3 := (1/L)*x^3 + (B/(2*L))*x^2 + (Delta/(8*L))*x + v^3/L;
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    D := J![x-1, Evaluate(h6, Q!1)];
    E := J![q, h3 mod q];
    ordD := Order(D);
    ordE := Order(E);
    fI, scale := IntegralModel(f);
    CI := HyperellipticCurve(fI);
    JI := Jacobian(CI);
    G, phi := TorsionSubgroup(JI);
    return true, Invariants(G), ordD, ordE, simple, pcert, fI;
end function;

function PrimitivePolynomial(f)
    if f eq 0 then
        return f;
    end if;
    den := LCM([Denominator(c) : c in Coefficients(f)]);
    g := Parent(f)!(den*f);
    nums := [Z!c : c in Coefficients(g)];
    cont := GCD([Abs(n) : n in nums | n ne 0]);
    if cont gt 1 then
        g := Parent(f)!(g/cont);
    end if;
    return g;
end function;

function GenericSliceSolutions(b0, v0)
    R<M,U> := PolynomialRing(Q, 2);
    K := FieldOfFractions(R);
    PA<a> := PolynomialRing(K);
    b := K!b0;
    v := K!v0;

    c1 := 2*a + 6;
    c2 := a^2 + 2*b - 15;
    c3 := 2*a*b + 22;
    c4 := 2*a + b^2 - 15;
    c5 := 2*b + 6;
    B := c5*M + 3*U;
    Delta := 4*c4*M + 12*(U^2 + v^2) - B^2;
    F3 := B*Delta + 16*v^3 - 8*c3*M - 8*U^3 - 48*U*v^2;
    F2 := Delta^2 + 64*B*v^3 - 64*c2*M - 192*(U^2*v^2 + v^4);
    F1 := Delta*v^3 - 4*c1*M - 12*U*v^4;
    D := Coefficient(F1, 1);
    N := -Coefficient(F1, 0);
    if D eq 0 then
        return [];
    end if;

    function SubNum(P)
        d := Degree(P);
        value := K!0;
        for i in [0..d] do
            value +:= Coefficient(P,i)*N^i*D^(d-i);
        end for;
        return PrimitivePolynomial(R!Numerator(value));
    end function;

    G2 := SubNum(F2);
    G3 := SubNum(F3);
    boundary := M*(U^2 - 4*v0^2)*(R!(N + (b0+2)*D));
    I := ideal<R | G2, G3>;
    try
        I := Saturation(I, ideal<R | boundary>);
    catch e
        return [];
    end try;
    if Dimension(I) ne 0 then
        return [];
    end if;
    sols := [];
    for pt in Variety(I) do
        M0 := Q!pt[1];
        U0 := Q!pt[2];
        if M0 eq 0 then
            continue;
        end if;
        ok, L0 := IsSquare(M0);
        if not ok then
            continue;
        end if;
        Nnum := R!Numerator(N);
        Nden := R!Denominator(N);
        Dnum := R!Numerator(D);
        Dden := R!Denominator(D);
        a0 := (Evaluate(Nnum, <M0,U0>)/Evaluate(Nden, <M0,U0>)) /
              (Evaluate(Dnum, <M0,U0>)/Evaluate(Dden, <M0,U0>));
        Append(~sols, <a0,b0,L0,U0,v0>);
        if L0 ne 0 then
            Append(~sols, <a0,b0,-L0,U0,v0>);
        end if;
    end for;
    return sols;
end function;

params := RationalParametersOfHeight(height);
checked := 0;
sols_total := 0;
verified := 0;
hits := [];

print "Contact-6 core slice search";
print "height", height, "parameter_count", #params, "simple_only", simple_only;

for b in params do
    for v in params do
        if v eq 0 or v^3 eq 1 then
            continue;
        end if;
        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "progress", checked, "slice_solutions", sols_total,
                  "verified", verified, "hits", #hits;
        end if;
        sols := GenericSliceSolutions(b,v);
        sols_total +:= #sols;
        for S in sols do
            a := S[1]; b0 := S[2]; L := S[3]; U := S[4]; v0 := S[5];
            ok, invs, ordD, ordE, simple, pcert, f :=
                VerifyCandidate(a,b0,L,U,v0);
            if ok then
                verified +:= 1;
            end if;
            if ok and ordD eq 6 and ordE eq 3 and Has66(invs) then
                Append(~hits, <a,b0,L,U,v0,invs,ordD,ordE,simple,pcert,f>);
                print "HIT66", "a", a, "b", b0, "L", L, "U", U, "v", v0,
                      "invs", invs, "simple", simple, "pcert", pcert;
                print " f", f;
                if #hits ge max_hits then
                    break b;
                end if;
            end if;
        end for;
    end for;
end for;

print "Done";
print "checked_slices", checked;
print "slice_solutions", sols_total;
print "verified", verified;
print "hits", #hits;
for H in hits do
    print "H", H;
end for;

quit;
