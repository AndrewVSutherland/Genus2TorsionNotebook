//////////////////////////////////////////////////////////////////////
//  Finite check: on eps=+1, do the extra-root/cubic-contact points satisfy
//  one of the explicit extra-involution conditions?
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned prime_bound then
    prime_bound := 19;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;

Z := Integers();

function ExtraRootA(F, eps, r, b)
    return (eps*(r-1)^3 - 1 - r^3 - b*r^2)/r;
end function;

function ContactData(F, eps, r, b, L, U, v)
    PF<x> := PolynomialRing(F);
    a := ExtraRootA(F, eps, r, b);
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
    return a, F1, F2, F3, f, h6, q;
end function;

function QuadraticRoots(F, y0, y1, y2)
    PB<z> := PolynomialRing(F);
    c2 := (y2 - 2*y1 + y0)/2;
    c1 := y1 - y0 - c2;
    c0 := y0;
    poly := PB![c0, c1, c2];
    if poly eq 0 then
        return [], true;
    end if;
    return [rt[1] : rt in Roots(poly)], false;
end function;

function IsGood(F, eps, r, b, L, U, v)
    a,F1,F2,F3,f,h6,q := ContactData(F, eps, r, b, L, U, v);
    if F1 ne 0 or F2 ne 0 or F3 ne 0 then
        return false, [], a, f;
    end if;
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false, [], a, f;
    end if;
    if Evaluate(h6, F!1) eq 0 or Evaluate(f, r) ne 0 then
        return false, [], a, f;
    end if;
    if Discriminant(q) eq 0 or Degree(GCD(q, f)) gt 0 then
        return false, [], a, f;
    end if;
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    A, phi := AbelianGroup(J);
    return true, Invariants(A), a, f;
end function;

function AutoA(b,r)
    return [
        2*b^3*r^4 + 30*b^2*r^4 - 36*b^2*r^3 + 8*b^2*r^2 + 126*b*r^4 - 296*b*r^3 + 216*b*r^2 - 48*b*r + 162*r^4 - 612*r^3 + 816*r^2 - 464*r + 96,
        -b^4*r^4 - 12*b^3*r^4 + 12*b^3*r^3 - 2*b^3*r^2 - 54*b^2*r^4 + 108*b^2*r^3 - 66*b^2*r^2 + 12*b^2*r - 108*b*r^4 + 324*b*r^3 - 342*b*r^2 + 152*b*r - 24*b - 81*r^4 + 324*r^3 - 470*r^2 + 300*r - 72,
        b^4*r^3 + 6*b^3*r^3 - 6*b^3*r^2 - 14*b^2*r^2 + 12*b^2*r - 54*b*r^3 + 102*b*r^2 - 48*b*r - 81*r^3 + 270*r^2 - 284*r + 96
    ];
end function;

function AutoB(b,r)
    return [
        12*b^3*r^4 - 12*b^3*r^3 + 4*b^3*r^2 + 108*b^2*r^4 - 188*b^2*r^3 + 120*b^2*r^2 - 36*b^2*r + 4*b^2 + 324*b*r^4 - 804*b*r^3 + 732*b*r^2 - 296*b*r + 48*b + 324*r^4 - 1044*r^3 + 1224*r^2 - 612*r + 108,
        -6*b^3*r^4 + 6*b^3*r^3 - 2*b^3*r^2 - 54*b^2*r^4 + 94*b^2*r^3 - 66*b^2*r^2 + 24*b^2*r - 4*b^2 - 162*b*r^4 + 378*b*r^3 - 342*b*r^2 + 144*b*r - 24*b - 162*r^4 + 450*r^3 - 470*r^2 + 216*r - 36,
        6*b^2*r^2 - 6*b^2*r + 2*b^2 + 24*b*r^3 - 24*b*r^2 + 4*b*r + 72*r^3 - 142*r^2 + 90*r - 18
    ];
end function;

function IsAuto(b,r)
    A := AutoA(b,r);
    B := AutoB(b,r);
    return (&and[e eq 0 : e in A]) or (&and[e eq 0 : e in B]);
end function;

print "Finite automorphism-condition check on eps=+1 extra-root branch";
print "prime_bound", prime_bound;
for p in [p : p in PrimesUpTo(prime_bound) | p notin {2,3}] do
    F := GF(p);
    good := 0;
    auto := 0;
    nonauto_samples := [];
    for r in F do
        if r eq 0 or r eq 1 then
            continue;
        end if;
        eps := F!1;
        for L in F do
            if L eq 0 then
                continue;
            end if;
            for U in F do
                for v in F do
                    if v eq 0 or U^2 - 4*v^2 eq 0 then
                        continue;
                    end if;
                    a0,y0,F20,F30,f0,h0,q0 := ContactData(F,eps,r,F!0,L,U,v);
                    a1,y1,F21,F31,f1,h1,q1 := ContactData(F,eps,r,F!1,L,U,v);
                    a2,y2,F22,F32,f2,h2,q2 := ContactData(F,eps,r,F!2,L,U,v);
                    roots, iszero := QuadraticRoots(F,y0,y1,y2);
                    if iszero then
                        continue;
                    end if;
                    for b in roots do
                        ok, invs, aa, ff := IsGood(F,eps,r,b,L,U,v);
                        if not ok then
                            continue;
                        end if;
                        good +:= 1;
                        if IsAuto(b,r) then
                            auto +:= 1;
                        elif #nonauto_samples lt 5 then
                            Append(~nonauto_samples, <Z!r,Z!aa,Z!b,Z!L,Z!U,Z!v,invs>);
                        end if;
                    end for;
                end for;
            end for;
        end for;
    end for;
    print "p", p, "good_eps_plus", good, "auto_condition", auto,
          "nonauto", good-auto, "nonauto_samples", nonauto_samples;
end for;

quit;
