//////////////////////////////////////////////////////////////////////
//  Finite-field simplicity diagnostic for the contact-6 extra-root
//  cubic-contact equations.
//
//  This tests whether the algebraic [2,6,6]-producing equations are wholly
//  contained in a decomposable/Humbert locus.  If a smooth good finite point
//  on the same equations has irreducible Frobenius quartic, then the full
//  finite variety is not contained in the decomposable locus modulo that p.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned prime_bound then
    prime_bound := 19;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;

Z := Integers();

function Has66(invs)
    return #[n : n in invs | (Z!n) mod 6 eq 0] ge 2;
end function;

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

function IsGoodContact(F, eps, r, b, L, U, v)
    a,F1,F2,F3,f,h6,q := ContactData(F, eps, r, b, L, U, v);
    if F1 ne 0 or F2 ne 0 or F3 ne 0 then
        return false, [], a, f, false, 0, f;
    end if;
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false, [], a, f, false, 0, f;
    end if;
    if Evaluate(h6, F!1) eq 0 or Evaluate(f, r) ne 0 then
        return false, [], a, f, false, 0, f;
    end if;
    if Discriminant(q) eq 0 or Degree(GCD(q, f)) gt 0 then
        return false, [], a, f, false, 0, f;
    end if;
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    A, phi := AbelianGroup(J);
    Lp := LPolynomial(C);
    fac := Factorization(Lp);
    irreducible := #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4;
    return true, Invariants(A), a, f, irreducible, Lp, f;
end function;

print "Finite simplicity diagnostic for contact-6 extra-root cubic-contact equations";
print "prime_bound", prime_bound;

for p in [p : p in PrimesUpTo(prime_bound) | p notin {2,3}] do
    F := GF(p);
    contact := 0;
    good := 0;
    finite66 := 0;
    irreducible_good := 0;
    irreducible_66 := 0;
    reducible_66 := 0;
    sample_irred := [];
    sample_red66 := [];

    for eps in [F!1, -F!1] do
        for r in F do
            if r eq 0 or r eq 1 then
                continue;
            end if;
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
                            a,F1,F2,F3,f,h6,q := ContactData(F,eps,r,b,L,U,v);
                            if F1 ne 0 or F2 ne 0 or F3 ne 0 then
                                continue;
                            end if;
                            contact +:= 1;
                            ok, invs, aa, ff, irred, Lp, fpoly :=
                                IsGoodContact(F,eps,r,b,L,U,v);
                            if not ok then
                                continue;
                            end if;
                            good +:= 1;
                            if irred then
                                irreducible_good +:= 1;
                                if #sample_irred lt 3 then
                                    Append(~sample_irred,
                                           <Z!eps,Z!r,Z!aa,Z!b,Z!L,Z!U,Z!v,invs,Lp>);
                                end if;
                            end if;
                            if Has66(invs) then
                                finite66 +:= 1;
                                if irred then
                                    irreducible_66 +:= 1;
                                else
                                    reducible_66 +:= 1;
                                    if #sample_red66 lt 3 then
                                        Append(~sample_red66,
                                               <Z!eps,Z!r,Z!aa,Z!b,Z!L,Z!U,Z!v,invs,Lp>);
                                    end if;
                                end if;
                            end if;
                        end for;
                    end for;
                end for;
            end for;
        end for;
    end for;

    print "p", p,
          "contact", contact,
          "good", good,
          "finite66", finite66,
          "irreducible_good", irreducible_good,
          "irreducible_66", irreducible_66,
          "reducible_66", reducible_66;
    print " sample_irred", sample_irred;
    print " sample_red66", sample_red66;
end for;

quit;
