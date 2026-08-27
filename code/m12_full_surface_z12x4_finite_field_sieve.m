//////////////////////////////////////////////////////////////////////
//  Finite-field sieve on the full two-dimensional extra-Weierstrass
//  surface inside M(12).
//
//  We split T+1 by writing
//      a = (1-z^2)/(4*(r+1)),
//  and then ask for a rational root u of Q4 = W/(T+1).  Such a root
//  gives an independent rational 2-torsion class.  This script tests,
//  over finite fields, whether any such class is divisible by 2.
//
//  If no good F_p point works, any rational Z/12 x Z/4 example on this
//  full surface must reduce to the boundary modulo p.
//
//  Typical run from torsion_jac:
//      magma code/m12_full_surface_z12x4_finite_field_sieve.m
//////////////////////////////////////////////////////////////////////

prime_list := [5,7,11,13,17,19,23,29,31,37,41,43];

function IsDivisibleBy2Finite(J, D)
    G, phi := AbelianGroup(J);
    a := D @@ phi;
    coords := Eltseq(a);
    invs := Invariants(G);
    for i in [1..#coords] do
        if GCD(2, invs[i]) ne 1 and coords[i] mod 2 ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

for p in prime_list do
    F := GF(p);
    P<x> := PolynomialRing(F);
    PX<X> := PolynomialRing(F);

    good := 0;
    qroots := 0;
    classes := 0;
    bad := 0;
    divs := [];

    for r in F do
        if r eq -1 then
            continue;
        end if;
        for z in F do
            if z^2 eq 1 then
                continue;
            end if;

            a := (1-z^2)/(4*(r+1));
            if a eq 0 then
                continue;
            end if;

            T := a*x^2 - x + r;
            h := (x-r)*(T+1);
            W := h^2 + 4*a*x^2*T*(T+1);
            if Degree(W) ne 6 or Discriminant(W) eq 0 then
                bad +:= 1;
                continue;
            end if;

            rootsT := Roots(T+1);
            if #rootsT lt 2 then
                bad +:= 1;
                continue;
            end if;

            Q4 := ExactQuotient(W, T+1);
            rootsQ := [ rt[1] : rt in Roots(Q4) | rt[2] eq 1 ];
            if #rootsQ eq 0 then
                continue;
            end if;

            good +:= 1;
            qroots +:= #rootsQ;

            for wd in rootsT do
                w := wd[1];
                f5 := PX!0;
                for i in [0..Degree(W)] do
                    for j in [0..i] do
                        f5 +:= Coefficient(W, i)*Binomial(i,j)*w^(i-j)*X^(6-j);
                    end for;
                end for;
                if Degree(f5) ne 5 or Discriminant(f5) eq 0 then
                    continue;
                end if;

                C := HyperellipticCurve(f5);
                J := Jacobian(C);
                for u in rootsQ do
                    if u eq w then
                        continue;
                    end if;
                    beta := 1/(u-w);
                    Tbeta := J![X-beta, F!0];
                    if Tbeta eq J!0 then
                        continue;
                    end if;
                    classes +:= 1;
                    if IsDivisibleBy2Finite(J, Tbeta) then
                        Append(~divs, <r,z,w,u,beta>);
                        if #divs le 10 then
                            print "  divisible", p, <r,z,w,u,beta>;
                        end if;
                    end if;
                end for;
            end for;
        end for;
    end for;

    print "p", p, "good_rz_with_qroot", good, "qroots", qroots,
          "classes", classes, "divisible", #divs, "bad", bad;
    if classes gt 0 and #divs eq 0 then
        print "NO GOOD FULL-SURFACE HALVES MOD", p;
        break;
    end if;
end for;
quit;
