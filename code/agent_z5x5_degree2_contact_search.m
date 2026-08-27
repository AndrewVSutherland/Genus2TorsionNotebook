//////////////////////////////////////////////////////////////////////
//  Smoke tests for the degree-2/contact-5 route.
//
//  This uses the normal form from agent_z5x5_degree2_contact_symbolic.m:
//
//      U = x^2+s*x+t,
//      H = x^5 + ... + m,
//      f = H^2 - U^5.
//
//  If the quartic tail of f is h^2, then
//
//      f = h^2 - K*x^5
//
//  gives a rational contact-5 class D0, while y-H gives
//  5*D2=0 for D2=[U,H mod U].  This script verifies those classes
//  over Q on the fastest rational branch, and over small finite fields.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned prime_bound then
    prime_bound := 23;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;

if not assigned full_field_search then
    full_field_search := true;
elif Type(full_field_search) eq MonStgElt then
    full_field_search := (full_field_search eq "true") or (full_field_search eq "True");
end if;

if not assigned rational_height then
    rational_height := 3;
elif Type(rational_height) eq MonStgElt then
    rational_height := StringToInteger(rational_height);
end if;

function PolynomialSquareRoot(Qpoly)
    P := Parent(Qpoly);
    if Qpoly eq 0 then
        return true, P!0;
    end if;

    ok, lcroot := IsSquare(LeadingCoefficient(Qpoly));
    if not ok then
        return false, P!0;
    end if;

    h := P!lcroot;
    for pair in Factorization(Qpoly) do
        if pair[2] mod 2 ne 0 then
            return false, P!0;
        end if;
        h *:= pair[1]^(pair[2] div 2);
    end for;

    if h^2 ne Qpoly then
        h := -h;
    end if;
    if h^2 ne Qpoly then
        return false, P!0;
    end if;
    return true, h;
end function;

function Degree2Model(k, s0, t0, m0)
    P<x> := PolynomialRing(k);

    U := x^2 + s0*x + t0;
    U5 := U^5;
    a4 := Coefficient(U5, 9)/2;
    a3 := (Coefficient(U5, 8) - a4^2)/2;
    a2 := (Coefficient(U5, 7) - 2*a4*a3)/2;
    a1 := (Coefficient(U5, 6) - 2*a4*a2 - a3^2)/2;
    H := x^5 + a4*x^4 + a3*x^3 + a2*x^2 + a1*x + m0;
    f := H^2 - U5;

    if Degree(f) gt 5 then
        return false, P!0, P!0, P!0, P!0, "degree(f)>5";
    end if;
    if Degree(f) ne 5 then
        return false, P!0, P!0, P!0, P!0, "degree(f)!=5";
    end if;

    tail := &+[Coefficient(f, i)*x^i : i in [0..4]];
    ok, h := PolynomialSquareRoot(tail);
    if not ok then
        return false, P!0, P!0, P!0, P!0, "quartic tail not square";
    end if;

    if Discriminant(f) eq 0 then
        return false, P!0, P!0, P!0, P!0, "singular";
    end if;

    if Discriminant(U) eq 0 then
        return false, P!0, P!0, P!0, P!0, "U not squarefree";
    end if;

    return true, f, h, U, H mod U, "ok";
end function;

function RelationPairs(J, D0, D2)
    rels := [];
    for a in [0..4] do
        for b in [0..4] do
            if a eq 0 and b eq 0 then
                continue;
            end if;
            if a*D0 + b*D2 eq J!0 then
                Append(~rels, <a,b>);
            end if;
        end for;
    end for;
    return rels;
end function;

function CheckClasses(f, h, U, V)
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    P<x> := Parent(f);

    D0 := J![x, Evaluate(h, 0)];
    D2 := J![U, V];

    if D0 eq J!0 then
        return false, J!0, J!0, [], "D0=0";
    end if;
    if D2 eq J!0 then
        return false, J!0, J!0, [], "D2=0";
    end if;
    if 5*D0 ne J!0 then
        return false, J!0, J!0, [], "5*D0!=0";
    end if;
    if 5*D2 ne J!0 then
        return false, J!0, J!0, [], "5*D2!=0";
    end if;

    rels := RelationPairs(J, D0, D2);
    return (#rels eq 0), D0, D2, rels, "checked";
end function;

procedure RationalBranchSmoke()
    Q := Rationals();
    samples := [Q!1, Q!-1, Q!2, Q!-2];

    print "# Rational branch t=5*s^2/4, m=45*s^5/32";
    for s0 in samples do
        t0 := 5*s0^2/4;
        m0 := 45*s0^5/32;
        ok, f, h, U, V, reason := Degree2Model(Q, s0, t0, m0);
        printf "s=%o t=%o m=%o model=%o", s0, t0, m0, reason;
        if not ok then
            print "";
            continue;
        end if;

        independent, D0, D2, rels, creason := CheckClasses(f, h, U, V);
        printf " class_check=%o independent=%o rels=%o\n", creason, independent, rels;
        printf "  f=%o\n", f;
        printf "  h=%o\n", h;
        printf "  U=%o\n", U;
        printf "  V=H mod U=%o\n", V;
    end for;

    print "";
    print "# Boundary branch t=s^2/4 from the resultant";
    for s0 in [Q!1, Q!-1] do
        for m0 in [s0^5/32, 3*s0^5/64] do
            t0 := s0^2/4;
            ok, f, h, U, V, reason := Degree2Model(Q, s0, t0, m0);
            printf "s=%o t=%o m=%o model=%o\n", s0, t0, m0, reason;
        end for;
    end for;
end procedure;

procedure SmallRationalAffineSmoke()
    Q := Rationals();
    tested := 0;
    square_smooth := 0;
    independent_hits := 0;

    print "";
    printf "# Small affine rational smoke: integer s,t,m in [-%o,%o]\n",
        rational_height, rational_height;
    for si in [-rational_height..rational_height] do
        for ti in [-rational_height..rational_height] do
            for mi in [-rational_height..rational_height] do
                tested +:= 1;
                ok, f, h, U, V, reason := Degree2Model(Q, Q!si, Q!ti, Q!mi);
                if not ok then
                    continue;
                end if;
                square_smooth +:= 1;
                independent, D0, D2, rels, creason := CheckClasses(f, h, U, V);
                if independent then
                    independent_hits +:= 1;
                    printf "rational_full_hit s=%o t=%o m=%o\n", si, ti, mi;
                    printf "  f=%o\n  h=%o\n  U=%o\n  V=%o\n", f, h, U, V;
                else
                    printf "rational_model s=%o t=%o m=%o rels=%o\n", si, ti, mi, rels;
                end if;
            end for;
        end for;
    end for;
    printf "tested=%o square_smooth=%o independent_hits=%o\n",
        tested, square_smooth, independent_hits;
end procedure;

procedure BranchFiniteFieldSmoke()
    print "";
    print "# Finite-field smoke on branch t=5*s^2/4, m=45*s^5/32";
    for p in PrimesUpTo(prime_bound) do
        if p in {2,5} then
            continue;
        end if;

        Fp := GF(p);
        hits := 0;
        for s0 in Fp do
            if s0 eq 0 then
                continue;
            end if;
            t0 := 5*s0^2/4;
            m0 := 45*s0^5/32;
            ok, f, h, U, V, reason := Degree2Model(Fp, s0, t0, m0);
            if not ok then
                continue;
            end if;
            independent, D0, D2, rels, creason := CheckClasses(f, h, U, V);
            if independent then
                hits +:= 1;
                printf "p=%o branch_hit s=%o #J=%o f=%o U=%o V=%o\n",
                    p, Integers()!s0, #Jacobian(HyperellipticCurve(f)), f, U, V;
                break;
            end if;
        end for;
        if hits eq 0 then
            printf "p=%o branch_no_independent_hit\n", p;
        end if;
    end for;
end procedure;

procedure FullFiniteFieldSmoke()
    if not full_field_search then
        return;
    end if;

    print "";
    print "# Full finite-field scan of small (s,t,m)";
    for p in PrimesUpTo(prime_bound) do
        if p in {2,5} then
            continue;
        end if;

        Fp := GF(p);
        tested := 0;
        square_models := 0;
        smooth_models := 0;
        independent_hits := 0;

        for s0 in Fp do
            if independent_hits gt 0 then
                break;
            end if;
            for t0 in Fp do
                if independent_hits gt 0 then
                    break;
                end if;
                for m0 in Fp do
                    tested +:= 1;
                    ok, f, h, U, V, reason := Degree2Model(Fp, s0, t0, m0);
                    if not ok then
                        continue;
                    end if;
                    square_models +:= 1;
                    smooth_models +:= 1;

                    independent, D0, D2, rels, creason := CheckClasses(f, h, U, V);
                    if independent then
                        independent_hits +:= 1;
                        printf "p=%o full_hit s=%o t=%o m=%o #J=%o\n",
                            p, Integers()!s0, Integers()!t0, Integers()!m0,
                            #Jacobian(HyperellipticCurve(f));
                        printf "  f=%o\n  h=%o\n  U=%o\n  V=%o\n", f, h, U, V;
                        break;
                    end if;
                end for;
            end for;
        end for;

        printf "p=%o tested=%o square_smooth=%o independent_hits=%o\n",
            p, tested, smooth_models, independent_hits;
    end for;
end procedure;

RationalBranchSmoke();
SmallRationalAffineSmoke();
BranchFiniteFieldSmoke();
FullFiniteFieldSmoke();

quit;
