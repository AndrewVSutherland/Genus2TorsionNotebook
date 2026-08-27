//////////////////////////////////////////////////////////////////////
//  Full degree-2 Mumford norm probe for the Z/5 x Z/5 lane.
//
//  One contact-5 class is fixed at x=0:
//
//      f = h^2 - K*x^5,       h = 1 + h1*x + h2*x^2.
//
//  For a second degree-2 class, use the full norm identity
//
//      A^2 - B^2*f = Lambda*U^5,
//      U = x^2 + s*x + t,     deg(A)<=5, deg(B)<=2.
//
//  On the generic chart A monic, Lambda=1.  The coefficients a4..a0 of
//  A are forced by canceling degrees 9..5.  The script records the five
//  remaining residual equations and runs bounded rational / finite-field
//  probes with B nonconstant, so it does not revisit the B=1 slice.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned do_symbolic then
    do_symbolic := true;
elif Type(do_symbolic) eq MonStgElt then
    do_symbolic := (do_symbolic eq "true") or (do_symbolic eq "True");
end if;

if not assigned rational_height then
    rational_height := 1;
elif Type(rational_height) eq MonStgElt then
    rational_height := StringToInteger(rational_height);
end if;

if not assigned prime_bound then
    prime_bound := 31;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;

if not assigned finite_trials then
    finite_trials := 300000;
elif Type(finite_trials) eq MonStgElt then
    finite_trials := StringToInteger(finite_trials);
end if;

if not assigned cantor_trials then
    cantor_trials := 20000;
elif Type(cantor_trials) eq MonStgElt then
    cantor_trials := StringToInteger(cantor_trials);
end if;

if not assigned random_seed then
    random_seed := 20260702;
elif Type(random_seed) eq MonStgElt then
    random_seed := StringToInteger(random_seed);
end if;
SetSeed(random_seed);

function ForcedAData(h, kap, U, B)
    P := Parent(h);
    x := P.1;
    f := h^2 - kap*x^5;
    A := x^5;
    E := A^2 - B^2*f - U^5;
    avec := [];

    for deg in [9,8,7,6,5] do
        a := -Coefficient(E, deg)/2;
        Append(~avec, a);
        A +:= a*x^(deg-5);
        E := A^2 - B^2*f - U^5;
    end for;

    return A, E, f, avec;
end function;

function Residuals(E)
    return [Coefficient(E, i) : i in [0..4]];
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
    P := Parent(f);
    x := P.1;
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    D0 := J![x, Evaluate(h, 0)];
    D2 := J![U, V];

    if D0 eq J!0 or D2 eq J!0 then
        return false, [], "zero class";
    end if;
    if 5*D0 ne J!0 then
        return false, [], "5*D0 != 0";
    end if;
    if 5*D2 ne J!0 then
        return false, [], "5*D2 != 0";
    end if;

    rels := RelationPairs(J, D0, D2);
    return (#rels eq 0), rels, "checked";
end function;

function FullNormModel(k, h1v, h2v, kapv, sv, tv, b0v, b1v, b2v)
    P<x> := PolynomialRing(k);
    h := 1 + h1v*x + h2v*x^2;
    U := x^2 + sv*x + tv;
    B := b0v + b1v*x + b2v*x^2;
    A, E, f, avec := ForcedAData(h, kapv, U, B);

    if &or[Coefficient(E, i) ne 0 : i in [0..4]] then
        return false, P!0, P!0, P!0, P!0, P!0, [], "residual nonzero";
    end if;
    if kapv eq 0 or Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false, P!0, P!0, P!0, P!0, P!0, [], "bad contact curve";
    end if;
    if b1v eq 0 and b2v eq 0 then
        return false, P!0, P!0, P!0, P!0, P!0, [], "constant B slice";
    end if;
    if Discriminant(U) eq 0 or GCD(U, f) ne 1 then
        return false, P!0, P!0, P!0, P!0, P!0, [], "bad U";
    end if;
    if GCD(B, U) ne 1 then
        return false, P!0, P!0, P!0, P!0, P!0, [], "B not invertible mod U";
    end if;

    Binv := InverseMod(B, U);
    V := (-A*Binv) mod U;
    if (V^2 - f) mod U ne 0 then
        return false, P!0, P!0, P!0, P!0, P!0, [], "V recovery failed";
    end if;

    independent, rels, reason := CheckClasses(f, h, U, V);
    if not independent then
        return false, A, B, f, U, V, rels, reason;
    end if;
    return true, A, B, f, U, V, rels, "independent";
end function;

procedure SymbolicReport()
    Q := Rationals();
    R<h1,h2,kap,s,t,b0,b1,b2> := PolynomialRing(Q, 8, "grevlex");
    Fr := FieldOfFractions(R);
    P<x> := PolynomialRing(Fr);

    h := 1 + h1*x + h2*x^2;
    U := x^2 + s*x + t;
    B := b0 + b1*x + b2*x^2;
    A, E, f, avec := ForcedAData(h, kap, U, B);
    res := [R!Numerator(Coefficient(E, i)) : i in [0..4]];

    print "# full Mumford norm normal form";
    print "# Normalizations: h(0)=1, A monic, Lambda=1.";
    print "# B is left as b0+b1*x+b2*x^2; the probes below require (b1,b2)!=(0,0).";
    print "h =", h;
    print "U =", U;
    print "B =", B;
    print "A =", A;
    print "";
    print "# residual equations E_i=coeff_x^i(A^2-B^2*f-U^5), i=0..4";
    for i in [0..4] do
        printf "E%o: degree=%o total_degree=%o terms=%o\n",
            i, Degree(res[i+1]), TotalDegree(res[i+1]), #Terms(res[i+1]);
    end for;

    print "";
    print "# immediate boundary checks";
    print "B=0 gives A^2=U^5, so U must be a double-root boundary for a polynomial A.";
    print "constant B is the old contact-type slice after rescaling; it is deliberately not probed here.";

    Rlin<h1l,h2l,kapl,sl,tl,b0l,b1l> := PolynomialRing(Q, 7, "grevlex");
    phi_lin := hom<R -> Rlin | h1l,h2l,kapl,sl,tl,b0l,b1l,0>;
    lin_res := [phi_lin(F) : F in res];
    print "";
    print "# linear-B branch b2=0 residual sizes";
    for i in [0..4] do
        printf "L%o: degree=%o total_degree=%o terms=%o\n",
            i, Degree(lin_res[i+1]), TotalDegree(lin_res[i+1]), #Terms(lin_res[i+1]);
    end for;
end procedure;

procedure RationalIntegerProbe()
    Q := Rationals();
    H := rational_height;
    tested := 0;
    bad_U_hits := 0;
    bad_curve_hits := 0;
    bad_B_mod_U_hits := 0;
    dependent_hits := 0;
    other_boundary_hits := 0;
    independent_hits := 0;

    printf "# rational integer probe on normalized full norm chart, box [-%o,%o]\n", H, H;
    for h1i in [-H..H] do
    for h2i in [-H..H] do
    for kapi in [-H..H] do
        if kapi eq 0 then
            continue;
        end if;
    for si in [-H..H] do
    for ti in [-H..H] do
    for b0i in [-H..H] do
    for b1i in [-H..H] do
    for b2i in [-H..H] do
        if b1i eq 0 and b2i eq 0 then
            continue;
        end if;
        tested +:= 1;
        ok, A, B, f, U, V, rels, reason := FullNormModel(
            Q, Q!h1i, Q!h2i, Q!kapi, Q!si, Q!ti, Q!b0i, Q!b1i, Q!b2i);
        if ok then
            independent_hits +:= 1;
            printf "RATIONAL_FULL_HIT h1=%o h2=%o K=%o s=%o t=%o b0=%o b1=%o b2=%o\n",
                h1i,h2i,kapi,si,ti,b0i,b1i,b2i;
            printf "  f=%o\n  U=%o\n  V=%o\n  A=%o\n  B=%o\n", f, U, V, A, B;
        elif reason ne "residual nonzero" then
            if reason eq "bad U" then
                bad_U_hits +:= 1;
            elif reason eq "bad contact curve" then
                bad_curve_hits +:= 1;
            elif reason eq "B not invertible mod U" then
                bad_B_mod_U_hits +:= 1;
            elif reason eq "checked" or reason eq "zero class" or reason eq "5*D2 != 0" then
                dependent_hits +:= 1;
                printf "rational_dependent_norm h1=%o h2=%o K=%o s=%o t=%o b0=%o b1=%o b2=%o reason=%o rels=%o\n",
                    h1i,h2i,kapi,si,ti,b0i,b1i,b2i,reason,rels;
            else
                other_boundary_hits +:= 1;
                printf "rational_norm_boundary h1=%o h2=%o K=%o s=%o t=%o b0=%o b1=%o b2=%o reason=%o rels=%o\n",
                    h1i,h2i,kapi,si,ti,b0i,b1i,b2i,reason,rels;
            end if;
        end if;
    end for;
    end for;
    end for;
    end for;
    end for;
    end for;
    end for;
    end for;

    printf "rational_tested=%o independent_hits=%o bad_U_boundary=%o bad_curve=%o bad_B_mod_U=%o dependent_norm=%o other_boundary=%o\n",
        tested, independent_hits, bad_U_hits, bad_curve_hits,
        bad_B_mod_U_hits, dependent_hits, other_boundary_hits;
end procedure;

procedure FiniteFieldNormProbe()
    print "# finite-field random probe of full norm equations with B nonconstant";
    for p in PrimesUpTo(prime_bound) do
        if p in {2,5} then
            continue;
        end if;

        F := GF(p);
        hit := false;
        boundary_hits := 0;
        for trial in [1..finite_trials] do
            h1v := Random(F);
            h2v := Random(F);
            kapv := Random(F);
            if kapv eq 0 then
                continue;
            end if;
            sv := Random(F);
            tv := Random(F);
            b0v := Random(F);
            b1v := Random(F);
            b2v := Random(F);
            if b1v eq 0 and b2v eq 0 then
                continue;
            end if;

            ok, A, B, f, U, V, rels, reason := FullNormModel(
                F, h1v, h2v, kapv, sv, tv, b0v, b1v, b2v);

            if ok then
                printf "p=%o FULL_NORM_HIT trial=%o h1=%o h2=%o K=%o s=%o t=%o b0=%o b1=%o b2=%o degB=%o #J=%o\n",
                    p, trial, Integers()!h1v, Integers()!h2v, Integers()!kapv,
                    Integers()!sv, Integers()!tv, Integers()!b0v, Integers()!b1v,
                    Integers()!b2v, Degree(B), #Jacobian(HyperellipticCurve(f));
                printf "  h=%o\n  f=%o\n  U=%o\n  V=%o\n  A=%o\n  B=%o\n",
                    1 + h1v*Parent(f).1 + h2v*Parent(f).1^2, f, U, V, A, B;
                hit := true;
                break;
            elif reason ne "residual nonzero" then
                boundary_hits +:= 1;
            end if;
        end for;

        if not hit then
            printf "p=%o no_full_norm_hit_in_%o_random_trials boundary_hits=%o\n",
                p, finite_trials, boundary_hits;
        end if;
    end for;
end procedure;

procedure FiniteFieldCantorProbe()
    print "# finite-field Cantor probe on one-contact family";
    print "# Randomly samples curves and degree-2 Mumford divisors; this is equivalent to [5][U,V]=0.";
    for p in PrimesUpTo(prime_bound) do
        if p in {2,5} then
            continue;
        end if;

        F := GF(p);
        P<x> := PolynomialRing(F);
        hit := false;
        tested := 0;

        for trial in [1..cantor_trials] do
            h1v := Random(F);
            h2v := Random(F);
            kapv := Random(F);
            if kapv eq 0 then
                continue;
            end if;
            h := 1 + h1v*x + h2v*x^2;
            f := h^2 - kapv*x^5;
            if Degree(f) ne 5 or Discriminant(f) eq 0 then
                continue;
            end if;

            C := HyperellipticCurve(f);
            J := Jacobian(C);
            if (#J mod 25) ne 0 then
                continue;
            end if;
            D0 := J![x, F!1];
            if 5*D0 ne J!0 then
                continue;
            end if;

            sv := Random(F);
            tv := Random(F);
            U := x^2 + sv*x + tv;
            if Discriminant(U) eq 0 or GCD(U, f) ne 1 then
                continue;
            end if;
            v1 := Random(F);
            v0 := Random(F);
            V := v1*x + v0;
            tested +:= 1;
            if (V^2 - f) mod U ne 0 then
                continue;
            end if;

            D2 := J![U, V];
            if D2 eq J!0 or 5*D2 ne J!0 then
                continue;
            end if;
            rels := RelationPairs(J, D0, D2);
            if #rels eq 0 then
                printf "p=%o CANTOR_HIT trial=%o h1=%o h2=%o K=%o s=%o t=%o v1=%o v0=%o #J=%o\n",
                    p, trial, Integers()!h1v, Integers()!h2v, Integers()!kapv,
                    Integers()!sv, Integers()!tv, Integers()!v1, Integers()!v0, #J;
                printf "  h=%o\n  f=%o\n  U=%o\n  V=%o\n", h, f, U, V;
                hit := true;
                break;
            end if;
        end for;

        if not hit then
            printf "p=%o no_cantor_hit_in_%o_curve_divisor_trials tested_UV=%o\n",
                p, cantor_trials, tested;
        end if;
    end for;
end procedure;

if do_symbolic then
    SymbolicReport();
    print "";
end if;

RationalIntegerProbe();
print "";
FiniteFieldNormProbe();
print "";
FiniteFieldCantorProbe();

quit;
