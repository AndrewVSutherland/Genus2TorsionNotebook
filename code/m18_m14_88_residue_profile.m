//////////////////////////////////////////////////////////////////////
//  Residue profile for the reduced [8,8] square-subcover.
//
//  This is a finite-prime diagnostic for the algebraic second-halving
//  equations, not a rational height search.  It works on the rational
//  square branch w=s^2, V=R^2*s^2.
//
//  The key simplification is that the reduced second equation is linear
//  in tau after the x^2 coefficient is used:
//
//      E2 = (P(z) + 4*D(z)*tau)/4.
//
//  Thus finite-field testing only has to enumerate z and the first-cover
//  tangent signs; tau is determined unless D=0.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

prime_list := [7,11,13,17,19,23,29,31,37,41,43];
if assigned primes then
    if Type(primes) eq MonStgElt then
        prime_list := [StringToInteger(s) : s in Split(primes, ",")];
    else
        prime_list := primes;
    end if;
end if;

max_pairs_print := 80;
if assigned max_pairs then
    if Type(max_pairs) eq MonStgElt then
        max_pairs_print := StringToInteger(max_pairs);
    else
        max_pairs_print := max_pairs;
    end if;
end if;

Z := Integers();

function GoodHyperellipticPolynomial(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function BoundaryLabels(F, R, w)
    labels := [];
    if R eq 0 then Append(~labels, "R"); end if;
    if w eq 0 then Append(~labels, "w"); end if;
    if w - 1 eq 0 then Append(~labels, "w-1"); end if;
    if w + 1 eq 0 then Append(~labels, "w+1"); end if;
    if R - 1 eq 0 then Append(~labels, "R-1"); end if;
    if R + 1 eq 0 then Append(~labels, "R+1"); end if;
    if R - w eq 0 then Append(~labels, "R-w"); end if;
    if R + w eq 0 then Append(~labels, "R+w"); end if;
    if R*w - 3*R + 3*w - 1 eq 0 then Append(~labels, "Lplus"); end if;
    if R*w + 3*R + 3*w + 1 eq 0 then Append(~labels, "Lminus"); end if;
    if 2*R^2 - R*w^2 + R - 2*w^2 eq 0 then Append(~labels, "Q"); end if;
    if R^4 - 2*R^3 + R^2*w^2 - R^2 + 2*R*w^2 - w^2 eq 0 then
        Append(~labels, "Quartic");
    end if;
    return labels;
end function;

function FamilyData(F, R, s)
    P<X> := PolynomialRing(F);
    w := s^2;
    t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
    A := X^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*X + R^4;
    B := (R + 2 + 4*t)*X^2 + (R^2 + 4*R + 1 + 8*t)*X
         + (2*R^2 + R + 4*t);
    h := A*B;
    return X*h, h, Coefficient(h, 1), Coefficient(h, 2),
           Coefficient(h, 3), Coefficient(h, 4);
end function;

function SquareRootsFinite(F, a)
    P<T> := PolynomialRing(F);
    return [rt[1] : rt in Roots(T^2 - a)];
end function;

function FirstCoverSolutionsFast(F, R, s)
    f, h, c1, c2, c3, c4 := FamilyData(F, R, s);
    P<U> := PolynomialRing(F);
    V := R^2*s^2;
    FU := 4*(c3 - 2*c4*U)*(c1 - 2*c4*U*V)
          - (c2 - c4*(U^2 + 2*V))^2;

    sols := [];
    for rt in Roots(FU) do
        U0 := rt[1];
        M2 := c3 - 2*c4*U0;
        N2 := c1 - 2*c4*U0*V;
        C := c2 - c4*(U0^2 + 2*V); // C = 2*M*N.

        for M0 in SquareRootsFinite(F, M2) do
            if M0 ne 0 then
                N0 := C/(2*M0);
                if N0^2 eq N2 then
                    Append(~sols, <U0, M0, N0, c4>);
                end if;
            else
                if C ne 0 then
                    continue;
                end if;
                for N0 in SquareRootsFinite(F, N2) do
                    Append(~sols, <U0, M0, N0, c4>);
                end for;
            end if;
        end for;
    end for;
    return sols;
end function;

function TauSolves(F, A, B, C)
    // A*tau^2 + B*tau + C has an F-rational root.
    if A eq 0 then
        if B eq 0 then
            return C eq 0, F!0;
        end if;
        return true, -C/B;
    end if;
    P<T> := PolynomialRing(F);
    rts := Roots(A*T^2 + B*T + C);
    if #rts eq 0 then
        return false, _;
    end if;
    return true, rts[1][1];
end function;

function ReducedSecondSolvableFast(F, R, s, sol)
    U := sol[1];
    M := sol[2];
    N := sol[3];
    c := sol[4];
    V := R^2*s^2;

    for eta in [F!1, -F!1] do
        mu := eta*R*s;
        A := 2*mu - U;
        for z in F do
            if z eq 0 then
                continue;
            end if;

            D := A - c*z^2;
            Pz := c^2*z^4 - 4*M*c*z^3
                  + (4*M^2 + 2*U*c)*z^2
                  + (4*M*U - 8*N)*z + U^2 - 4*V;
            B := 2*z*(M*mu - N) - mu*A - c*mu*z^2;
            C := V*c*z^2;

            if D ne 0 then
                tau := -Pz/(4*D);
                if A*tau^2 + B*tau + C eq 0 then
                    return true, <U,M,N,z,tau,eta>;
                end if;
            elif Pz eq 0 then
                ok, tau := TauSolves(F, A, B, C);
                if ok then
                    return true, <U,M,N,z,tau,eta>;
                end if;
            end if;
        end for;
    end for;
    return false, _;
end function;

print "Reduced [8,8] residue profile";
print "primes", prime_list;

for p in prime_list do
    if p eq 2 then
        print "p", p, "skipped_char2";
        continue;
    end if;
    F := GF(p);
    P<X> := PolynomialRing(F);

    open_good := 0;
    first_pairs := 0;
    first_solutions := 0;
    reduced_pairs := 0;
    reduced_first_solutions := 0;
    allowed_pairs := [];
    samples := [];

    for R in F do
        for s in F do
            w := s^2;
            if #BoundaryLabels(F, R, w) ne 0 then
                continue;
            end if;
            f, h, c1, c2, c3, c4 := FamilyData(F, R, s);
            if not GoodHyperellipticPolynomial(f) then
                continue;
            end if;
            open_good +:= 1;

            fs := FirstCoverSolutionsFast(F, R, s);
            if #fs eq 0 then
                continue;
            end if;
            first_pairs +:= 1;
            first_solutions +:= #fs;

            pair_ok := false;
            for sol in fs do
                ok, witness := ReducedSecondSolvableFast(F, R, s, sol);
                if ok then
                    reduced_first_solutions +:= 1;
                    pair_ok := true;
                    if #samples lt max_pairs_print then
                        Append(~samples, <Z!R, Z!s, witness>);
                    end if;
                end if;
            end for;
            if pair_ok then
                reduced_pairs +:= 1;
                if #allowed_pairs lt max_pairs_print then
                    Append(~allowed_pairs, <Z!R, Z!s>);
                end if;
            end if;
        end for;
    end for;

    print "p", p,
          "open_good", open_good,
          "first_pairs", first_pairs,
          "first_solutions", first_solutions,
          "reduced_pairs", reduced_pairs,
          "reduced_first_solutions", reduced_first_solutions;
    print "  allowed_pairs", allowed_pairs;
    print "  samples", samples;
end for;

quit;
