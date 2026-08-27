//////////////////////////////////////////////////////////////////////
//  Finite-field verifier for the reduced [8,8] halving equations.
//
//  This solves the reduced equations from m18_m14_88_reduced_conditions.m
//  over small finite fields and compares the result with the intrinsic
//  condition T_x in 4*J(F_p) on the good open chart.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

prime_list := [7,11,13,17];
if assigned primes then
    if Type(primes) eq MonStgElt then
        prime_list := [StringToInteger(s) : s in Split(primes, ",")];
    else
        prime_list := primes;
    end if;
end if;

max_samples := 12;
if assigned samples then
    if Type(samples) eq MonStgElt then
        max_samples := StringToInteger(samples);
    else
        max_samples := samples;
    end if;
end if;

Z := Integers();

function GoodHyperellipticPolynomial(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function DivisibleByNFinite(J, D, n)
    G, phi := AbelianGroup(J);
    a := D @@ phi;
    coords := Eltseq(a);
    invs := Invariants(G);
    for i in [1..#coords] do
        g := GCD(n, invs[i]);
        if coords[i] mod g ne 0 then
            return false;
        end if;
    end for;
    return true;
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

function FamilyPolynomialFinite(F, R, w)
    P<X> := PolynomialRing(F);
    t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
    A := X^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*X + R^4;
    B := (R + 2 + 4*t)*X^2 + (R^2 + 4*R + 1 + 8*t)*X
         + (2*R^2 + R + 4*t);
    return X*A*B, A*B, Coefficient(A*B, 4);
end function;

function FirstSolutions(F, R, s)
    P<X> := PolynomialRing(F);
    w := s^2;
    f, h, c4 := FamilyPolynomialFinite(F, R, w);
    V := R^2*s^2;
    sols := [];
    for U in F do
        a0 := X^2 + U*X + V;
        for M in F do
            for N in F do
                first := h - X*(M*X+N)^2 - c4*a0^2;
                if &and [Coefficient(first, i) eq 0 : i in [1..4]] then
                    Append(~sols, <U,M,N>);
                end if;
            end for;
        end for;
    end for;
    return sols;
end function;

function ReducedSecondSolvable(F, R, s, first_sols)
    P<X> := PolynomialRing(F);
    w := s^2;
    f, h, c4 := FamilyPolynomialFinite(F, R, w);
    V := R^2*s^2;
    eta_vals := [F!1, -F!1];
    for sol in first_sols do
        U := sol[1];
        M := sol[2];
        N := sol[3];
        a0 := X^2 + U*X + V;
        ell0 := X*(M*X+N);
        for z in F do
            if z eq 0 then
                continue;
            end if;
            for tau in F do
                for eta in eta_vals do
                    a := U/2 + M*z + tau - (c4/2)*z^2;
                    b := eta*R*s*tau;
                    q := X^2 + a*X + b;
                    second := z^2*X*c4*a0 - 2*z*ell0*(X+tau)
                              - a0*(X+tau)^2 + q^2;
                    if Coefficient(second, 1) eq 0
                       and Coefficient(second, 2) eq 0 then
                        return true, <U,M,N,z,tau,eta,a,b>;
                    end if;
                end for;
            end for;
        end for;
    end for;
    return false, _;
end function;

print "REDUCED [8,8] finite-field verifier";
print "primes", prime_list;

for p in prime_list do
    if p eq 2 then
        print "p", p, "skipped_char2";
        continue;
    end if;
    F := GF(p);
    P<X> := PolynomialRing(F);
    good_square_cover_bases := 0;
    intrinsic_88 := 0;
    reduced_88 := 0;
    reduced_without_intrinsic := 0;
    intrinsic_without_reduced := 0;
    first_total := 0;
    samples := [];
    misses := [];

    for delta in [F!1, -F!1] do
        for R in F do
            for s in F do
                w := delta*s^2;
                if #BoundaryLabels(F, R, w) ne 0 then
                    continue;
                end if;
                f, h, c4 := FamilyPolynomialFinite(F, R, w);
                if not GoodHyperellipticPolynomial(f) then
                    continue;
                end if;
                good_square_cover_bases +:= 1;

                C := HyperellipticCurve(f);
                J := Jacobian(C);
                Tx := J![X, F!0];
                intr_ok := DivisibleByNFinite(J, Tx, 4);
                if intr_ok then
                    intrinsic_88 +:= 1;
                end if;

                fs := FirstSolutions(F, R, s);
                if #fs gt 0 then
                    first_total +:= 1;
                end if;
                reduced, sol := ReducedSecondSolvable(F, R, s, fs);
                if reduced then
                    reduced_88 +:= 1;
                    if #samples lt max_samples then
                        Append(~samples, <Z!delta, Z!R, Z!s, Z!w, sol>);
                    end if;
                end if;

                if reduced and not intr_ok then
                    reduced_without_intrinsic +:= 1;
                    if #misses lt max_samples then
                        Append(~misses, <"reduced_not_intrinsic", Z!delta, Z!R, Z!s, Z!w, sol>);
                    end if;
                elif intr_ok and not reduced then
                    intrinsic_without_reduced +:= 1;
                    if #misses lt max_samples then
                        Append(~misses, <"intrinsic_not_reduced", Z!delta, Z!R, Z!s, Z!w>);
                    end if;
                end if;
            end for;
        end for;
    end for;

    print "p", p,
          "good_square_cover_bases", good_square_cover_bases,
          "first_cover_bases", first_total,
          "intrinsic_Tx_fourdiv", intrinsic_88,
          "reduced_88", reduced_88,
          "reduced_without_intrinsic", reduced_without_intrinsic,
          "intrinsic_without_reduced", intrinsic_without_reduced;
    print "  samples", samples;
    print "  mismatches", misses;
end for;

quit;
