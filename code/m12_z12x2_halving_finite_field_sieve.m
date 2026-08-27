//////////////////////////////////////////////////////////////////////
//  Finite-field sieve for halving the independent 2-torsion class in
//  the one-parameter M(12) family a=(1-r)/4.
//
//  If a good residue r mod p has the independent rational 2-torsion
//  class divisible by 2 over Q, its reduction must be divisible by 2
//  in J(F_p).  Thus a prime with no good divisible residues obstructs
//  nonboundary rational points at that prime.
//
//  Typical run from torsion_jac:
//      magma code/m12_z12x2_halving_finite_field_sieve.m
//////////////////////////////////////////////////////////////////////

prime_list := [3,5,7,11,13,17,19,23,29,31,37,41,43,47];

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
    P<X> := PolynomialRing(F);
    good := [];
    bad := [];
    divs := [];

    if p eq 2 then
        continue;
    end if;

    for r in F do
        if r in {F!0, F!1, F!2} then
            Append(~bad, <r, "degenerate-denominator">);
            continue;
        end if;

        a := (1-r)/4;
        T := a*X^2 - X + r;
        h := (X-r)*(T+1);
        W := h^2 + 4*a*X^2*T*(T+1);

        f5 := P!0;
        for i in [0..Degree(W)] do
            for j in [0..i] do
                f5 +:= Coefficient(W, i)*Binomial(i,j)*2^(i-j)*X^(6-j);
            end for;
        end for;

        if Degree(f5) ne 5 then
            Append(~bad, <r, "degree">);
            continue;
        end if;
        if Discriminant(f5) eq 0 then
            Append(~bad, <r, "disc">);
            continue;
        end if;

        beta_ind := (2-r)/(4*(r-1));
        beta_div := (1-r)/(4*r);
        C := HyperellipticCurve(f5);
        J := Jacobian(C);
        Tind := J![X-beta_ind, F!0];
        Tdiv := J![X-beta_div, F!0];

        known := IsDivisibleBy2Finite(J, Tdiv);
        independent := IsDivisibleBy2Finite(J, Tind);
        Append(~good, <r, known, independent>);
        if independent then
            Append(~divs, r);
            print "  independent divisible residue", p, r;
        end if;
    end for;

    print "p", p, "good", #good, "bad", #bad, "independent_divisible", #divs;
    print "  good data <r, known_6D_divisible, independent_divisible>:", good;
    print "  bad data:", bad;
    if #good gt 0 and #divs eq 0 then
        print "NO GOOD INDEPENDENT HALVES MOD", p;
        break;
    end if;
end for;
quit;
