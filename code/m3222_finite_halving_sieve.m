//////////////////////////////////////////////////////////////////////
//  Finite-field check for divisibility by 2 of the distinguished
//  order-8 class in the odd M_1(8,2,2) family.
//
//  If a prime has no nonsingular (u,v) over F_p for which Q is
//  divisible by 2 in J(F_p), then rational examples over the open
//  family must reduce to the boundary modulo that prime.
//
//  Typical run from torsion_jac:
//      magma code/m3222_finite_halving_sieve.m
//////////////////////////////////////////////////////////////////////

prime_list := [3,5,7,11,13,17,19,23,29,31,37,41,43,47];
max_print := 8;

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

    checked := 0;
    nonsingular := 0;
    q_order8 := 0;
    divisible := 0;
    examples := [];

    for u in F do
        for v in F do
            checked +:= 1;
            if u eq 0 or v eq 0 or u eq v then
                continue;
            end if;
            if u eq 1 or v eq 1 or u+v+1 eq 0 or u+v+2 eq 0 then
                continue;
            end if;

            qtilde := -X^2
                + (u^2*v - u^2 + u*v^2 - u - v^2 - v - 2)*X
                - (u^2 + u*v + v^2 + u + v + 1);
            f := ((1-u)*X + 1)*((1-v)*X + 1)*((u+v+2)*X + 1)*qtilde;
            yq := u*v*(u+v+1);

            if Degree(f) ne 5 or Discriminant(f) eq 0 or yq eq 0 then
                continue;
            end if;
            if yq^2 ne Evaluate(f, F!-1) then
                continue;
            end if;

            nonsingular +:= 1;
            C := HyperellipticCurve(f);
            J := Jacobian(C);
            DQ := J![X + F!1, yq];
            if Order(DQ) eq 8 then
                q_order8 +:= 1;
            end if;
            if IsDivisibleBy2Finite(J, DQ) then
                divisible +:= 1;
                if #examples lt max_print then
                    Append(~examples, <u,v,Order(DQ),#J>);
                end if;
            end if;
        end for;
    end for;

    print "p", p,
          "checked", checked,
          "nonsingular", nonsingular,
          "Q_order8", q_order8,
          "Q_divisible_by_2", divisible;
    if #examples gt 0 then
        print "  examples", examples;
    else
        print "  NO OPEN FINITE-FIELD HALVES";
    end if;
end for;

quit;
