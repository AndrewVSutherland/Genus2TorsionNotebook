//////////////////////////////////////////////////////////////////////
//  Finite-field diagnostic for halving the extra split-root 2-torsion
//  in the [2,4,8] one-split family.
//
//  Counts good affine residues where:
//    - the M_1(8,4) first halving holds;
//    - exactly one of A,B splits;
//    - one of the split rational branch-point classes is 2-divisible.
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

function GoodPolynomial(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function IsDoubleInFiniteJacobian(J, D)
    G, phi := AbelianGroup(J);
    d := D @@ phi;
    H := sub<G | [2*G.i : i in [1..Ngens(G)]]>;
    return d in H;
end function;

for p in prime_list do
    if p eq 2 then
        continue;
    end if;
    F := GF(p);
    P<x> := PolynomialRing(F);
    total := 0;
    good := 0;
    first := 0;
    one_split := 0;
    root_tests := 0;
    root_halvable := 0;
    base_halvable := 0;

    for R in F do
        for w in F do
            if R eq 0 or w eq 0 or w eq 1 or w eq -1 then
                continue;
            end if;
            total +:= 1;
            den := 4*(w^2-1);
            if den eq 0 then
                continue;
            end if;
            t := (2*R^2 + (1-w^2)*R - 2*w^2)/den;
            A := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
            B := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x
                 + (2*R^2 + R + 4*t);
            f := x*A*B;
            if not GoodPolynomial(f) then
                continue;
            end if;
            good +:= 1;

            C := HyperellipticCurve(f);
            J := Jacobian(C);
            Tx := J![x, F!0];
            ok_first := IsDoubleInFiniteJacobian(J, Tx);
            if not ok_first then
                continue;
            end if;
            first +:= 1;

            discA := Discriminant(A);
            discB := Discriminant(B);
            if discA eq 0 or discB eq 0 then
                continue;
            end if;
            splitA := IsSquare(discA);
            splitB := IsSquare(discB);
            if not ((splitA or splitB) and not (splitA and splitB)) then
                continue;
            end if;
            one_split +:= 1;

            roots := [];
            if splitA then
                roots cat:= [rt[1] : rt in Roots(A)];
            end if;
            if splitB then
                roots cat:= [rt[1] : rt in Roots(B)];
            end if;
            any_base := false;
            for alpha in roots do
                root_tests +:= 1;
                Dalpha := J![x-alpha, F!0];
                ok := IsDoubleInFiniteJacobian(J, Dalpha);
                if ok then
                    root_halvable +:= 1;
                    any_base := true;
                end if;
            end for;
            if any_base then
                base_halvable +:= 1;
            end if;
        end for;
    end for;

    print "p", p,
          "total", total,
          "good", good,
          "first", first,
          "one_split", one_split,
          "root_tests", root_tests,
          "root_halvable", root_halvable,
          "base_halvable", base_halvable;
end for;

quit;
