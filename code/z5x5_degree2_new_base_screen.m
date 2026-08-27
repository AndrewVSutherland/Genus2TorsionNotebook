//////////////////////////////////////////////////////////////////////
// Actual finite-group screen for the degree-2 [5,5] contact search.
//
// Typical run:
//   magma -b p1:=7 p2:=11 code/z5x5_degree2_new_base_screen.m
//
// A CELL line is printed exactly when the smooth genus-2 curve
//   y^2 = (1+a*x+b*x^2)^2-k*x^5
// has J(F_p)[5] of rank at least two.  This uses the invariant factors of
// the finite Jacobian, rather than the weaker test 25 | #J(F_p).
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned p1 then
    p1 := 7;
elif Type(p1) eq MonStgElt then
    p1 := StringToInteger(p1);
end if;
if not assigned p2 then
    p2 := 11;
elif Type(p2) eq MonStgElt then
    p2 := StringToInteger(p2);
end if;

Z := Integers();

function FiveRank(invs)
    return #[n : n in invs | (Z!n) mod 5 eq 0];
end function;

procedure ScreenPrime(p)
    K := GF(p);
    P<x> := PolynomialRing(K);
    tested := 0;
    smooth := 0;
    rank_counts := AssociativeArray();
    allowed := 0;

    for a in K do
        for b in K do
            for k in K do
                if k eq 0 then
                    continue;
                end if;
                tested +:= 1;
                f := (1+a*x+b*x^2)^2-k*x^5;
                if Degree(f) ne 5 or Discriminant(f) eq 0 then
                    continue;
                end if;
                smooth +:= 1;
                C := HyperellipticCurve(f);
                A, phi := AbelianGroup(Jacobian(C));
                invs := Invariants(A);
                rank := FiveRank(invs);
                if IsDefined(rank_counts, rank) then
                    rank_counts[rank] +:= 1;
                else
                    rank_counts[rank] := 1;
                end if;
                if rank ge 2 then
                    allowed +:= 1;
                    printf "CELL %o %o %o %o %o %o\n",
                           p, Z!a, Z!b, Z!k, invs, #A;
                end if;
            end for;
        end for;
    end for;

    printf "SUMMARY p=%o tested=%o smooth=%o allowed_rank_ge_2=%o rank_counts=%o\n",
           p, tested, smooth, allowed,
           Sort([<r, rank_counts[r]> : r in Keys(rank_counts)]);
end procedure;

print "# degree-2 [5,5] actual finite-group residue screen";
print "# columns: CELL p a b k invariants order";
for p in [p1, p2] do
    require IsPrime(p) and p notin {2,5}: "p must be prime and different from 2,5";
    ScreenPrime(p);
end for;
