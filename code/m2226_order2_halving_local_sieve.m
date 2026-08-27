//////////////////////////////////////////////////////////////////////
//  Local sieve for the first route on M(2,2,2,6): can any rational
//  2-torsion class be divided by 2?
//
//  For each unordered pair of branch points, this checks the necessary
//  squareclass conditions over F_p away from the bad locus.  If a class
//  has no nondegenerate F_p point for some p, then it has no rational
//  point with good reduction at p.  Boundary reductions need a separate
//  p-adic/blowup check.
//
//  Typical run from torsion_jac:
//      magma code/m2226_order2_halving_local_sieve.m
//////////////////////////////////////////////////////////////////////

prime_list := [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101];

function Branches(F, ss, mm, nn)
    s := F!ss;
    m := F!mm;
    n := F!nn;
    A := [ F!1, F!1, F!1, F!2, F!2 ];
    B := [
        2*s^2 - s*n,
        2*s^2 + s*m - 2*s*n - m*n,
        2*s^2 + s*m - s*n - m*n,
        -m*n,
        4*s^2 - 4*s*n - m*n
    ];
    return [ <F!1, F!0> ] cat [ < -A[i], B[i] > : i in [1..5] ];
end function;

function DetP(P, Qp)
    return P[1]*Qp[2] - Qp[1]*P[2];
end function;

function BranchesDistinct(br)
    for i in [1..5] do
        for j in [i+1..6] do
            if DetP(br[i], br[j]) eq 0 then
                return false;
            end if;
        end for;
    end for;
    return true;
end function;

function OrderedPairHalvesOverFp(br, i_inf, i_zero)
    remaining := [ k : k in [1..6] | k ne i_inf and k ne i_zero ];
    ref := remaining[1];
    for idx in [2..#remaining] do
        k := remaining[idx];
        num := DetP(br[k], br[i_zero]) * DetP(br[ref], br[i_inf]);
        den := DetP(br[k], br[i_inf]) * DetP(br[ref], br[i_zero]);
        if den eq 0 or num eq 0 then
            return false;
        end if;
        if not IsSquare(num/den) then
            return false;
        end if;
    end for;
    return true;
end function;

function ClassHasPointOverFp(p, unordered)
    F := GF(p);
    i := unordered[1];
    j := unordered[2];

    points_checked := 0;
    for ss in F do
        for mm in F do
            for nn in F do
                if ss eq 0 and mm eq 0 and nn eq 0 then
                    continue;
                end if;
                // Normalize projectively by first nonzero coordinate.
                if ss ne 0 and ss ne 1 then
                    continue;
                end if;
                if ss eq 0 and mm ne 0 and mm ne 1 then
                    continue;
                end if;
                if ss eq 0 and mm eq 0 and nn ne 1 then
                    continue;
                end if;

                br := Branches(F, ss, mm, nn);
                if not BranchesDistinct(br) then
                    continue;
                end if;
                points_checked +:= 1;

                if OrderedPairHalvesOverFp(br, i, j) or OrderedPairHalvesOverFp(br, j, i) then
                    return true, [ ss, mm, nn ], points_checked;
                end if;
            end for;
        end for;
    end for;

    return false, [], points_checked;
end function;

classes := [];
for i in [1..5] do
    for j in [i+1..6] do
        Append(~classes, [i,j]);
    end for;
end for;

for cls in classes do
    print "class", cls;
    obstructed := false;
    for p in prime_list do
        ok, pt, checked := ClassHasPointOverFp(p, cls);
        if checked eq 0 then
            print "  p", p, "skipped: no good reductions";
            continue;
        elif ok then
            print "  p", p, "has point", pt, "checked", checked;
        else
            print "  p", p, "NO POINT", "checked", checked;
            obstructed := true;
            break;
        end if;
    end for;
    if not obstructed then
        print "  no obstruction found in prime list";
    end if;
end for;
quit;
