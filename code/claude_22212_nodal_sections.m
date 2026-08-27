// Hyperplane sections of S through many nodes: genus <= 1 sections would
// be [2,2,2,12] infinitude candidates.
SetColumns(0);
SetMemoryLimit(6*10^9);
Q := Rationals();
// node list
nodes := [];
for c in CartesianPower([-1,0,1], 5) do
    v := [c[i] : i in [1..5]];
    if &and[x eq 0 : x in v] then continue; end if;
    if v[1]^2+v[2]^2+v[3]^2 ne v[4]^2+v[5]^2 then continue; end if;
    if v[1]^4+v[2]^4+v[3]^4 ne v[4]^4+v[5]^4 then continue; end if;
    neg := [-x : x in v];
    if exists(t){ w : w in nodes | w eq neg } then continue; end if;
    Append(~nodes, v);
end for;
printf "nodes: %o\n", #nodes;

// enumerate hyperplanes spanned by 4-subsets; count nodes on each
seen := {};
best := [];
S4 := Subsets({1..#nodes}, 4);
for ss in S4 do
    idx := Setseq(ss);
    M := Matrix(Q, 4, 5, [ nodes[i] : i in idx ]);
    if Rank(M) ne 4 then continue; end if;
    ker := KernelMatrix(Transpose(M));
    h := Eltseq(ker[1]);
    d := LCM([Denominator(x) : x in h]);
    h := [Q!(x*d) : x in h];
    g := GCD([Integers()!x : x in h]);
    h := [x/g : x in h];
    if h[1] lt 0 or (h[1] eq 0 and (h[2] lt 0 or (h[2] eq 0 and (h[3] lt 0 or (h[3] eq 0 and (h[4] lt 0 or (h[4] eq 0 and h[5] lt 0))))))) then
        h := [-x : x in h];
    end if;
    if h in seen then continue; end if;
    Include(~seen, h);
    // count nodes on hyperplane (mod +-, both signs count as incidences of the projective node)
    cnt := #[ v : v in nodes | &+[h[k]*v[k] : k in [1..5]] eq 0 ];
    if cnt ge 7 then Append(~best, <cnt, h>); end if;
end for;
printf "distinct spanned hyperplanes: %o\n", #seen;
Sort(~best, func<x,y | y[1] - x[1]>);
printf "hyperplanes with >=7 nodes: %o\n", #best;
for i in [1..Min(20, #best)] do
    printf "  H %o nodes: %o\n", best[i][1], best[i][2];
end for;
printf "NODAL_HYPS_DONE\n";
quit;
