// Decompose sections of S by hyperplanes through >=8 nodes (excluding
// degeneracy-form hyperplanes); hunt nondegenerate components of genus <= 1.
SetColumns(0);
SetMemoryLimit(8*10^9);
Q := Rationals();
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
seen := {};
cands := [];
for ss in Subsets({1..#nodes}, 4) do
    idx := Setseq(ss);
    M := Matrix(Q, 4, 5, [ nodes[i] : i in idx ]);
    if Rank(M) ne 4 then continue; end if;
    h := Eltseq(KernelMatrix(Transpose(M))[1]);
    d := LCM([Denominator(x) : x in h]);
    hh := [Integers()!(x*d) : x in h];
    g := GCD(hh); hh := [x div g : x in hh];
    for k in [1..5] do
        if hh[k] ne 0 then if hh[k] lt 0 then hh := [-x : x in hh]; end if; break; end if;
    end for;
    if hh in seen then continue; end if;
    Include(~seen, hh);
    nz := #[x : x in hh | x ne 0];
    if nz le 2 then continue; end if;    // degeneracy-form hyperplanes x_i, x_i+-x_j
    cnt := #[ v : v in nodes | &+[hh[k]*v[k] : k in [1..5]] eq 0 ];
    if cnt ge 8 then Append(~cands, <cnt, hh>); end if;
end for;
Sort(~cands, func<x,y | y[1] - x[1]>);
printf "candidate hyperplanes (>=8 nodes, honest): %o\n", #cands;

P3<y1,y2,y3,y4> := ProjectiveSpace(Q, 3);
nan := 0;
for cd in cands do
    if nan ge 248 then printf "  (cap reached)\n"; break; end if;
    nan +:= 1;
    hh := cd[2];
    // solve h.x = 0 for the last nonzero coordinate
    piv := Max([ k : k in [1..5] | hh[k] ne 0 ]);
    X := []; j := 0;
    ys := [y1,y2,y3,y4];
    for k in [1..5] do
        if k ne piv then j +:= 1; X[k] := ys[j]; end if;
    end for;
    X[piv] := -&+[ Q!hh[k]/hh[piv] * X[k] : k in [1..5] | k ne piv ];
    Q2 := X[1]^2+X[2]^2+X[3]^2-X[4]^2-X[5]^2;
    F4 := X[1]^4+X[2]^4+X[3]^4-X[4]^4-X[5]^4;
    Se := Scheme(P3, [Q2, F4]);
    comps := IrreducibleComponents(Se);
    printf "H=%o (%o nodes): %o components:", hh, cd[1], #comps;
    for co in comps do
        if Dimension(co) ne 1 then continue; end if;
        cor := ReducedSubscheme(co);
        dg := Degree(cor);
        // degeneracy check: any of the 25 forms identically zero on component?
        Icor := Ideal(cor);
        forms := [ X[i] : i in [1..5] ];
        for i1 in [1..4] do for j1 in [i1+1..5] do
            Append(~forms, X[i1]-X[j1]); Append(~forms, X[i1]+X[j1]);
        end for; end for;
        bad := exists(t){ L : L in forms | L in Icor };
        if bad then printf " [deg %o DEGEN]", dg; continue; end if;
        Cu := Curve(cor);
        gg := Genus(Cu);
        printf " [deg %o genus %o NONDEG]", dg, gg; if gg eq 1 then try pts := PointSearch(Curve(cor), 200); if #pts gt 0 then E := EllipticCurve(Curve(cor), pts[1]); Em := MinimalModel(E); okr, rr := RankBounds(Em); printf " aInv=%o tors=%o rk=[%o,%o]", aInvariants(Em), Invariants(TorsionSubgroup(Em)), okr, rr; end if; catch ee printf " ERR"; end try; end if;
        if gg le 1 then printf " <<< CANDIDATE H=%o", hh; end if;
    end for;
    printf "\n";
end for;
printf "SECTIONS_DONE\n";
quit;
