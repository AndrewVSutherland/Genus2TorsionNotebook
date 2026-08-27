// lane_gainhunt.m — hunt for torsion GAINS above the generic glued image, by
// imposing non-generic factorization on the glued sextic (extra rational
// Weierstrass structure = extra rational 2-torsion).  Families X1(N),
// N in {6,8,10,12}, larger heights than lane_x1fam.
// Targets: [2,2,24] ((6,8) with >=2 rational roots), [10,10] ((10,10) with
// split quartic), [2,60] ((10,12)), [2,40] ((8,10)), [4,24] ((8,12)),
// [2,2,6,6] ((6,6) fully-split), [12,12] ((12,12), funnel-all).
// Usage: magma -b lane_gainhunt.m > ../logs/gainhunt.log
SetColumns(0);
SetMemoryLimit(4*10^9);
load "split_lab.m";  // run from product/code/

QT<t> := FunctionField(Rationals());

CANDS := AssociativeArray();
CANDS[6]  := [* [* t^2+t, t *], [* t-t^2, t *], [* t^2-t, -t *] *];
CANDS[8]  := [* [* (2*t-1)*(t-1), (2*t-1)*(t-1)/t *] *];
CANDS[10] := [*
  [* t^3*(t-1)*(2*t-1)/(t^2-3*t+1)^2, -t*(t-1)*(2*t-1)/(t^2-3*t+1) *],
  [* t^3*(t-1)*(2*t-1)/(t^2-3*t+1)^2,  t*(t-1)*(2*t-1)/(t^2-3*t+1) *]
*];
CANDS[12] := [*
  [* -t*(2*t-1)*(3*t^2-3*t+1)*(2*t^2-2*t+1)/(t-1)^4, -t*(2*t-1)*(3*t^2-3*t+1)/(t-1)^3 *],
  [*  t*(2*t-1)*(3*t^2-3*t+1)*(2*t^2-2*t+1)/(t-1)^4,  t*(2*t-1)*(3*t^2-3*t+1)/(t-1)^3 *]
*];

function TryCurve(bv, cv)
    if bv eq 0 then return false, 0; end if;
    E := 0;
    try E := EllipticCurve([1-cv, -bv, -bv, 0, 0]); catch e return false, 0; end try;
    return true, E;
end function;

function VerifyFamily(N, pair)
    nok := 0;
    for tv in [Rationals()| 5, 7/2, -3/4 ] do
        ok := true; bv := 0; cv := 0;
        try bv := Evaluate(pair[1], tv); cv := Evaluate(pair[2], tv); catch e ok := false; end try;
        if not ok then continue; end if;
        ok, E := TryCurve(bv, cv);
        if not ok then continue; end if;
        P := E![0,0];
        if Order(P) eq N then nok +:= 1; end if;
    end for;
    return nok ge 2;
end function;

FAM := AssociativeArray();
for N in [6,8,10,12] do
    for pair in CANDS[N] do
        if VerifyFamily(N, pair) then FAM[N] := pair; printf "FAMILY N=%o verified\n", N; break; end if;
    end for;
    if not IsDefined(FAM, N) then printf "FAMILY N=%o FAILED\n", N; end if;
end for;

function HeightRats(H)
    S := [];
    for a in [1..H], b in [1..H] do
        if GCD(a,b) eq 1 then Append(~S, a/b); Append(~S, -a/b); end if;
    end for;
    return S;
end function;

HN := AssociativeArray();
HN[6] := 24; HN[8] := 24; HN[10] := 32; HN[12] := 32;

pool := [* *];
jseen := {};
for N in Sort(Setseq(Keys(FAM))) do
    pair := FAM[N];
    cnt := 0;
    for tv in HeightRats(HN[N]) do
        ok := true; bv := 0; cv := 0;
        try bv := Evaluate(pair[1], tv); cv := Evaluate(pair[2], tv); catch e ok := false; end try;
        if not ok then continue; end if;
        ok, E := TryCurve(bv, cv);
        if not ok then continue; end if;
        P := 0;
        try P := E![0,0]; catch e continue; end try;
        if Order(P) ne N then continue; end if;
        j := jInvariant(E);
        if j in jseen then continue; end if;
        Include(~jseen, j);
        Tinv := Invariants(TorsionSubgroup(E));
        Append(~pool, <N, tv, E, Tinv, j>);
        cnt +:= 1;
    end for;
    printf "POOL N=%o: %o curves\n", N, cnt;
end for;

buckets := AssociativeArray();
for i in [1..#pool] do
    rec := pool[i];
    fE := RQx!HyperellipticPolynomials(WeierstrassModel(rec[3]));
    rts := Roots(fE);
    if #rts eq 3 then key := "Q";
    elif #rts eq 1 then
        quad := fE div (RQx.1 - rts[1][1]);
        dd := Discriminant(quad);
        d0 := SquarefreeFactorization(Numerator(dd)*Denominator(dd));
        key := Sprintf("D%o", d0);
    else
        continue;   // irreducible cubics irrelevant for these N
    end if;
    if not IsDefined(buckets, key) then buckets[key] := []; end if;
    Append(~(buckets[key]), i);
end for;

// pair policy: <allowed, funnel-all?>; otherwise require non-generic pattern
ALLOWED := { [6,6],[6,8],[6,10],[6,12],[8,10],[8,12],[10,10],[10,12],[12,12] };
FUNNELALL := { [10,10],[10,12],[12,12] };

function FactorPatternOK(g, isQ)
    // degree-only use of Factorization (safe); true = non-generic, worth exact
    fac := Factorization(g);
    degs := Sort([Degree(t[1]) : t in fac]);
    nlin := #[d : d in degs | d eq 1];
    if isQ then
        return nlin ge 2;         // generic class-(a) glue: three quadratics
    else
        return #degs ge 3;        // generic class-(b) glue: quadratic x quartic
    end if;
end function;

seen := {};
nstat := AssociativeArray();
for st in ["skip","abort","known","hit","fail"] do nstat[st] := 0; end for;
npair := 0; nfun := 0;
t0 := Cputime();
for key in Sort(Setseq(Keys(buckets))) do
    idxs := buckets[key];
    if #idxs lt 2 then continue; end if;
    isQ := key eq "Q";
    printf "BUCKET %o size %o Ns=%o\n", key, #idxs, {* pool[i][1] : i in idxs *};
    for ii in [1..#idxs] do
        for jj in [ii+1..#idxs] do
            a := pool[idxs[ii]]; b := pool[idxs[jj]];
            NN := Sort([a[1], b[1]]);
            if not NN in ALLOWED then continue; end if;
            if a[5] eq b[5] then continue; end if;
            npair +:= 1;
            oddT := OddPartInvs([Integers()|d : d in a[4]] cat [Integers()|d : d in b[4]]);
            L := [];
            try L := Genus2Elliptic2(a[3], b[3]); catch e L := []; end try;
            for k in [1..#L] do
                gk := "";
                try gk := Sprintf("%o", G2Invariants(L[k])); catch e gk := "bad"; end try;
                if gk eq "bad" or gk in seen then continue; end if;
                Include(~seen, gk);
                okg, g := IntegralSextic(L[k]);
                if not okg then continue; end if;
                if not (NN in FUNNELALL) and not FactorPatternOK(g, isQ) then continue; end if;
                nfun +:= 1;
                st := Funnel(L[k], Sprintf("gain|%o(t=%o)|%o(t=%o)|%o", a[1], a[2], b[1], b[2], k) : OddInvs := oddT);
                nstat[st] +:= 1;
            end for;
            if npair mod 200 eq 0 then
                printf "PROGRESS gain %o pairs, funneled %o, aborts %o known %o hits %o, %o s\n",
                    npair, nfun, nstat["abort"], nstat["known"], nstat["hit"], Cputime()-t0;
            end if;
        end for;
    end for;
end for;
printf "SEARCH_DONE gainhunt pairs %o funneled %o aborts %o known %o hits %o fails %o %o s\n",
    npair, nfun, nstat["abort"], nstat["known"], nstat["hit"], nstat["fail"], Cputime()-t0;
quit;
