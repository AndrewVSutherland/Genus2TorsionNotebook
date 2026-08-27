// lane_x1fam.m — X1(N) family sweep (N in {7,8,9,10,12}), Tate normal form
// E(b,c): y^2 + (1-c)xy - by = x^3 - bx^2 with P=(0,0) of order N.
// Kubert parametrizations are runtime-verified (candidate variants tried at
// sample points; a family is skipped if none verifies).  Curves are bucketed
// by 2-torsion field; bucket collisions are glued (Genus2Elliptic2) and
// funneled.  Targets: [10,10], [12,12], [2,60], [2,40], [4,24], [9,9], [63].
// Usage: magma -b lane_x1fam.m > ../logs/x1fam.log
SetColumns(0);
SetMemoryLimit(4*10^9);
load "split_lab.m";  // run from product/code/

QT<t> := FunctionField(Rationals());

// candidate (b,c) parametrizations per N
CANDS := AssociativeArray();
CANDS[7]  := [* [* t^3-t^2, t^2-t *] *];
CANDS[8]  := [* [* (2*t-1)*(t-1), (2*t-1)*(t-1)/t *] *];
CANDS[9]  := [* [* t^2*(t-1)*(t^2-t+1), t^2*(t-1) *] *];
CANDS[10] := [*
  [* t^3*(t-1)*(2*t-1)/(t^2-3*t+1)^2, -t*(t-1)*(2*t-1)/(t^2-3*t+1) *],
  [* t^3*(t-1)*(2*t-1)/(t^2-3*t+1)^2,  t*(t-1)*(2*t-1)/(t^2-3*t+1) *],
  [* -t^3*(t-1)*(2*t-1)/(t^2-3*t+1)^2, -t*(t-1)*(2*t-1)/(t^2-3*t+1) *]
*];
CANDS[12] := [*
  [* -t*(2*t-1)*(3*t^2-3*t+1)*(2*t^2-2*t+1)/(t-1)^4, -t*(2*t-1)*(3*t^2-3*t+1)/(t-1)^3 *],
  [*  t*(2*t-1)*(3*t^2-3*t+1)*(2*t^2-2*t+1)/(t-1)^4,  t*(2*t-1)*(3*t^2-3*t+1)/(t-1)^3 *],
  [*  t*(2*t-1)*(3*t^2-3*t+1)*(2*t^2-2*t+1)/(t-1)^4, -t*(2*t-1)*(3*t^2-3*t+1)/(t-1)^3 *],
  [* -t*(2*t-1)*(3*t^2-3*t+1)*(2*t^2-2*t+1)/(t-1)^4,  t*(2*t-1)*(3*t^2-3*t+1)/(t-1)^3 *]
*];

function TryCurve(bv, cv)
    // Tate normal form; returns ok, E
    if bv eq 0 then return false, 0; end if;
    E := 0;
    try
        E := EllipticCurve([1-cv, -bv, -bv, 0, 0]);
    catch e
        return false, 0;
    end try;
    return true, E;
end function;

function VerifyFamily(N, pair)
    nok := 0;
    for tv in [Rationals()| 5, 7/2, -3/4 ] do
        ok := true; bv := 0; cv := 0;
        try
            bv := Evaluate(pair[1], tv); cv := Evaluate(pair[2], tv);
        catch e
            ok := false;
        end try;
        if not ok then continue; end if;
        ok, E := TryCurve(bv, cv);
        if not ok then continue; end if;
        P := E![0,0];
        if Order(P) eq N then nok +:= 1; end if;
    end for;
    return nok ge 2;
end function;

FAM := AssociativeArray();
for N in [7,8,9,10,12] do
    for pair in CANDS[N] do
        if VerifyFamily(N, pair) then
            FAM[N] := pair;
            printf "FAMILY N=%o verified\n", N;
            break;
        end if;
    end for;
    if not IsDefined(FAM, N) then printf "FAMILY N=%o FAILED verification -- skipped\n", N; end if;
end for;

function HeightRats(H)
    S := [];
    for a in [1..H], b in [1..H] do
        if GCD(a,b) eq 1 then Append(~S, a/b); Append(~S, -a/b); end if;
    end for;
    return S;
end function;

HN := AssociativeArray();
HN[7] := 20; HN[8] := 20; HN[9] := 20; HN[10] := 20; HN[12] := 20;

// collect curves: rec = <N, tval, E, torsinvs, jinv>
pool := [* *];
jseen := {};
for N in Sort(Setseq(Keys(FAM))) do
    pair := FAM[N];
    cnt := 0;
    for tv in HeightRats(HN[N]) do
        ok := true; bv := 0; cv := 0;
        try
            bv := Evaluate(pair[1], tv); cv := Evaluate(pair[2], tv);
        catch e
            ok := false;   // pole of the parametrization
        end try;
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

// bucket by 2-torsion field: rational-root cubics -> quadratic disc kernel;
// irreducible cubics -> disc kernel (stage 1), refined by field iso at pair time
buckets := AssociativeArray();
for i in [1..#pool] do
    rec := pool[i];
    fE := RQx!HyperellipticPolynomials(WeierstrassModel(rec[3]));
    rts := Roots(fE);
    if #rts eq 3 then
        key := "Q";
    elif #rts eq 1 then
        quad := fE div (RQx.1 - rts[1][1]);
        dd := Discriminant(quad);
        d0 := SquarefreeFactorization(Numerator(dd)*Denominator(dd));
        key := Sprintf("D%o", d0);
    else
        dd := Discriminant(fE);
        d0 := SquarefreeFactorization(Numerator(dd)*Denominator(dd));
        key := Sprintf("C%o", d0);
    end if;
    if not IsDefined(buckets, key) then buckets[key] := []; end if;
    Append(~(buckets[key]), i);
end for;

ALLOWED := { [7,7],[7,9],[9,9],[8,10],[8,12],[10,10],[10,12],[12,12],[8,8] };

seen := {};
nstat := AssociativeArray();
for st in ["skip","abort","known","hit","fail"] do nstat[st] := 0; end for;
npair := 0;
t0 := Cputime();
for key in Sort(Setseq(Keys(buckets))) do
    idxs := buckets[key];
    if #idxs lt 2 then continue; end if;
    printf "BUCKET %o size %o Ns=%o\n", key, #idxs, [ pool[i][1] : i in idxs ];
    for ii in [1..#idxs] do
        for jj in [ii+1..#idxs] do
            a := pool[idxs[ii]]; b := pool[idxs[jj]];
            NN := Sort([a[1], b[1]]);
            if not NN in ALLOWED then continue; end if;
            if a[5] eq b[5] then continue; end if;
            // for irreducible-cubic buckets confirm field isomorphism
            if key[1] eq "C" then
                fa := RQx!HyperellipticPolynomials(WeierstrassModel(a[3]));
                fb := RQx!HyperellipticPolynomials(WeierstrassModel(b[3]));
                if not IsIsomorphic(NumberField(fa), NumberField(fb)) then continue; end if;
                printf "CUBICMATCH N=%o t=%o and N=%o t=%o\n", a[1], a[2], b[1], b[2];
            end if;
            npair +:= 1;
            oddT := OddPartInvs([Integers()|d : d in a[4]] cat [Integers()|d : d in b[4]]);
            L := [];
            try L := Genus2Elliptic2(a[3], b[3]); catch e L := []; end try;
            for k in [1..#L] do
                gk := "";
                try gk := Sprintf("%o", G2Invariants(L[k])); catch e gk := "bad"; end try;
                if gk eq "bad" or gk in seen then continue; end if;
                Include(~seen, gk);
                st := Funnel(L[k], Sprintf("x1fam|%o(t=%o)|%o(t=%o)|%o", a[1], a[2], b[1], b[2], k) : OddInvs := oddT);
                nstat[st] +:= 1;
            end for;
            if npair mod 25 eq 0 then
                printf "PROGRESS x1fam %o pairs, curves %o, aborts %o known %o hits %o, %o s\n",
                    npair, #seen, nstat["abort"], nstat["known"], nstat["hit"], Cputime()-t0;
            end if;
        end for;
    end for;
end for;
printf "SEARCH_DONE x1fam pairs %o curves %o aborts %o known %o hits %o fails %o %o s\n",
    npair, #seen, nstat["abort"], nstat["known"], nstat["hit"], nstat["fail"], Cputime()-t0;
quit;
