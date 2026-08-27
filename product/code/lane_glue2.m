// lane_glue2.m — sweep of (2,2)-gluings Genus2Elliptic2(E1,E2) over pairs of
// LMFDB elliptic curves with matching 2-torsion fields, aimed at split-table
// gaps: [8,8], [2,2,24], [10,10], [2,6,12], [2,40], [4,24], [2,60], [12,12], [9,9], ...
// Usage: magma -b Part:=1 NParts:=2 lane_glue2.m > ../logs/glue2_p1.log
SetColumns(0);
if not assigned Part then Part := 1; elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;
if not assigned NParts then NParts := 2; elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;
SetMemoryLimit(4*10^9);
load "split_lab.m";  // run from product/code/
load "../data/seeds_magma.m";

// ---------- seed selection ----------
function Cap(t)
    if t eq [Integers()|2,2] then return 30;
    elif t eq [2,4] then return 60;
    elif t eq [2,6] then return 96;
    elif t eq [2,8] then return 100;
    elif t eq [5] then return 300;
    elif t eq [6] then return 300;
    elif t eq [7] then return 300;
    elif t eq [8] then return 300;
    elif t eq [9] then return 300;
    elif t eq [10] then return 300;
    elif t eq [12] then return 300;
    else return 0;
    end if;
end function;

counts := AssociativeArray();
kept := [* *];
for s in SEEDS do
    t := s[4];
    c := Cap(t);
    if c eq 0 then continue; end if;
    key := Sprintf("%o", t);
    if not IsDefined(counts, key) then counts[key] := 0; end if;
    if counts[key] ge c then continue; end if;
    counts[key] +:= 1;
    Append(~kept, s);
end for;
printf "kept %o seeds\n", #kept;

// ---------- bucket by 2-torsion splitting field ----------
R<x> := PolynomialRing(Rationals());
buckets := AssociativeArray();
t0 := Cputime();
for s in kept do
    E := EllipticCurve([Rationals()!a : a in s[5]]);
    fE := HyperellipticPolynomials(WeierstrassModel(E));
    rts := Roots(fE);
    if #rts eq 3 then
        key := "Q";
    elif #rts eq 1 then
        quad := fE div (x - rts[1][1]);
        dd := Discriminant(quad);
        d0 := SquarefreeFactorization(Numerator(dd)*Denominator(dd));
        key := Sprintf("D%o", d0);
    else
        K := SplittingField(fE);
        key := Sprintf("S%o", Coefficients(DefiningPolynomial(polredabs(K))));
    end if;
    ord := IsEmpty(s[4]) select 1 else &*s[4];
    rec := <s[1], s[4], ord, jInvariant(E), E>;
    if not IsDefined(buckets, key) then buckets[key] := [* *]; end if;
    Append(~(buckets[key]), rec);
end for;
printf "bucketing done in %o s; %o buckets\n", Cputime()-t0, #Keys(buckets);

// ---------- pair generation ----------
function PairTypeKey(t1, t2)
    s1 := Sprintf("%o", t1); s2 := Sprintf("%o", t2);
    return s1 le s2 select s1 cat "|" cat s2 else s2 cat "|" cat s1;
end function;

TYPECAP := AssociativeArray();
TYPECAP["[ 2, 6 ]|[ 2, 6 ]"] := 300;
TYPECAP["[ 2, 4 ]|[ 2, 6 ]"] := 300;
TYPECAP["[ 2, 4 ]|[ 2, 4 ]"] := 0;   // dropped

pairs := [* *];
typecount := AssociativeArray();
for key in Sort(Setseq(Keys(buckets))) do
    mem := buckets[key];
    n := #mem;
    if n lt 2 then continue; end if;
    printf "BUCKET %o size %o : %o\n", key, n, [* m[1] : m in mem *];
    thresh := key eq "Q" select 96 else 40;
    // deterministic order: as loaded (conductor order within torsion group)
    for i in [1..n] do
        for j in [i+1..n] do
            a := mem[i]; b := mem[j];
            if a[4] eq b[4] then continue; end if;   // same j-invariant
            sc := a[3]*b[3];
            if sc lt thresh then continue; end if;
            tk := PairTypeKey(a[2], b[2]);
            cap := IsDefined(TYPECAP, tk) select TYPECAP[tk] else 10^9;
            if not IsDefined(typecount, tk) then typecount[tk] := 0; end if;
            if typecount[tk] ge cap then continue; end if;
            typecount[tk] +:= 1;
            Append(~pairs, <a, b, sc>);
        end for;
    end for;
end for;
printf "total pairs: %o\n", #pairs;
for tk in Sort(Setseq(Keys(typecount))) do printf "PAIRTYPE %o : %o\n", tk, typecount[tk]; end for;

// ---------- main loop ----------
seen := {};
nstat := AssociativeArray();
for st in ["skip","abort","known","hit","fail"] do nstat[st] := 0; end for;
done := 0;
t0 := Cputime();
for idx in [1..#pairs] do
    if (idx mod NParts) ne (Part mod NParts) then continue; end if;
    pr := pairs[idx];
    a := pr[1]; b := pr[2];
    oddT := OddPartInvs(a[2] cat b[2]);   // odd torsion of glued J is exactly this
    L := [];
    try
        L := Genus2Elliptic2(a[5], b[5]);
    catch e
        L := [];
    end try;
    for k in [1..#L] do
        C := L[k];
        gk := "";
        try gk := Sprintf("%o", G2Invariants(C)); catch e gk := Sprintf("new%o.%o", idx, k); end try;
        if gk in seen then continue; end if;
        Include(~seen, gk);
        st := Funnel(C, Sprintf("glue2|%o|%o|%o", a[1], b[1], k) : OddInvs := oddT);
        nstat[st] +:= 1;
    end for;
    done +:= 1;
    if done mod 20 eq 0 then
        printf "PROGRESS part %o: %o pairs, %o curves, aborts %o known %o hits %o fails %o, %o s\n",
            Part, done, #seen, nstat["abort"], nstat["known"], nstat["hit"], nstat["fail"], Cputime()-t0;
    end if;
end for;
printf "SEARCH_DONE part %o: %o pairs, curves %o, aborts %o known %o hits %o fails %o skips %o, %o s\n",
    Part, done, #seen, nstat["abort"], nstat["known"], nstat["hit"], nstat["fail"], nstat["skip"], Cputime()-t0;
quit;
