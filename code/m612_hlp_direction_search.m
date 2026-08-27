//////////////////////////////////////////////////////////////////////
// Search small transverse coefficient directions through the split HLP
// [6,12] seed
//
//   F_0 = 187392*x^6 - 118767*x^4 - 118767*x^2 + 187392,
//   F_t = F_0 + t*(g_0 + ... + g_6*x^6).
//
// Primitive vectors g_i in [-Box..Box] are quotiented by sign.  The two
// independent Humbert-4 normal pairings must both be nonzero.  For each
// finite field, projectively congruent directions are cached.  A direct
// point count over F_p and F_{p^2} first tests 72 | #J(F_p); only those
// rare fibers are sent to AbelianGroup for the exact [6,12] test.
//
// Stages:
//   1. all directions at p=5,7;
//   2. the best StageTop at p=11,13;
//   3. the best FinalTop at p=17,19,23.
//
// Finite singular/degree-drop residues are reported as BAD and never
// counted as genuine allowed residues.  The projective parameter infinity
// is recorded separately and is not scored.
//
// Usage:
//   magma -b Box:=2 StageTop:=5000 FinalTop:=250 \
//       output_file:="data/m612_hlp_direction_search_b2.txt" \
//       code/m612_hlp_direction_search.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Z := Integers();

if not assigned Box then Box := 2;
elif Type(Box) eq MonStgElt then Box := StringToInteger(Box); end if;
if not assigned StageTop then StageTop := 5000;
elif Type(StageTop) eq MonStgElt then StageTop := StringToInteger(StageTop); end if;
if not assigned FinalTop then FinalTop := 250;
elif Type(FinalTop) eq MonStgElt then FinalTop := StringToInteger(FinalTop); end if;
if not assigned ReportTop then ReportTop := 20;
elif Type(ReportTop) eq MonStgElt then ReportTop := StringToInteger(ReportTop); end if;
if not assigned Progress then Progress := 1000;
elif Type(Progress) eq MonStgElt then Progress := StringToInteger(Progress); end if;

writeOutput := assigned output_file;
out := 0;
if writeOutput then out := Open(output_file,"w"); end if;

procedure Put(s)
    print s;
    if writeOutput then fprintf out,"%o\n",s; end if;
end procedure;

function Has612(inv)
    return #[n:n in inv | (Z!n) mod 6 eq 0] ge 2 and
           #[n:n in inv | (Z!n) mod 12 eq 0] ge 1;
end function;

function FirstNonzero(v)
    for i in [1..#v] do
        if v[i] ne 0 then return i,v[i]; end if;
    end for;
    return 0,0;
end function;

function PrimitiveSignDirection(v)
    if &and[a eq 0:a in v] then return false; end if;
    c := GCD([Abs(Z!a):a in v]);
    if c ne 1 then return false; end if;
    i,a := FirstNonzero(v);
    return a gt 0;
end function;

function IsTransverse(v)
    nminus := 1298*v[2] + 2423*v[4] + 1298*v[6];
    nplus  := -649*v[1] - 3072*v[3] + 3072*v[5] + 649*v[7];
    return nminus ne 0 and nplus ne 0;
end function;

function CountC(f)
    k := BaseRing(Parent(f));
    n := 0;
    for a in k do
        z := Evaluate(f,a);
        if z eq 0 then n +:= 1;
        elif IsSquare(z) then n +:= 2;
        end if;
    end for;
    if Degree(f) eq 5 then
        n +:= 1;
    elif Degree(f) eq 6 and IsSquare(LeadingCoefficient(f)) then
        n +:= 2;
    end if;
    return n;
end function;

// status:  1 exact [6,12]; 0 smooth but killed; -1 bad; -2 unresolved.
// The returned invariants are nonempty only for the rare 72-divisible
// fibers on which AbelianGroup was computed.
function FiberStatus(f,p,k2,P2)
    if Degree(f) notin {5,6} or not IsSquarefree(f) then
        return -1,[],0;
    end if;
    f2 := P2![k2!Coefficient(f,i):i in [0..Degree(f)]];
    n1 := CountC(f);
    n2 := CountC(f2);
    a1 := p+1-n1;
    a2 := (n2-p^2-1+a1^2) div 2;
    nJ := 1-a1+a2-p*a1+p^2;
    if nJ mod 72 ne 0 then return 0,[],nJ; end if;
    try
        A,mp := AbelianGroup(Jacobian(HyperellipticCurve(f)));
        inv := Invariants(A);
        return (Has612(inv) select 1 else 0),inv,nJ;
    catch e
        return -2,[],nJ;
    end try;
end function;

function NormalizedDirection(v,k)
    w := [k!a:a in v];
    i,s := FirstNonzero(w);
    assert i ne 0;
    n := [a/s:a in w];
    key := Sprint([Z!a:a in n]);
    return key,s,n;
end function;

// Return base masks in the parameter u for the normalized direction n.
function BaseLineMask(n,p,k,P,k2,P2,seedInfo)
    x := P.1;
    f0 := 187392*x^6-118767*x^4-118767*x^2+187392;
    gn := &+[n[i+1]*x^i:i in [0..6]];
    allow := []; bad := []; unresolved := []; details := [];
    seedStatus,seedInv,seedNJ := Explode(seedInfo);
    if seedStatus eq 1 then
        Append(~allow,0); Append(~details,<0,seedInv,seedNJ>);
    elif seedStatus eq -1 then Append(~bad,0);
    elif seedStatus eq -2 then Append(~unresolved,0);
    end if;
    for uz in [1..p-1] do
        u := k!uz;
        status,inv,nJ := FiberStatus(f0+u*gn,p,k2,P2);
        z := Z!u;
        if status eq 1 then
            Append(~allow,z);
            Append(~details,<z,inv,nJ>);
        elif status eq -1 then Append(~bad,z);
        elif status eq -2 then Append(~unresolved,z);
        end if;
    end for;
    return <Sort(allow),Sort(bad),Sort(unresolved),details>;
end function;

// Convert a normalized-u mask to the actual parameter t, where u=s*t.
function RescaleMask(mask,s,k)
    return Sort([Z!((k!u)/s):u in mask]);
end function;

RF := recformat<g, masks, bads, unresolved, details, counts, goods, score1, score2,
                score3, score4>;

function NewRecord(g)
    return rec<RF|g:=g,masks:=[],bads:=[],unresolved:=[],details:=[],counts:=[],
                  goods:=[],score1:=0,score2:=0,score3:=0,score4:=0>;
end function;

function RecomputeScore(r,primes)
    // Lexicographic generosity score:
    //   score1 = number of primes with a genuine nonzero allowed residue;
    //   score2 = total number of genuine nonzero allowed residues;
    //   score3 = sum of nonzero-allowed densities, scaled by 10^6;
    //   score4 = negative total number of bad/unresolved finite residues.
    // t=0 is the split seed (and is actually a bad fiber at p=5), so it
    // never contributes to the nonsplit generosity score.
    extras := [#[t:t in r`masks[i] | t ne 0]:i in [1..#primes]];
    nongood := [#[t:t in (r`bads[i] cat r`unresolved[i]) | t ne 0]
                :i in [1..#primes]];
    denoms := [(primes[i]-1)-nongood[i]:i in [1..#primes]];
    r`score1 := #[n:n in extras | n gt 0];
    r`score2 := &+extras;
    r`score3 := &+[(denoms[i] eq 0) select 0 else
                   ((1000000*extras[i]) div denoms[i])
                   :i in [1..#primes]];
    r`score4 := -&+[#r`bads[i]+#r`unresolved[i]:i in [1..#primes]];
    return r;
end function;

function Better(a,b)
    // Return true when a should precede b.
    if a`score1 ne b`score1 then return a`score1 gt b`score1; end if;
    if a`score2 ne b`score2 then return a`score2 gt b`score2; end if;
    if a`score3 ne b`score3 then return a`score3 gt b`score3; end if;
    if a`score4 ne b`score4 then return a`score4 gt b`score4; end if;
    return Sprint(a`g) lt Sprint(b`g);
end function;

function SortRecords(R)
    Sort(~R,func<a,b| Better(a,b) select -1 else
                         (Better(b,a) select 1 else 0)>);
    return R;
end function;

function AddPrimeData(records,p)
    k := GF(p); P<x> := PolynomialRing(k);
    k2 := GF(p^2); P2<X> := PolynomialRing(k2);
    cache := AssociativeArray();
    hits := 0; misses := 0; rare := 0; unresolvedTotal := 0;
    fseed := 187392*x^6-118767*x^4-118767*x^2+187392;
    seedStatus,seedInv,seedNJ := FiberStatus(fseed,p,k2,P2);
    Put(Sprintf(" SEED p=%o status=%o invariants=%o order=%o",
                p,seedStatus,seedInv,seedNJ));
    seedInfo := <seedStatus,seedInv,seedNJ>;
    for j in [1..#records] do
        key,s,n := NormalizedDirection(records[j]`g,k);
        if not IsDefined(cache,key) then
            cache[key] := BaseLineMask(n,p,k,P,k2,P2,seedInfo);
            misses +:= 1;
        else
            hits +:= 1;
        end if;
        entry := cache[key];
        amask := RescaleMask(entry[1],s,k);
        bmask := RescaleMask(entry[2],s,k);
        umask := RescaleMask(entry[3],s,k);
        details := [<Z!((k!d[1])/s),d[2],d[3]>:d in entry[4]];
        Sort(~details,func<a,b| a[1] lt b[1] select -1 else
                                (a[1] gt b[1] select 1 else 0)>);
        r := records[j];
        Append(~r`masks,amask);
        Append(~r`bads,bmask);
        Append(~r`unresolved,umask);
        Append(~r`details,details);
        Append(~r`counts,#amask);
        Append(~r`goods,p-#bmask-#umask);
        records[j] := r;
        rare +:= #amask;
        unresolvedTotal +:= #umask;
        if Progress gt 0 and j mod Progress eq 0 then
            Put(Sprintf(" PRIME_PROGRESS p=%o records=%o/%o cache=%o",p,j,#records,misses));
        end if;
    end for;
    Put(Sprintf(" PRIME_DONE p=%o records=%o unique_projective=%o cache_hits=%o total_allowed_incidence=%o unresolved=%o",
                p,#records,misses,hits,rare,unresolvedTotal));
    return records;
end function;

function ScoreAndSort(records,primes)
    for i in [1..#records] do records[i] := RecomputeScore(records[i],primes); end for;
    return SortRecords(records);
end function;

procedure ReportRecords(label,records,primes,nreport)
    Put(Sprintf("REPORT %o primes=%o records=%o showing=%o",label,primes,#records,
                Min(nreport,#records)));
    for j in [1..Min(nreport,#records)] do
        r := records[j];
        Put(Sprintf("RANK %o G %o transverse_pairings <%o,%o> score <%o,%o,%o,%o>",
                    j,r`g,
                    1298*r`g[2]+2423*r`g[4]+1298*r`g[6],
                    -649*r`g[1]-3072*r`g[3]+3072*r`g[5]+649*r`g[7],
                    r`score1,r`score2,r`score3,r`score4));
        for i in [1..#primes] do
            Put(Sprintf(" MASK p=%o good=%o allowed=%o bad=%o unresolved=%o infinity=SEPARATE",
                        primes[i],r`goods[i],r`masks[i],r`bads[i],r`unresolved[i]));
            Put(Sprintf("  ALLOWED_RECORDS %o",r`details[i]));
        end for;
    end for;
end procedure;

Put(Sprintf("M612_HLP_DIRECTION_SEARCH Box=%o StageTop=%o FinalTop=%o",Box,StageTop,FinalTop));

vals := [-Box..Box];
dirs := [];
if assigned candidate_file then
    seenDirs := {};
    for line in Split(Read(candidate_file),"\n") do
        if #line eq 0 or line[1] eq "#" then continue; end if;
        g := [StringToInteger(s):s in Split(line,",")];
        assert #g eq 7 and PrimitiveSignDirection(g) and IsTransverse(g);
        key := Sprint(g);
        if key notin seenDirs then Include(~seenDirs,key); Append(~dirs,g); end if;
    end for;
else
    for g0,g1,g2,g3,g4,g5,g6 in vals do
        g := [g0,g1,g2,g3,g4,g5,g6];
        if PrimitiveSignDirection(g) and IsTransverse(g) then Append(~dirs,g); end if;
    end for;
end if;
Put(Sprintf("DIRECTIONS primitive_sign_transverse=%o",#dirs));

records := [NewRecord(g):g in dirs];
primes := [];
for p in [5,7] do
    records := AddPrimeData(records,p);
    Append(~primes,p);
end for;
records := ScoreAndSort(records,primes);
ReportRecords("STAGE1",records,primes,Min(ReportTop,50));
if #records gt StageTop then records := records[1..StageTop]; end if;

for p in [11,13] do
    records := AddPrimeData(records,p);
    Append(~primes,p);
end for;
records := ScoreAndSort(records,primes);
ReportRecords("STAGE2",records,primes,Min(ReportTop,50));
if #records gt FinalTop then records := records[1..FinalTop]; end if;

for p in [17,19,23] do
    records := AddPrimeData(records,p);
    Append(~primes,p);
end for;
records := ScoreAndSort(records,primes);
ReportRecords("FINAL",records,primes,ReportTop);

// Recompute the original delta F=1+x on the identical seven-prime test,
// whether or not it survived the staged ranking, for a fair baseline.
baseline := [NewRecord([1,1,0,0,0,0,0])];
basePrimes := [];
for p in [5,7,11,13,17,19,23] do
    baseline := AddPrimeData(baseline,p);
    Append(~basePrimes,p);
end for;
baseline := ScoreAndSort(baseline,basePrimes);
ReportRecords("BASELINE_1_PLUS_X",baseline,basePrimes,1);

Put("M612_HLP_DIRECTION_SEARCH_DONE");
if writeOutput then delete out; end if;
quit;
