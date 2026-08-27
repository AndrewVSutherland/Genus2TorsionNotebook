// lane_2266_sigma_fib2.m — COMPLETION pass for the rank >= 2 fibers of the
// 2026-08-14 org-A shards (codex round 3, PR #16): the original
// lane_2266_sigma_fib.m capped the cover-point search at TWO generators and
// never verified their independence, so rank-3/4 fiber boxes certainly
// missed directions and rank-2 boxes were unverified.  This lane reprocesses
// exactly the 6444 rank>=2 fibers (u-lists parsed from the committed shard
// logs into ../data/fib2_usets.m): per fiber it gathers up to CovMax
// nontorsion points from ALL TwoDescent covers, greedily selects a
// VERIFIED-INDEPENDENT subset (height-pairing regulator > RegEps) of size
// up to min(rlo, 4), enumerates the full coefficient box over that set
// (NB 10/10/5/4 for 1/2/3/4 generators) plus torsion translates, tests the
// third square condition at every point, and funnels survivors.  Fibers
// where the independent set stays smaller than rlo are flagged GENSHORT
// (counted, listed - honestly incomplete rather than silently truncated).
// sigma index map: si2=[2,1,3]=(12)  si3=[3,2,1]=(13)  si4=[1,3,2]=(23)
//                  si5=[2,3,1]=(123) si6=[3,1,2]=(132)
// Usage: cd product/code &&
//   magma -b Sig:=2 lane_2266_sigma_fib2.m > ../logs/sigfib2_si2.log
// params: Sig:=<2..6>  H:=150  CovBound:=3000  CovMax:=12  MemGB:=5
//   Skip:=0  TimeCap:=7200  RegEps (1e-6 default, via RegEpsE:=6)
SetColumns(0);
if not assigned Sig then Sig := 2; elif Type(Sig) eq MonStgElt then Sig := StringToInteger(Sig); end if;
if not assigned H then H := 150; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned CovBound then CovBound := 3000; elif Type(CovBound) eq MonStgElt then CovBound := StringToInteger(CovBound); end if;
if not assigned CovMax then CovMax := 12; elif Type(CovMax) eq MonStgElt then CovMax := StringToInteger(CovMax); end if;
if not assigned MemGB then MemGB := 5; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned Skip then Skip := 0; elif Type(Skip) eq MonStgElt then Skip := StringToInteger(Skip); end if;
if not assigned TimeCap then TimeCap := 7200; elif Type(TimeCap) eq MonStgElt then TimeCap := StringToInteger(TimeCap); end if;
if not assigned RegEpsE then RegEpsE := 6; elif Type(RegEpsE) eq MonStgElt then RegEpsE := StringToInteger(RegEpsE); end if;
if not assigned MaxFunnel then MaxFunnel := 60; elif Type(MaxFunnel) eq MonStgElt then MaxFunnel := StringToInteger(MaxFunnel); end if;
if not assigned SweepB then SweepB := 1000; elif Type(SweepB) eq MonStgElt then SweepB := StringToInteger(SweepB); end if;
if not assigned UFile then UFile := ""; end if;
if not assigned SatB then SatB := 100; elif Type(SatB) eq MonStgElt then SatB := StringToInteger(SatB); end if;
error if Sig lt 2 or Sig gt 6, "Sig must be 2..6";
SetMemoryLimit(MemGB*10^9);
SetClassGroupBounds("GRH");
load "split_lab.m";       // run from product/code/
load "../data/fib2_usets.m";

RQ := Rationals();
Org := 1;   // org-A throughout, as in the original shards

function E26(t)
    r1 := (-2*t+10)/((t+3)*(t-3));
    r2 := (-t^3+7*t^2-11*t+5)/(4*(t+3)*(t-3)^2);
    r3 := (-2*t^2+4*t-2)/((t+3)^2*(t-3));
    f := (RQx.1 - r1)*(RQx.1 - r2)*(RQx.1 - r3);
    return EllipticCurve(RQx!f);
end function;
ok26 := 0;
for tv in [RQ| 7, 9/2, -4/3 ] do
    try if Invariants(TorsionSubgroup(E26(tv))) eq [2,6] then ok26 +:= 1; end if; catch e; end try;
end for;
printf "VERIFY E26: %o/3\n", ok26;
error if ok26 lt 2, "universal curve reconstruction failed";

function SFrat(x)
    n := Numerator(x)*Denominator(x);
    s := Sign(n); n := Abs(n);
    a := SquarefreeFactorization(n);
    return s*a;
end function;
function HeightRats(H)
    S := [];
    for a in [1..H], b in [1..H] do
        if GCD(a,b) eq 1 then Append(~S, a/b); Append(~S, -a/b); end if;
    end for;
    return S;
end function;
function KO(k, a, b)
    if a lt b then
        return a eq 1 select (b eq 2 select k[1] else k[2]) else k[3];
    else
        return -(b eq 1 select (a eq 2 select k[1] else k[2]) else k[3]);
    end if;
end function;
function AV(a, b, x)
    if a lt b then
        if a eq 1 then
            return b eq 2 select (x+3)*(x-5) else 2*(x-3);
        else
            return -(x-1)*(x-9);
        end if;
    else
        return -AV(b, a, x);
    end if;
end function;
SIGMAS := [ [1,2,3],[2,1,3],[3,2,1],[1,3,2],[2,3,1],[3,1,2] ];
PAIRS  := [ [1,2], [1,3], [2,3] ];
EXCL := {RQ|3,-3,1,5,9};
HeightOf := func<q | Max(Abs(Numerator(q)), Abs(Denominator(q)))>;

sg := SIGMAS[Sig];
Test := 3; pTest := PAIRS[Test];
prog := Sprintf("../logs/sigmaF%ob.progress", Sig);
System(Sprintf("rm -f %o", prog));
if UFile ne "" then USETS := eval Read(UFile); end if;   // override u-lists (e.g. the NOBASE redo)
USET := USETS[Sig-1];
printf "FIB2 Sig=%o sg=%o fibers=%o ufile=%o\n", Sig, sg, #USET, UFile eq "" select "default" else UFile;

function FiberData(u0, sg)
    d  := 2*AV(sg[1], sg[3], u0);
    cA := AV(sg[1], sg[2], u0);
    cT := AV(sg[2], sg[3], u0);
    if d eq 0 or cA eq 0 or cT eq 0 then return false, RQx!0, d, cA, cT; end if;
    s := RQx.1;
    quart := cA*(s^2+6*d)*(s^2-2*d);
    return true, quart, d, cA, cT;
end function;

// join population (same as the original shard: drop = Test)
vals := [ v : v in HeightRats(H) | not v in EXCL ];
NV := #vals;
kk := [];
for v in vals do
    Append(~kk, [ SFrat((v+3)*(v-5)), SFrat(2*(v-3)), SFrat(-(v-1)*(v-9)) ]);
end for;
keep := [ j : j in [1..3] | j ne Test ];
M := AssociativeArray();
for ui in [1..NV] do
    key := < KO(kk[ui], sg[PAIRS[j][1]], sg[PAIRS[j][2]]) : j in keep >;
    if IsDefined(M, key) then Append(~M[key], ui); else M[key] := [ui]; end if;
end for;
fibpts := AssociativeArray();
for ti in [1..NV] do
    key := < kk[ti][keep[1]], kk[ti][keep[2]] >;
    if not IsDefined(M, key) then continue; end if;
    for ui in M[key] do
        if ui eq ti then continue; end if;
        u0 := vals[ui]; tv := vals[ti];
        if IsDefined(fibpts, u0) then Include(~fibpts[u0], tv); else fibpts[u0] := {tv}; end if;
    end for;
end for;
printf "JOIN fibers=%o\n", #Keys(fibpts);

seen := {};
nstat := AssociativeArray();
for st in ["skip","abort","known","hit","fail"] do nstat[st] := 0; end for;
procedure GlueFunnel(EA, EB, tag, oddT, ~seen, ~nstat)
    L := [];
    try L := Genus2Elliptic2(EA, EB); catch e L := []; end try;
    for k in [1..#L] do
        gk := "";
        try gk := Sprintf("%o", G2Invariants(L[k])); catch e gk := "bad"; end try;
        if gk eq "bad" or gk in seen then continue; end if;
        Include(~seen, gk);
        st := Funnel(L[k], Sprintf("%o|%o", tag, k) : OddInvs := oddT);
        nstat[st] +:= 1;
    end for;
end procedure;

nproc := 0; nshort := 0; ncomplete := 0; nmwskip := 0; nsurv := 0; ndeck := 0; nf := 0; nnobase := 0;
nsatfix := 0; nsatwarn := 0; nrbopen := 0;
shortlist := [];
funneled := {};
t0 := Cputime(); tR0 := Realtime();
RegEps := RealField(30)!10.0^(-RegEpsE);

for fi in [1..#USET] do
    if fi le Skip then continue; end if;
    if Realtime()-tR0 gt TimeCap then printf "TIMECAP %o s at fiber %o/%o\n", TimeCap, fi, #USET; break; end if;
    u0 := USET[fi];
    System(Sprintf("echo 'sigF%ob FIB %o/%o u=%o short=%o surv=%o %os' >> %o",
        Sig, fi, #USET, u0, nshort, nsurv, Round(Realtime()-tR0), prog));
    ok, quart, d, cA, cT := FiberData(u0, sg);
    if not ok then continue; end if;
    nproc +:= 1;
    C := HyperellipticCurve(quart);
    tset := IsDefined(fibpts, u0) select fibpts[u0] else {RQ|};
    // base point: affine known point, else infinity
    P0 := 0; havebase := false;
    for tv in tset do
        oks, sv := IsSquare((tv-3)*d);
        if not oks then continue; end if;
        oky, Yv := IsSquare(Evaluate(quart, sv));
        if not oky then continue; end if;
        P0 := C![sv, Yv]; havebase := true; break;
    end for;
    if not havebase then
        ptsinf := PointsAtInfinity(C);
        if #ptsinf gt 0 then P0 := Rep(ptsinf); havebase := true; end if;
    end if;
    if not havebase then
        // fiber was populated by the ORIGINAL direct sweep, not the join:
        // recover an affine base point the same way (codex round 4)
        spts := [];
        try spts := Points(C : Bound := SweepB); catch e spts := []; end try;
        for P in spts do
            if P[3] ne 0 then P0 := P; havebase := true; break; end if;
        end for;
    end if;
    if not havebase then printf "NOBASE u=%o\n", u0; nnobase +:= 1; continue; end if;
    okE := true; E := 0; mE := 0; Em := 0; phi := 0;
    try
        E, mE := EllipticCurve(C, P0);
        Em, phi := MinimalModel(E);
    catch e okE := false; end try;
    if not okE then printf "MWSKIP u=%o Ebuildfail\n", u0; nmwskip +:= 1; continue; end if;
    phiI := Inverse(phi); mEi := Inverse(mE);
    rlo := -1; rhi := -2;
    try rlo, rhi := RankBounds(Em); catch e rlo := -1; rhi := -2; end try;
    if rlo lt 0 then
        printf "MWSKIP u=%o RBfail\n", u0;
        nmwskip +:= 1; continue;
    end if;
    // unequal bounds: attempt to CONCLUDE rank = rhi by exhibiting rhi
    // independent points; if the cover search delivers them the rank is
    // pinned exactly, else the fiber is honestly open (codex round 5)
    target0 := rlo eq rhi select rlo else rhi;
    // torsion
    empts := [];
    okT := true;
    try
        T, mT := TorsionSubgroup(Em);
        for g in T do Append(~empts, mT(g)); end for;
    catch e okT := false; end try;
    if not okT then printf "MWSKIP u=%o torsfail\n", u0; nmwskip +:= 1; continue; end if;
    // gather candidate points from ALL covers
    cands := [];
    okD := true; S := []; mps := [];
    try S, mps := TwoDescent(Em); catch e okD := false; end try;
    if okD then
        for k in [1..#S] do
            if #cands ge CovMax then break; end if;
            cpts := [];
            try cpts := Points(S[k] : Bound := CovBound); catch e cpts := []; end try;
            for cp in cpts do
                if #cands ge CovMax then break; end if;
                okm := true; EP := 0;
                try EP := mps[k](cp); catch e okm := false; end try;
                if not okm then continue; end if;
                if Order(EP) ne 0 then continue; end if;
                if exists{ Q : Q in cands | Q eq EP or Q eq -EP } then continue; end if;
                Append(~cands, EP);
            end for;
        end for;
    end if;
    // greedy verified-independent subset up to min(target0, 4)
    target := Min(target0, 4);
    gens := [];
    for P in cands do
        if #gens ge target then break; end if;
        trial := gens cat [P];
        okp := true; det := RealField(30)!0;
        try
            HM := HeightPairingMatrix(trial);
            det := Determinant(HM);
        catch e okp := false; end try;
        if okp and det gt RegEps then Append(~gens, P); end if;
    end for;
    // SATURATE the independent set (codex round 5; adoption logic fixed in
    // the round-6 follow-up: the old "#gsat eq #gens" guard silently
    // discarded saturation results whenever Saturation returned extra
    // (e.g. torsion) points): the box over a finite-index sublattice
    // misses whole cosets, so any enlargement must be ADOPTED before the
    // box is enumerated.
    if #gens gt 0 then
        oks := true; gsel := gens;
        regb := RealField(30)!0;
        try
            regb := RealField(30)!Determinant(HeightPairingMatrix(gens));
            gsat := Saturation(gens, SatB);
            sel := [];
            for P in gsat do
                if Order(P) ne 0 then continue; end if;
                if #sel ge #gens then break; end if;
                trial := sel cat [P];
                det := RealField(30)!Determinant(HeightPairingMatrix(trial));
                if det gt RealField(30)!10^-8 then Append(~sel, P); end if;
            end for;
            if #sel eq #gens then gsel := sel; else oks := false; end if;
        catch e oks := false; end try;
        if oks then
            rega := RealField(30)!Determinant(HeightPairingMatrix(gsel));
            if regb gt rega*(RealField(30)!3/2) then
                nsatfix +:= 1;
                printf "SATFIX u=%o index~%o (enlargement ADOPTED - box runs on the saturated lattice)\n", u0, Round(Sqrt(regb/rega));
            end if;
            gens := gsel;
        else
            printf "SATWARN u=%o saturation failed - box over unsaturated span\n", u0;
            nsatwarn +:= 1;
        end if;
    end if;
    concl := rlo eq rhi or #gens eq rhi;   // rank pinned exactly?
    if not concl then
        nrbopen +:= 1;
        printf "RBOPEN u=%o bounds=[%o,%o] indep=%o - rank UNRESOLVED, fiber incomplete\n", u0, rlo, rhi, #gens;
    end if;
    if #gens lt Min((concl select (rlo eq rhi select rlo else rhi) else rhi), 4) or not concl then
        if concl then
            nshort +:= 1;
            Append(~shortlist, <u0, target0, #gens>);
            printf "GENSHORT u=%o target=%o indep=%o cands=%o\n", u0, target0, #gens, #cands;
        end if;
    else
        ncomplete +:= 1;
    end if;
    printf "GENS2 u=%o bounds=[%o,%o] cands=%o indep=%o\n", u0, rlo, rhi, #cands, #gens;
    // full box over the verified-independent generators + torsion
    ng := #gens;
    nb := ng le 2 select 10 else (ng eq 3 select 5 else 4);
    boxes := [ [Integers()|] ];
    for gi in [1..ng] do
        boxes := [ b cat [n] : n in [-nb..nb], b in boxes ];
    end for;
    tvals := tset;
    for b in boxes do
        base := Em!0;
        for gi in [1..ng] do base +:= b[gi]*gens[gi]; end for;
        for tp in empts do
            Q := base + tp;
            okq := true; P := 0;
            try P := mEi(phiI(Q)); catch e okq := false; end try;
            if not okq or P[3] eq 0 then continue; end if;
            sv := P[1]/P[3];
            tv := 3 + sv^2/d;
            if tv in EXCL or HeightOf(tv) gt 10^24 then continue; end if;
            Include(~tvals, tv);
        end for;
    end for;
    for tv in tvals do
        if not IsSquare(AV(pTest[1], pTest[2], tv) * cT) then continue; end if;
        nsurv +:= 1;
        printf "FIBSURV Sig=%o t=%o u=%o\n", Sig, tv, u0;
        if HeightOf(tv) gt 10^12 then printf "SKIPBIG t-height\n"; continue; end if;
        pr := [tv, u0];
        if pr in funneled then continue; end if;
        Include(~funneled, pr);
        EA := 0; EB := 0; okc := true;
        try EA := E26(tv); EB := E26(u0); catch e okc := false; end try;
        if not okc then continue; end if;
        if jInvariant(EA) eq jInvariant(EB) then ndeck +:= 1; printf "SKIPISO t=%o u=%o\n", tv, u0; continue; end if;
        nf +:= 1;
        if nf gt MaxFunnel then printf "FUNNEL CAP %o reached\n", MaxFunnel; continue; end if;
        printf "INSTANCE 2266sigfib2 t=%o u=%o si=%o\n", tv, u0, Sig;
        GlueFunnel(EA, EB, Sprintf("l2266sf2|si=%o|t=%o|u=%o", Sig, tv, u0), [Integers()|3,3], ~seen, ~nstat);
    end for;
end for;

printf "GENSHORT-LIST Sig=%o: %o\n", Sig, shortlist;
printf "SEARCH_DONE 2266sigfib2 Sig=%o fibers=%o proc=%o complete=%o genshort=%o nobase=%o rbopen=%o satfix=%o satwarn=%o mwskip=%o surv=%o deck=%o funneled=%o curves=%o aborts=%o known=%o hits=%o fails=%o %o s\n",
    Sig, #USET, nproc, ncomplete, nshort, nnobase, nrbopen, nsatfix, nsatwarn, nmwskip, nsurv, ndeck, nf, #seen,
    nstat["abort"], nstat["known"], nstat["hit"], nstat["fail"], Cputime()-t0;
System(Sprintf("echo 'sigF%ob DONE proc=%o complete=%o short=%o surv=%o' >> %o", Sig, nproc, ncomplete, nshort, nsurv, prog));
quit;
