// lane_2266_sigma_fib.m — per-fiber Mordell-Weil attack on one non-identity
// sigma-surface (2026-08-14 plan, section 1, steps 2-3).
// Fibration (uniform in sigma): condition C2 (pair {1,3}) always has t-side
// class 2(t-3); with cB := A(sg[1],sg[3])(u0) and d := 2*cB it reads
// 2(t-3)*cB = s^2, i.e. t = 3 + s^2/d.  Folding condition C_Org gives the
// genus-1 quartic fiber over u0:
//   Org=1: Y^2 = cA*(s^2+6d)*(s^2-2d),   cA := A(sg[1],sg[2])(u0)   (test C3)
//   Org=3: Y^2 = -cA*(s^2+2d)*(s^2-6d),  cA := A(sg[2],sg[3])(u0)   (test C1)
// Fibers are populated by (a) the in-process 2-of-3 hash join at height H
// (drop = the test-condition index) and (b) a direct Points sweep over u0 of
// height <= HU, plus (c) fibers with IsSquare(cA) (rational points at
// infinity - elliptic even without an affine point).  Per fiber: RankBounds
// FIRST; unequal/error -> MWSKIP + Points(C : Bound := PtsBound) fallback;
// rank 0 -> torsion only; conclusive rank >= 1 -> MordellWeilGroup +
// coefficient box.  Every t found is tested on the third condition exactly;
// survivors are funneled ([2,2,6,6] fires HIT).
// sigma index map: si2=[2,1,3]=(12)  si3=[3,2,1]=(13)  si4=[1,3,2]=(23)
//                  si5=[2,3,1]=(123) si6=[3,1,2]=(132)
// Usage: cd product/code &&
//   magma -b Sig:=2 lane_2266_sigma_fib.m > ../logs/sigfib_si2_o1.log
// params: Sig:=<2..6>  Org:=<1|3> (default 1)  H:=150  HU:=40
//   SweepBound:=1000  PtsBound:=100000  NB:=10  MaxFunnel:=60  MemGB:=5
//   Skip:=0 (resume: skip the first Skip fibers)  AlarmS:=9000 (hard cap)
SetColumns(0);
if not assigned Sig then Sig := 2; elif Type(Sig) eq MonStgElt then Sig := StringToInteger(Sig); end if;
if not assigned Org then Org := 1; elif Type(Org) eq MonStgElt then Org := StringToInteger(Org); end if;
if not assigned H then H := 150; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned HU then HU := 40; elif Type(HU) eq MonStgElt then HU := StringToInteger(HU); end if;
if not assigned SweepBound then SweepBound := 1000; elif Type(SweepBound) eq MonStgElt then SweepBound := StringToInteger(SweepBound); end if;
if not assigned PtsBound then PtsBound := 100000; elif Type(PtsBound) eq MonStgElt then PtsBound := StringToInteger(PtsBound); end if;
if not assigned NB then NB := 10; elif Type(NB) eq MonStgElt then NB := StringToInteger(NB); end if;
if not assigned MaxFunnel then MaxFunnel := 60; elif Type(MaxFunnel) eq MonStgElt then MaxFunnel := StringToInteger(MaxFunnel); end if;
if not assigned MemGB then MemGB := 5; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned Skip then Skip := 0; elif Type(Skip) eq MonStgElt then Skip := StringToInteger(Skip); end if;
if not assigned MaxFib then MaxFib := 0; elif Type(MaxFib) eq MonStgElt then MaxFib := StringToInteger(MaxFib); end if;
if not assigned CovBound then CovBound := 3000; elif Type(CovBound) eq MonStgElt then CovBound := StringToInteger(CovBound); end if;
if not assigned TimeCap then TimeCap := 7200; elif Type(TimeCap) eq MonStgElt then TimeCap := StringToInteger(TimeCap); end if;
if not assigned AlarmS then AlarmS := 9000; elif Type(AlarmS) eq MonStgElt then AlarmS := StringToInteger(AlarmS); end if;
error if Sig lt 2 or Sig gt 6, "Sig must be 2..6";
error if Org ne 1 and Org ne 3, "Org must be 1 or 3";
SetMemoryLimit(MemGB*10^9);
SetClassGroupBounds("GRH");
Alarm(AlarmS);
load "split_lab.m";  // run from product/code/

RQ := Rationals();

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

// ordered pair-class VALUE at x (not the kernel)
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
// condition indices: C2 is the parametrized one; Org is folded; the third is tested
Test := Org eq 1 select 3 else 1;
pOrg := PAIRS[Org]; pTest := PAIRS[Test];
prog := Sprintf("../logs/sigmaF%oo%o.progress", Sig, Org);
System(Sprintf("rm -f %o", prog));
printf "FIBER LANE Sig=%o sg=%o Org=%o Test=%o H=%o HU=%o\n", Sig, sg, Org, Test, H, HU;

// quartic fiber over u0; returns ok, quart(s), d, cA, cTest
function FiberData(u0, sg, Org)
    d  := 2*AV(sg[1], sg[3], u0);
    cA := Org eq 1 select AV(sg[1], sg[2], u0) else AV(sg[2], sg[3], u0);
    cT := Org eq 1 select AV(sg[2], sg[3], u0) else AV(sg[1], sg[2], u0);
    if d eq 0 or cA eq 0 or cT eq 0 then return false, RQx!0, d, cA, cT; end if;
    s := RQx.1;
    quart := Org eq 1 select cA*(s^2+6*d)*(s^2-2*d) else -cA*(s^2+2*d)*(s^2-6*d);
    return true, quart, d, cA, cT;
end function;

// ---- populate fibers ----
// (a) in-process join at height H for this sigma, drop = Test (keep C2 & C_Org)
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
fibpts := AssociativeArray();   // u0 -> set of t-values with (C2 & C_Org)
njoin := 0;
for ti in [1..NV] do
    key := < kk[ti][keep[1]], kk[ti][keep[2]] >;
    if not IsDefined(M, key) then continue; end if;
    for ui in M[key] do
        if ui eq ti then continue; end if;
        njoin +:= 1;
        u0 := vals[ui]; tv := vals[ti];
        if IsDefined(fibpts, u0) then Include(~fibpts[u0], tv); else fibpts[u0] := {tv}; end if;
    end for;
end for;
printf "JOIN Sig=%o drop=%o matches=%o fibers=%o\n", Sig, Test, njoin, #Keys(fibpts);

// (b) direct sweep over small u0: Points on the quartic + points at infinity
uSweep := [ v : v in HeightRats(HU) | not v in EXCL ];
nswF := 0; nswP := 0; ninfF := 0;
t0c := Cputime();
for u0 in uSweep do
    ok, quart, d, cA, cT := FiberData(u0, sg, Org);
    if not ok then continue; end if;
    newt := {};
    if IsSquare(cA) then ninfF +:= 1; if not IsDefined(fibpts, u0) then fibpts[u0] := {}; end if; end if;
    C := HyperellipticCurve(quart);
    pts := [];
    try pts := Points(C : Bound := SweepBound); catch e pts := []; end try;
    for P in pts do
        if P[3] eq 0 then continue; end if;
        sv := P[1]/P[3];
        tv := 3 + sv^2/d;
        if tv in EXCL then continue; end if;
        Include(~newt, tv);
    end for;
    if #newt gt 0 then
        nswF +:= 1; nswP +:= #newt;
        if IsDefined(fibpts, u0) then fibpts[u0] join:= newt; else fibpts[u0] := newt; end if;
    end if;
end for;
printf "SWEEP HU=%o fibers-with-points=%o pts=%o infsq-fibers=%o %o s\n", HU, nswF, nswP, ninfF, Cputime()-t0c;

// multi-point fibers first (rank hints), then by ascending height
fibus := Sort(Setseq(Keys(fibpts)),
    func<a,b | ma eq mb select HeightOf(a) - HeightOf(b) else mb - ma
               where ma := #fibpts[a] ge 2 select 1 else 0
               where mb := #fibpts[b] ge 2 select 1 else 0>);
nmulti := #[ u : u in fibus | #fibpts[u] ge 2 ];
printf "FIBERS Sig=%o total=%o multi=%o\n", Sig, #fibus, nmulti;
System(Sprintf("echo 'sigF%oo%o POPULATED %o fibers (%o multi)' >> %o", Sig, Org, #fibus, nmulti, prog));

// ---- per-fiber MW ----
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

nproc := 0; nmwskip := 0; nrank0 := 0; nrankpos := 0; nsurv := 0; nf := 0; ndeck := 0;
t0 := Cputime();
tR0 := Realtime();
funneled := {};

for fi in [1..#fibus] do
    if fi le Skip then continue; end if;
    if Realtime()-tR0 gt TimeCap then printf "TIMECAP %o s at fiber %o/%o\n", TimeCap, fi, #fibus; break; end if;
    u0 := fibus[fi];
    System(Sprintf("echo 'sigF%oo%o FIB %o/%o u=%o skip=%o surv=%o %os' >> %o",
        Sig, Org, fi, #fibus, u0, nmwskip, nsurv, Round(Cputime()-t0), prog));
    ok, quart, d, cA, cT := FiberData(u0, sg, Org);
    if not ok then continue; end if;
    if MaxFib gt 0 and nproc ge MaxFib then printf "MAXFIB %o reached at fiber %o\n", MaxFib, fi; break; end if;
    nproc +:= 1;
    C := HyperellipticCurve(quart);
    // gather all t-values to test: start from the known ones
    tset := fibpts[u0];
    // choose a base point: prefer an affine known point, else infinity
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
    // discipline (2026-08-14, after the u=-1 probe): ALWAYS MinimalModel before
    // any descent (RankBounds on the raw glued model stalled); NEVER call
    // unbounded MordellWeilGroup - rank 0 needs only TorsionSubgroup, rank >= 1
    // uses bounded TwoDescent-cover point searches (index-in-MW risk accepted).
    mwdone := false;
    if havebase then
        okE := true; E := 0; mE := 0; Em := 0; phi := 0;
        try
            E, mE := EllipticCurve(C, P0);
            Em, phi := MinimalModel(E);
        catch e okE := false; end try;
        if okE then
            // map Em-point back to a t-value; returns ok, tv
            phiI := Inverse(phi); mEi := Inverse(mE);
            rlo := -1; rhi := -2;
            try rlo, rhi := RankBounds(Em); catch e rlo := -1; rhi := -2; end try;
            if rlo lt 0 or rlo ne rhi then
                printf "MWSKIP u=%o bounds=[%o,%o]\n", u0, rlo, rhi;
                nmwskip +:= 1;
            else
                // collect Em-points: torsion + (if rank >= 1) bounded cover search
                empts := [];
                okT := true;
                try
                    T, mT := TorsionSubgroup(Em);
                    for g in T do Append(~empts, mT(g)); end for;
                catch e okT := false; end try;
                gens := [];
                genfail := false;
                if rlo ge 1 then
                    okD := true; S := []; mps := [];
                    try S, mps := TwoDescent(Em); catch e okD := false; end try;
                    if okD then
                        for k in [1..#S] do
                            if #gens ge 2 then break; end if;
                            cpts := [];
                            try cpts := Points(S[k] : Bound := CovBound); catch e cpts := []; end try;
                            for cp in cpts do
                                okm := true; EP := 0;
                                try EP := mps[k](cp); catch e okm := false; end try;
                                if not okm then continue; end if;
                                if Order(EP) ne 0 then continue; end if;
                                if #gens ge 1 and (#gens ge 2 or EP eq gens[1] or EP eq -gens[1]) then continue; end if;
                                Append(~gens, EP);
                                if #gens ge 2 then break; end if;
                            end for;
                        end for;
                    end if;
                    if #gens eq 0 then
                        genfail := true;
                        printf "GENFAIL u=%o r=%o (no cover point <= %o)\n", u0, rlo, CovBound;
                        nmwskip +:= 1;
                    end if;
                end if;
                if okT and not genfail then
                    if rlo eq 0 then nrank0 +:= 1;
                    else nrankpos +:= 1; printf "MWRANK u=%o r=%o gens=%o\n", u0, rlo, #gens;
                    end if;
                    mwdone := true;
                    // enumerate box: n1*gens[1] + n2*gens[2] + torsion
                    boxes := [ [Integers()|] ];
                    for gi in [1..#gens] do
                        boxes := [ b cat [n] : n in [-NB..NB], b in boxes ];
                    end for;
                    for b in boxes do
                        base := Em!0;
                        for gi in [1..#gens] do base +:= b[gi]*gens[gi]; end for;
                        for tp in empts do
                            Q := base + tp;
                            okq := true; P := 0;
                            try P := mEi(phiI(Q)); catch e okq := false; end try;
                            if not okq or P[3] eq 0 then continue; end if;
                            sv := P[1]/P[3];
                            tv := 3 + sv^2/d;
                            if tv in EXCL or HeightOf(tv) gt 10^24 then continue; end if;
                            Include(~tset, tv);
                        end for;
                    end for;
                end if;
            end if;
        else
            printf "MWSKIP u=%o Ebuildfail\n", u0;
            nmwskip +:= 1;
        end if;
    end if;
    if not mwdone then
        // fallback: deeper direct point search on the quartic
        pts := [];
        try pts := Points(C : Bound := PtsBound); catch e pts := []; end try;
        for P in pts do
            if P[3] eq 0 then continue; end if;
            sv := P[1]/P[3];
            tv := 3 + sv^2/d;
            if tv in EXCL or HeightOf(tv) gt 10^24 then continue; end if;
            Include(~tset, tv);
        end for;
    end if;
    // test the third condition at every gathered t
    for tv in tset do
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
        printf "INSTANCE 2266sigfib t=%o u=%o si=%o\n", tv, u0, Sig;
        GlueFunnel(EA, EB, Sprintf("l2266sf|si=%o|t=%o|u=%o", Sig, tv, u0), [Integers()|3,3], ~seen, ~nstat);
    end for;
end for;

printf "SEARCH_DONE 2266sigfib Sig=%o Org=%o fibers=%o proc=%o mwskip=%o rank0=%o rankpos=%o surv=%o deck=%o funneled=%o curves=%o aborts=%o known=%o hits=%o fails=%o %o s\n",
    Sig, Org, #fibus, nproc, nmwskip, nrank0, nrankpos, nsurv, ndeck, nf, #seen,
    nstat["abort"], nstat["known"], nstat["hit"], nstat["fail"], Cputime()-t0;
System(Sprintf("echo 'sigF%oo%o DONE %o fibers %o surv %o funneled' >> %o", Sig, Org, #fibus, nsurv, nf, prog));
quit;
