// lane_2266_sigma_sat.m — FULL-SATURATION CERTIFICATE pass (codex round 6
// follow-up, owner-approved): for each certified-complete fiber of the
// saturated completion run, prove the index of the generator lattice in
// E(Q) is exactly 1.
// Logic: Saturation(gens, 100) already excludes any index with a prime
// factor <= 100, so a nontrivial index is >= 101.  If idx >= 101 then the
// Mordell-Weil height lattice has Gram determinant R_MW = R_sub/idx^2 <=
// R_sub/101^2, and Minkowski gives a nontorsion point with
//   hhat <= T := gamma_r * (R_sub/101^2)^(1/r)
// (gamma_1..4 = 1, 2/sqrt(3), 2^(1/3), sqrt(2) — exact Hermite constants).
// CPSHeightBounds(Em) gives h(P) - hhat(P) <= ub (calibrated in
// probe_ht.m: Magma hhat ~ log H(x), CPS bounds h - hhat, Points' Bound
// caps max(|num|,|den|) of x).  So an EXHAUSTIVE search of naive height
// <= T + ub decides:
//   - no nontorsion point with hhat <= T          -> m1 > T -> idx <= 100
//     -> idx = 1.                                              [SATCERT]
//   - points found, m1 := exact minimum hhat; idx <= sqrt(R_sub*gamma^r/m1^r)
//     <= 100 -> idx = 1.                                       [SATCERTM]
//   - else: found points are membership-checked against <gens, tors>;
//     outside -> SATVIOLATION (lattice genuinely larger — investigate);
//     inside but index bound >= 101 -> SATOPENM (uncertified).
//   - search bound exceeds PtsCap -> SATOPENB (uncertified, bound logged).
// Usage: cd product/code &&
//   magma -b Sig:=2 lane_2266_sigma_sat.m > ../logs/sigsat_si2.log
// params: Sig:=<2..6>  H:=150  CovBound:=3000  CovMax:=12  RegEpsE:=6
//   SatB:=100  PtsCap:=10000000  MemGB:=5  Skip:=0  TimeCap:=7200
SetColumns(0);
if not assigned Sig then Sig := 2; elif Type(Sig) eq MonStgElt then Sig := StringToInteger(Sig); end if;
if not assigned H then H := 150; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned CovBound then CovBound := 3000; elif Type(CovBound) eq MonStgElt then CovBound := StringToInteger(CovBound); end if;
if not assigned CovMax then CovMax := 12; elif Type(CovMax) eq MonStgElt then CovMax := StringToInteger(CovMax); end if;
if not assigned RegEpsE then RegEpsE := 6; elif Type(RegEpsE) eq MonStgElt then RegEpsE := StringToInteger(RegEpsE); end if;
if not assigned SatB then SatB := 100; elif Type(SatB) eq MonStgElt then SatB := StringToInteger(SatB); end if;
if not assigned PtsCap then PtsCap := 2000000; elif Type(PtsCap) eq MonStgElt then PtsCap := StringToInteger(PtsCap); end if;
if not assigned PtsCapHard then PtsCapHard := 12000000; elif Type(PtsCapHard) eq MonStgElt then PtsCapHard := StringToInteger(PtsCapHard); end if;
if not assigned SatCap then SatCap := 20000; elif Type(SatCap) eq MonStgElt then SatCap := StringToInteger(SatCap); end if;
if not assigned SatFar then SatFar := 20000; elif Type(SatFar) eq MonStgElt then SatFar := StringToInteger(SatFar); end if;
if not assigned MemGB then MemGB := 5; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned Skip then Skip := 0; elif Type(Skip) eq MonStgElt then Skip := StringToInteger(Skip); end if;
if not assigned TimeCap then TimeCap := 7200; elif Type(TimeCap) eq MonStgElt then TimeCap := StringToInteger(TimeCap); end if;
error if Sig lt 2 or Sig gt 6, "Sig must be 2..6";
SetMemoryLimit(MemGB*10^9);
SetClassGroupBounds("GRH");
load "split_lab.m";       // run from product/code/  (for RQx; funnel unused)

RQ := Rationals();
USETS := eval Read("../data/sat_usets.m");
USET := USETS[Sig-1];

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
sg := SIGMAS[Sig];
Test := 3;
prog := Sprintf("../logs/sigmaF%os.progress", Sig);
System(Sprintf("rm -f %o", prog));
// kill-proof per-fiber results (the -b buffer dies with KILLed segments;
// accounting reads THIS file, appended unbuffered after every fiber)
resf := Sprintf("../logs/sigsatR_si%o.txt", Sig);
procedure Rec(resf, cls, u0)
    System(Sprintf("echo '%o u=%o' >> %o", cls, u0, resf));
end procedure;
printf "SATCERT LANE Sig=%o sg=%o fibers=%o PtsCap=%o\n", Sig, sg, #USET, PtsCap;

function FiberData(u0, sg)
    d  := 2*AV(sg[1], sg[3], u0);
    cA := AV(sg[1], sg[2], u0);
    if d eq 0 or cA eq 0 then return false, RQx!0, d; end if;
    s := RQx.1;
    return true, cA*(s^2+6*d)*(s^2-2*d), d;
end function;

// verified saturation adoption (codex-round-6 follow-up bug fix: the old
// "#gsat eq #gens" guard silently DISCARDED saturation results whenever
// Saturation returned extra (e.g. torsion) points, so enlargements could
// slip through unverified).  Returns ok, newgens, newreg, enlarged.
function SatAdopt(gens, B, Rsub, RR)
    r := #gens;
    ok := true; gs := gens;
    try gs := Saturation(gens, B); catch e ok := false; end try;
    if not ok then return false, gens, Rsub, false; end if;
    sel := [];
    for P in gs do
        if Order(P) ne 0 then continue; end if;
        if #sel ge r then break; end if;
        trial := sel cat [P];
        okp := true; det := RR!0;
        try det := RR!Determinant(HeightPairingMatrix(trial)); catch e okp := false; end try;
        if okp and det gt RR!10^-8 then Append(~sel, P); end if;
    end for;
    if #sel lt r then return false, gens, Rsub, false; end if;
    Rafter := RR!Determinant(HeightPairingMatrix(sel));
    if Rafter lt Rsub * RR!0.9 then
        return true, sel, Rafter, true;    // genuine enlargement adopted
    end if;
    return true, sel, Rafter, false;
end function;

// join population (base points)
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

RR := RealField(30);
GAMMA := [ RR!1, 2/Sqrt(RR!3), Root(RR!2, 3), Sqrt(RR!2) ];

ncert := 0; ncertm := 0; nopenb := 0; nopenm := 0; nviol := 0; nskip := 0; nsatfix := 0;
t0 := Cputime(); tR0 := Realtime();

for fi in [1..#USET] do
    if fi le Skip then continue; end if;
    if Realtime()-tR0 gt TimeCap then printf "TIMECAP %o s at fiber %o/%o\n", TimeCap, fi, #USET; break; end if;
    u0 := USET[fi];
    System(Sprintf("echo 'sigF%os FIB %o/%o u=%o cert=%o openB=%o openM=%o %os' >> %o",
        Sig, fi, #USET, u0, ncert+ncertm, nopenb, nopenm, Round(Realtime()-tR0), prog));
    ok, quart, d := FiberData(u0, sg);
    if not ok then nskip +:= 1; continue; end if;
    C := HyperellipticCurve(quart);
    tset := IsDefined(fibpts, u0) select fibpts[u0] else {RQ|};
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
        spts := [];
        try spts := Points(C : Bound := 1000); catch e spts := []; end try;
        for P in spts do
            if P[3] ne 0 then P0 := P; havebase := true; break; end if;
        end for;
    end if;
    if not havebase then printf "SATSKIP u=%o nobase\n", u0; Rec(resf, "SATSKIP", u0); nskip +:= 1; continue; end if;
    okE := true; E := 0; mE := 0; Em := 0; phi := 0;
    try
        E, mE := EllipticCurve(C, P0);
        Em, phi := MinimalModel(E);
    catch e okE := false; end try;
    if not okE then printf "SATSKIP u=%o Ebuild\n", u0; Rec(resf, "SATSKIP", u0); nskip +:= 1; continue; end if;
    rlo := -1; rhi := -2;
    try rlo, rhi := RankBounds(Em); catch e rlo := -1; rhi := -2; end try;
    // torsion
    tors := [];
    okT := true;
    try
        T0, mT := TorsionSubgroup(Em);
        for g in T0 do Append(~tors, mT(g)); end for;
    catch e okT := false; end try;
    if not okT then printf "SATSKIP u=%o tors\n", u0; Rec(resf, "SATSKIP", u0); nskip +:= 1; continue; end if;
    // gens: covers -> greedy independent -> Saturation(<=100)
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
    target0 := rlo eq rhi select rlo else rhi;
    target := Min(target0, 4);
    gens := [];
    RegEps := RR!10.0^(-RegEpsE);
    for P in cands do
        if #gens ge target then break; end if;
        trial := gens cat [P];
        okp := true; det := RR!0;
        try det := RR!Determinant(HeightPairingMatrix(trial)); catch e okp := false; end try;
        if okp and det gt RegEps then Append(~gens, P); end if;
    end for;
    concl := rlo eq rhi or #gens eq rhi;
    if not concl or #gens lt Min(target0, 4) then
        printf "SATSKIP u=%o not-complete-here (bounds=[%o,%o] gens=%o)\n", u0, rlo, rhi, #gens;
        Rec(resf, "SATSKIP", u0); nskip +:= 1; continue;
    end if;
    r := #gens;
    Rsub := RR!Determinant(HeightPairingMatrix(gens));
    okS, gens, Rsub, enl := SatAdopt(gens, SatB, Rsub, RR);
    if not okS then
        printf "SATSKIP u=%o saturation-at-%o failed\n", u0, SatB;
        Rec(resf, "SATSKIP", u0); nskip +:= 1; continue;
    end if;
    if enl then
        nsatfix +:= 1;
        printf "SATFIX u=%o - saturation at <=%o ENLARGED the lattice (old boxes ran on a sublattice)\n", u0, SatB;
        printf "SATFIXBOX u=%o\n", u0;   // needs a box rerun with corrected gens
        Rec(resf, "SATFIX", u0);
    end if;
    gam := GAMMA[r];
    T := gam * (Rsub/RR!(101^2))^(RR!1/r);
    // upper bound ub for h - hhat: min of CPS and Silverman (probe_ht.m
    // calibration; CPS degrades to ub ~ 10^1..10^2 on some fibers, where
    // Silverman's discriminant-scale bound is saner)
    ub := RR!10^9;
    try
        lb0, ub0 := CPSHeightBounds(Em);
        ub := Min(ub, Max(RR!ub0, RR!0));
    catch e; end try;
    try
        ub := Min(ub, Max(RR!SilvermanBound(Em), RR!0));
    catch e; end try;
    // v2 (codex round 6 follow-up): a CAPPED search is still exhaustive
    // below Teff := log(Bsearch) - ub; if no nontorsion point has
    // hhat <= Teff then m1 > Teff and idx <= sqrt(Rsub*gam^r/Teff^r) =: Bidx,
    // and Saturation(gens, Bidx) then excludes every possible prime factor
    // of idx, proving idx = 1 - no e^(T+ub) search needed.
    // v4 cost balance: Teff trades search cost e^(Teff+ub) against the
    // final Saturation cost at Bidx = sqrt(Rsub*gam^r/Teff^r).  A fixed
    // tiny Teff (v3) made Bidx ~ 10^3-10^4 on rank-3/4 fibers and the
    // closing Saturation ate 30-60 min/fiber; instead aim Bidx ~ 300
    // (sub-second saturation) and let the search absorb the rest:
    Btarget := RR!300;
    TeffN := (Rsub * gam^r / Btarget^2)^(RR!1/r);
    Tgoal := Min(T, Max(RR!1/50, TeffN));
    Bfull := Exp(Min(Tgoal + ub, RR!60));
    Bsearch := Bfull le RR!PtsCap select Ceiling(Bfull) else
               (Bfull le RR!PtsCapHard select Ceiling(Bfull) else PtsCapHard);
    Teff := Min(T, Log(RR!Bsearch) - ub);
    if Teff le RR!0.005 then
        // even the hard cap cannot reach positive exhaustive height:
        // record, and at least push the plain saturation much further
        oksat2 := true;
        try gs2 := Saturation(gens, SatFar); catch e oksat2 := false; end try;
        printf "SATOPENB u=%o r=%o Rsub=%o T=%o ub=%o Bfull=%o (Teff<=0; plain saturation pushed to %o: %o)\n",
            u0, r, Rsub, T, ub, Bfull, SatFar, oksat2 select "ok" else "failed";
        Rec(resf, "SATOPENB", u0); nopenb +:= 1; continue;
    end if;
    // a FAILED exhaustive search must never masquerade as an empty one
    // (codex): the certificate rests on exhaustiveness below Teff
    okP := true;
    spts := [];
    try spts := Points(Em : Bound := Bsearch); catch e okP := false; end try;
    if not okP then
        printf "SATOPENM u=%o exhaustive search FAILED (Points threw at B=%o)\n", u0, Bsearch;
        Rec(resf, "SATOPENM", u0); nopenm +:= 1; continue;
    end if;
    low := [];
    for P in spts do
        if P[3] eq 0 then continue; end if;
        if Order(P) ne 0 then continue; end if;
        hh := RR!Height(P);
        if hh le Teff*(1 + RR!10^-9) + RR!10^-12 then Append(~low, P); end if;
    end for;
    m1lb := Teff;
    if #low gt 0 then
        m1lb := Min([ RR!Height(P) : P in low ]);
        // sanity: a low point outside <gens, tors> would mean a missed
        // generator - loud marker, no certificate
        viol := false;
        HM := HeightPairingMatrix(gens);
        HMi := HM^-1;
        for P in low do
            vec := Matrix(RR, r, 1, [ HeightPairing(P, g) : g in gens ]);
            cf := HMi * vec;
            cint := [ Round(cf[i][1]) : i in [1..r] ];
            if exists{ i : i in [1..r] | Abs(cf[i][1] - cint[i]) gt RR!0.001 } then viol := true;
            else
                Q := Em!0;
                for i in [1..r] do Q +:= cint[i]*gens[i]; end for;
                if not exists{ t : t in tors | Q + t eq P } then viol := true; end if;
            end if;
            if viol then break; end if;
        end for;
        if viol then
            nviol +:= 1;
            printf "SATVIOLATION u=%o - low point OUTSIDE <gens,tors>: INVESTIGATE\n", u0;
            Rec(resf, "SATVIOLATION", u0); continue;
        end if;
    end if;
    Bidx := Floor(Sqrt(Rsub * gam^r / m1lb^r));
    if Bidx le 100 then
        ncert +:= 1;
        printf "SATCERT u=%o r=%o Rsub=%o Teff=%o ub=%o B=%o Bidx=%o (idx=1)\n",
            u0, r, Rsub, Teff, ub, Bsearch, Bidx;
        Rec(resf, "SATCERT", u0);
    elif Bidx le SatCap then
        okS2, gens2, R2, enl2 := SatAdopt(gens, Bidx, Rsub, RR);
        if okS2 and not enl2 then
            ncertm +:= 1;
            printf "SATCERT2 u=%o r=%o m1lb=%o Bidx=%o (saturated to Bidx -> idx=1)\n", u0, r, m1lb, Bidx;
            Rec(resf, "SATCERT2", u0);
        elif okS2 and enl2 then
            // an index in (100, Bidx] was real and is now adopted: the
            // certificate for the ENLARGED lattice needs (and gets) a
            // second round: recompute Bidx from the new regulator and
            // re-saturate; boxes must be rerun on the corrected lattice
            nsatfix +:= 1;
            printf "SATFIX u=%o - saturation at <=%o ENLARGED the lattice\n", u0, Bidx;
            printf "SATFIXBOX u=%o\n", u0;
            Rec(resf, "SATFIX", u0);
            Bidx2 := Floor(Sqrt(R2 * gam^r / m1lb^r));
            okS3, gens3, R3, enl3 := SatAdopt(gens2, Min(Bidx2, SatCap), R2, RR);
            if okS3 and not enl3 and Bidx2 le SatCap then
                ncertm +:= 1;
                printf "SATCERT2 u=%o r=%o m1lb=%o Bidx=%o (post-fix, idx=1)\n", u0, r, m1lb, Bidx2;
                Rec(resf, "SATCERT2", u0);
            else
                nopenm +:= 1;
                printf "SATOPENM u=%o post-fix undecided (Bidx2=%o okS3=%o enl3=%o)\n", u0, Bidx2, okS3, enl3;
                Rec(resf, "SATOPENM", u0);
            end if;
        else
            nopenm +:= 1;
            printf "SATOPENM u=%o r=%o m1lb=%o Bidx=%o (saturation call failed)\n", u0, r, m1lb, Bidx;
            Rec(resf, "SATOPENM", u0);
        end if;
    else
        nopenm +:= 1;
        printf "SATOPENM u=%o r=%o m1lb=%o Bidx=%o > SatCap\n", u0, r, m1lb, Bidx;
        Rec(resf, "SATOPENM", u0);
    end if;
end for;

printf "SEARCH_DONE satcert Sig=%o fibers=%o cert=%o certm=%o openB=%o openM=%o viol=%o skip=%o satfix=%o %o s\n",
    Sig, #USET, ncert, ncertm, nopenb, nopenm, nviol, nskip, nsatfix, Cputime()-t0;
System(Sprintf("echo 'sigF%os DONE cert=%o certm=%o openB=%o openM=%o viol=%o skip=%o' >> %o",
    Sig, ncert, ncertm, nopenb, nopenm, nviol, nskip, prog));
quit;
