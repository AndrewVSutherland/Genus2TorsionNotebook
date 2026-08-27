// lane_2266_sigma_lines.m — deep sweeps along the rational curves found ON
// the sigma-surfaces by the empirical section scan (sigma_section_scan.py,
// 2026-08-14): lines and deg<=2/deg<=2 curves t = N(u)/D(u) carrying many
// 2-of-3 surface points.  On such a curve at least one condition typically
// degenerates (square factors drop out), so the full [2,2,6,6] system
// reduces to 1-2 square conditions in the single variable u — sweepable to
// heights far beyond any (t,u)-join.  For each curve: (1) symbolic
// restriction of all three conditions, classes + degrees printed (RESIDUAL);
// (2) exact u-sweep to height HSweep testing the residual conditions;
// (3) any full pass -> guards (EXCL/pole/deck) -> GlueFunnel (HIT = [2,2,6,6]).
// sigma index map: si2=[2,1,3]=(12) si3=[3,2,1]=(13) si4=[1,3,2]=(23)
//                  si5=[2,3,1]=(123) si6=[3,1,2]=(132)
// Usage: cd product/code && magma -b lane_2266_sigma_lines.m > ../logs/lane2266_sigma_lines.log
//   optional: HSweep:=2000 MaxFunnel:=60 MemGB:=4
SetColumns(0);
if not assigned HSweep then HSweep := 2000; elif Type(HSweep) eq MonStgElt then HSweep := StringToInteger(HSweep); end if;
if not assigned MaxFunnel then MaxFunnel := 60; elif Type(MaxFunnel) eq MonStgElt then MaxFunnel := StringToInteger(MaxFunnel); end if;
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
load "split_lab.m";  // run from product/code/

RQ := Rationals();
QU<u> := FunctionField(RQ);
Pu<U> := PolynomialRing(RQ);

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

function SqClassFF(q)
    P := Numerator(q)*Denominator(q);
    P := Pu!P;
    error if P eq 0, "zero class";
    fac := Factorization(P);
    unit := P div &*[ f[1]^f[2] : f in fac ];
    error if Degree(unit) ne 0, "unit not constant";
    c := RQ!unit;
    sc := Sign(Numerator(c)*Denominator(c));
    n := Abs(Numerator(c)*Denominator(c));
    sq := SquarefreeFactorization(n);
    r := Pu!(sc*sq);
    for f in fac do
        if IsOdd(f[2]) then r *:= f[1]; end if;
    end for;
    return r;
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
DECKMAP := [ QU!u, (u+15)/(u-1), 6-u, (5*u-9)/(u-5), (u-21)/(u-5), (5*u-21)/(u-1) ];

function HeightRats(H)
    S := [];
    for a in [1..H], b in [1..H] do
        if GCD(a,b) eq 1 then Append(~S, a/b); Append(~S, -a/b); end if;
    end for;
    return S;
end function;

// curve entries: <si, [n0,n1,n2], [d0,d1,d2]> for t(u) = (n0+n1 u+n2 u^2)/(d0+d1 u+d2 u^2)
// from sigma_section_scan.py summary (support counts in comments)
CURVES := [
    // si3 heavy hitters
    <3, [RQ|4,1,0],   [RQ|1,0,0]>,      // t=u+4      [56]
    <3, [RQ|-4,1,0],  [RQ|1,0,0]>,      // t=u-4      [56]
    <3, [RQ|172,67,6],[RQ|43,6,0]>,     // DEG2       [58]
    <3, [RQ|261,-318,96],[RQ|0,-29,16]>,// DEG2       [50]
    <3, [RQ|1080,-1089,0],[RQ|720,-846,121]>, // DEG2 [49]
    // si3 lambda-ruling representatives (slope -k^2 lines; C2-identical)
    <3, [RQ|15,-4,0], [RQ|1,0,0]>,      // t=-4u+15   [16]
    <3, [RQ|15/4,-1/4,0],[RQ|1,0,0]>,   // t=-u/4+15/4 [16]
    <3, [RQ|58/3,-49/9,0],[RQ|1,0,0]>,  // [13]
    <3, [RQ|174/49,-9/49,0],[RQ|1,0,0]>,// [13]
    <3, [RQ|87/4,-25/4,0],[RQ|1,0,0]>,  // [8]
    <3, [RQ|39/4,-9/4,0],[RQ|1,0,0]>,   // [8]
    <3, [RQ|222/25,-49/25,0],[RQ|1,0,0]>, // [8]
    <3, [RQ|222/49,-25/49,0],[RQ|1,0,0]>, // [8]
    <3, [RQ|13/3,-4/9,0],[RQ|1,0,0]>,   // [8]
    <3, [RQ|87/25,-4/25,0],[RQ|1,0,0]>, // [8]
    <3, [RQ|30,-9,0],[RQ|1,0,0]>,       // [7]
    <3, [RQ|111,-36,0],[RQ|1,0,0]>,     // [7]
    // si3 constants (swapped-fibration high-rank hints)
    <3, [RQ|2,0,0], [RQ|1,0,0]>,        // t=2  [5]
    <3, [RQ|4,0,0], [RQ|1,0,0]>,        // t=4  [5]
    // si2
    <2, [RQ|8/5,7/15,0],[RQ|1,0,0]>,    // [6]
    <2, [RQ|-24/7,15/7,0],[RQ|1,0,0]>,  // [6]
    <2, [RQ|19/3,0,0],[RQ|1,0,0]>,      // t=19/3 [5]
    <2, [RQ|17,0,0], [RQ|1,0,0]>,       // t=17   [5]
    <2, [RQ|-113/24,13/24,0],[RQ|1,0,0]>, // [5]
    <2, [RQ|113/13,24/13,0],[RQ|1,0,0]>,  // [5]
    // si4
    <4, [RQ|8/5,7/15,0],[RQ|1,0,0]>,    // [6]
    <4, [RQ|-24/7,15/7,0],[RQ|1,0,0]>,  // [6]
    <4, [RQ|-235/64,-81/64,0],[RQ|1,0,0]>, // [5]
    <4, [RQ|-235/81,-64/81,0],[RQ|1,0,0]>, // [5]
    <4, [RQ|-11,0,0],[RQ|1,0,0]>,       // t=-11  [5]
    <4, [RQ|-1/3,0,0],[RQ|1,0,0]>,      // t=-1/3 [5]
    // si5
    <5, [RQ|9,438,-1071],[RQ|9,454,-255]>,  // DEG2 [40]
    <5, [RQ|99,336,-1323],[RQ|99,512,-315]>,// DEG2 [39]
    <5, [RQ|66/7,-15/7,0],[RQ|1,0,0]>,  // [6]
    <5, [RQ|22/5,-7/15,0],[RQ|1,0,0]>,  // [6]
    <5, [RQ|257/13,-24/13,0],[RQ|1,0,0]>, // [5]
    <5, [RQ|19/3,0,0],[RQ|1,0,0]>,      // [5]
    <5, [RQ|17,0,0],[RQ|1,0,0]>,        // [5]
    <5, [RQ|1231/245,114/245,0],[RQ|1,0,0]>, // [5]
    // si6
    <6, [RQ|33,59,20],[RQ|-3,-1,4]>,    // DEG2 [54]
    <6, [RQ|517,774,245],[RQ|-47,-2,49]>, // DEG2 [54]
    <6, [RQ|913,1438,465],[RQ|-83,-10,93]>, // DEG2 [54]
    <6, [RQ|66/7,-15/7,0],[RQ|1,0,0]>,  // [6]
    <6, [RQ|22/5,-7/15,0],[RQ|1,0,0]>,  // [6]
    <6, [RQ|257/24,-13/24,0],[RQ|1,0,0]>, // [5]
    <6, [RQ|-11,0,0],[RQ|1,0,0]>,       // [5]
    <6, [RQ|-1/3,0,0],[RQ|1,0,0]>,      // [5]
    <6, [RQ|-1231/114,245/114,0],[RQ|1,0,0]>  // [5]
];

uvals := [ v : v in HeightRats(HSweep) | not v in EXCL ];
printf "LINES sweep HSweep=%o uvals=%o curves=%o\n", HSweep, #uvals, #CURVES;

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

nf := 0;
funneled := {};
t0 := Cputime();
for ci in [1..#CURVES] do
    ent := CURVES[ci];
    si := ent[1]; sg := SIGMAS[si];
    Nc := ent[2]; Dc := ent[3];
    tofu := (Nc[1] + Nc[2]*u + Nc[3]*u^2) / (Dc[1] + Dc[2]*u + Dc[3]*u^2);
    // deck guard: skip curves that ARE deck graphs
    if &or[ tofu eq DECKMAP[k] : k in [1..6] ] then
        printf "CURVE %o si=%o t(u)=%o IS DECK - skipped\n", ci, si, tofu;
        continue;
    end if;
    // symbolic restriction of the three conditions
    clss := [Pu|]; trivial := [];
    okres := true;
    for j in [1..3] do
        p := PAIRS[j];
        okc := true; cls := Pu!1;
        try cls := SqClassFF(AV(p[1], p[2], tofu) * AV(sg[p[1]], sg[p[2]], QU!u)); catch e okc := false; end try;
        if not okc then okres := false; break; end if;
        Append(~clss, cls);
        Append(~trivial, cls eq Pu!1);
    end for;
    if not okres then printf "CURVE %o si=%o t(u)=%o RESTRICTFAIL\n", ci, si, tofu; continue; end if;
    printf "RESIDUAL %o si=%o t(u)=%o trivial=%o degs=%o\n", ci, si, tofu, trivial,
        [ Degree(clss[j]) : j in [1..3] ];
    if &and trivial then
        printf "PARAMETRIC FULL FAMILY curve %o si=%o t(u)=%o <<<<< CHECK\n", ci, si, tofu;
    end if;
    // residual conditions to sweep (nontrivial ones), cheapest first
    resj := [ j : j in [1..3] | not trivial[j] ];
    Sort(~resj, func<a,b | Degree(clss[a]) - Degree(clss[b])>);
    nfull := 0; ntest := 0;
    for uv in uvals do
        // pole / degeneracy guards
        dv := Dc[1] + Dc[2]*uv + Dc[3]*uv^2;
        if dv eq 0 then continue; end if;
        tv := (Nc[1] + Nc[2]*uv + Nc[3]*uv^2) / dv;
        if tv in EXCL then continue; end if;
        okall := true;
        for j in resj do
            w := Evaluate(clss[j], uv);
            if w eq 0 or not IsSquare(w) then okall := false; break; end if;
        end for;
        if not okall then continue; end if;
        // full verification on the ORIGINAL conditions (belt and braces)
        c1 := IsSquare(AV(1,2,tv) * AV(sg[1],sg[2], uv));
        c2 := IsSquare(AV(1,3,tv) * AV(sg[1],sg[3], uv));
        c3 := IsSquare(AV(2,3,tv) * AV(sg[2],sg[3], uv));
        if not (c1 and c2 and c3) then continue; end if;
        nfull +:= 1;
        printf "LINEFULL curve=%o si=%o t=%o u=%o ht_u=%o\n", ci, si, tv, uv, HeightOf(uv);
        if HeightOf(tv) gt 10^12 then printf "SKIPBIG\n"; continue; end if;
        pr := [tv, uv];
        if pr in funneled then continue; end if;
        Include(~funneled, pr);
        EA := 0; EB := 0; okc := true;
        try EA := E26(tv); EB := E26(uv); catch e okc := false; end try;
        if not okc then continue; end if;
        if jInvariant(EA) eq jInvariant(EB) then printf "SKIPISO t=%o u=%o\n", tv, uv; continue; end if;
        nf +:= 1;
        if nf gt MaxFunnel then printf "FUNNEL CAP %o reached\n", MaxFunnel; continue; end if;
        printf "INSTANCE 2266line t=%o u=%o si=%o curve=%o\n", tv, uv, si, ci;
        GlueFunnel(EA, EB, Sprintf("l2266ln|si=%o|t=%o|u=%o", si, tv, uv), [Integers()|3,3], ~seen, ~nstat);
    end for;
    printf "CURVE_DONE %o si=%o full=%o %o s\n", ci, si, nfull, Cputime()-t0;
    System(Sprintf("echo 'lines curve %o/%o si=%o full=%o' >> ../logs/sigmalines.progress", ci, #CURVES, si, nfull));
end for;
printf "SEARCH_DONE 2266lines HSweep=%o funneled=%o curves=%o aborts=%o known=%o hits=%o fails=%o %o s\n",
    HSweep, nf, #seen, nstat["abort"], nstat["known"], nstat["hit"], nstat["fail"], Cputime()-t0;
quit;
