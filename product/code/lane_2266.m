// lane_2266.m — [2,2,6,6] sweep (2026-08-13 plan, section 1, step 3).
// E26(t) x E26(u): a pair is funneled iff the FULL two-gain Prop-12 system
// holds for some sigma in S3, i.e. the three unordered-pair class-matching
// conditions (validated symbolically in lane_2266_derive.m; anchors PASSed):
//   pair classes of E26(t) mod squares:
//     A12 = (t+3)(t-5),  A13 = 2(t-3),  A23 = -(t-1)(t-9)
//   condition (sigma): KO(t;i,k) = KO(u;sigma(i),sigma(k)) as signed
//   squarefree kernels, for the three pairs (ordered; transposes negate).
// Any survivor has full rational J[2] (2-rank 4) and odd part [3,3]:
// J(Q)_tors >= [2,2,6,6] (order 144).
// Usage: cd product/code && magma -b lane_2266.m > ../logs/lane2266_H60.log
//   optional: H:=<int> (default 60), MaxFunnel:=<int> (default 200)
SetColumns(0);
if not assigned H then H := 60; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned MaxFunnel then MaxFunnel := 200; elif Type(MaxFunnel) eq MonStgElt then MaxFunnel := StringToInteger(MaxFunnel); end if;
SetMemoryLimit(MemGB*10^9);
load "split_lab.m";  // run from product/code/

RQ := Rationals();

function E26(t)
    r1 := (-2*t+10)/((t+3)*(t-3));
    r2 := (-t^3+7*t^2-11*t+5)/(4*(t+3)*(t-3)^2);
    r3 := (-2*t^2+4*t-2)/((t+3)^2*(t-3));
    f := (RQx.1 - r1)*(RQx.1 - r2)*(RQx.1 - r3);
    return EllipticCurve(RQx!f);
end function;

// runtime verification of the family
ok26 := 0;
for tv in [RQ| 7, 9/2, -4/3 ] do
    try if Invariants(TorsionSubgroup(E26(tv))) eq [2,6] then ok26 +:= 1; end if; catch e; end try;
end for;
printf "VERIFY E26: %o/3\n", ok26;
error if ok26 lt 2, "universal curve reconstruction failed";

function SFrat(x)  // signed squarefree part of a nonzero rational
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

EXCL := {RQ|3,-3,1,5,9};
vals := [ v : v in HeightRats(H) | not v in EXCL ];
NV := #vals;
printf "SWEEP 2266 H=%o values=%o\n", H, NV;

// kernels of the three pair classes, per value
kk := [];
for v in vals do
    Append(~kk, [ SFrat((v+3)*(v-5)), SFrat(2*(v-3)), SFrat(-(v-1)*(v-9)) ]);
end for;

// ordered kernel: pairs indexed (1,2)->1, (1,3)->2, (2,3)->3; transposes negate
function KO(k, a, b)
    if a lt b then
        return a eq 1 select (b eq 2 select k[1] else k[2]) else k[3];
    else
        return -(b eq 1 select (a eq 2 select k[1] else k[2]) else k[3]);
    end if;
end function;

SIGMAS := [ [1,2,3],[2,1,3],[3,2,1],[1,3,2],[2,3,1],[3,1,2] ];

// hash join: single AA keyed by <si, k1, k2, k3> of the sigma-scrambled u-keys
M := AssociativeArray();
for ui in [1..NV] do
    for si in [1..6] do
        sg := SIGMAS[si];
        key := <si, KO(kk[ui],sg[1],sg[2]), KO(kk[ui],sg[1],sg[3]), KO(kk[ui],sg[2],sg[3])>;
        if IsDefined(M, key) then M[key] := M[key] cat [ui]; else M[key] := [ui]; end if;
    end for;
end for;

pairs := {};           // unordered index pairs as sorted 2-sequences
sigof := AssociativeArray();   // pair -> sigma indices that matched
for ti in [1..NV] do
    for si in [1..6] do
        key := <si, kk[ti][1], kk[ti][2], kk[ti][3]>;
        if not IsDefined(M, key) then continue; end if;
        for ui in M[key] do
            if ui eq ti then continue; end if;
            pr := ti lt ui select [ti,ui] else [ui,ti];
            Include(~pairs, pr);
            if IsDefined(sigof, pr) then sigof[pr] := sigof[pr] cat [si]; else sigof[pr] := [si]; end if;
        end for;
    end for;
end for;
plist := Sort(Setseq(pairs));
printf "MATCHES 2266 pairs=%o\n", #plist;

// ---- funnel the survivors ----
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
t0 := Cputime();
for pr in plist do
    tv := vals[pr[1]]; uv := vals[pr[2]];
    EA := 0; EB := 0; okc := true;
    try EA := E26(tv); EB := E26(uv); catch e okc := false; end try;
    if not okc then continue; end if;
    if jInvariant(EA) eq jInvariant(EB) then continue; end if;
    nf +:= 1;
    if nf gt MaxFunnel then printf "FUNNEL CAP %o reached\n", MaxFunnel; break; end if;
    printf "INSTANCE 2266 t=%o u=%o sigmas=%o\n", tv, uv, sigof[pr];
    GlueFunnel(EA, EB, Sprintf("l2266|t=%o|u=%o", tv, uv), [Integers()|3,3], ~seen, ~nstat);
end for;
printf "SEARCH_DONE 2266 H=%o pairs=%o funneled=%o curves=%o aborts=%o known=%o hits=%o fails=%o %o s\n",
    H, #plist, nf, #seen, nstat["abort"], nstat["known"], nstat["hit"], nstat["fail"], Cputime()-t0;
quit;
