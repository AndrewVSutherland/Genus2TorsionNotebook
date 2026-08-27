// lane_2266_sigma.m — [2,2,6,6] via the five non-identity sigma-surfaces
// (2026-08-14 plan, section 1, step 1).  For each sigma in S3 \ {id} and each
// drop-choice d in {1,2,3}: 2-of-3 hash-join on the pair-class kernels
//   A12 = (t+3)(t-5),  A13 = 2(t-3),  A23 = -(t-1)(t-9)
// (ordered; transposes negate; condition j of the sigma-system is
//   KO(t; p[1], p[2]) = KO(u; sigma(p[1]), sigma(p[2])),  p = PAIRS[j]).
// Every surface point is written to ../data/sigma_si<si>_d<d>_pts.txt
// (columns: t u deck; deck=1 iff j-equal deck pair, -1 if curve build failed);
// the dropped third condition is tested on the spot (kernel equality on the
// dropped coordinate) and full passes are funneled: a genuine survivor prints
// HIT ... invs=[ 2, 2, 6, 6 ].
// sigma index map: si2=[2,1,3]=(12)  si3=[3,2,1]=(13)  si4=[1,3,2]=(23)
//                  si5=[2,3,1]=(123) si6=[3,1,2]=(132)
// Modes:
//   Census:=1   count-only, all six sigmas, lane_misc2-compatible counts
//               (ordered pairs incl. t=u, multiplicity-summed), no files/funnel
//   Scramble:=1 negative control: swap the two kept coordinates of the t-side
//               probe key (matches now require a codim-2 accident -> ~0)
// Usage: cd product/code && magma -b H:=150 lane_2266_sigma.m > ../logs/lane2266_sigma.log
//   optional: H:=<int> (default 150), MaxFunnel:=<int> (default 200),
//             Census:=<0|1>, Scramble:=<0|1>, MemGB:=<int>
SetColumns(0);
if not assigned H then H := 150; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned MaxFunnel then MaxFunnel := 200; elif Type(MaxFunnel) eq MonStgElt then MaxFunnel := StringToInteger(MaxFunnel); end if;
if not assigned Census then Census := 0; elif Type(Census) eq MonStgElt then Census := StringToInteger(Census); end if;
if not assigned Scramble then Scramble := 0; elif Type(Scramble) eq MonStgElt then Scramble := StringToInteger(Scramble); end if;
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
printf "SWEEP 2266sigma H=%o values=%o census=%o scramble=%o\n", H, NV, Census, Scramble;

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
PAIRS  := [ [1,2], [1,3], [2,3] ];

// lazy j-invariant cache: 0 = not computed, RQ value, or "bad"
jval := AssociativeArray();
function GetJ(idx, vals, jv)
    if IsDefined(jv, idx) then return jv[idx], jv; end if;
    r := "bad";
    try r := jInvariant(E26(vals[idx])); catch e r := "bad"; end try;
    jv[idx] := r;
    return r, jv;
end function;

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

silist := Census eq 1 select [1..6] else [2..6];
funneled := {};   // unordered value-index pairs already funneled (any si/drop)
nf := 0;
t0 := Cputime();

for si in silist do
    sg := SIGMAS[si];
    for drop in [1..3] do
        keep := [ j : j in [1..3] | j ne drop ];
        // u-side hash on the kept sigma-scrambled coordinates
        M := AssociativeArray();
        for ui in [1..NV] do
            key := < KO(kk[ui], sg[PAIRS[j][1]], sg[PAIRS[j][2]]) : j in keep >;
            if IsDefined(M, key) then Append(~M[key], ui); else M[key] := [ui]; end if;
        end for;
        nraw := 0;       // lane_misc2-compatible: ordered pairs incl. t=u
        lines := [];
        nfull := 0; ndeck := 0; ngenuine := 0;
        for ti in [1..NV] do
            if Scramble eq 1 then
                key := < kk[ti][keep[2]], kk[ti][keep[1]] >;
            else
                key := < kk[ti][keep[1]], kk[ti][keep[2]] >;
            end if;
            if not IsDefined(M, key) then continue; end if;
            nraw +:= #M[key];
            if Census eq 1 then continue; end if;
            for ui in M[key] do
                if ui eq ti then continue; end if;
                tv := vals[ti]; uv := vals[ui];
                jA, jval := GetJ(ti, vals, jval);
                jB, jval := GetJ(ui, vals, jval);
                deck := -1;
                if Type(jA) ne MonStgElt and Type(jB) ne MonStgElt then
                    deck := jA eq jB select 1 else 0;
                end if;
                Append(~lines, Sprintf("%o %o %o", tv, uv, deck));
                // dropped third condition: kernel equality on the dropped coord
                pd := PAIRS[drop];
                if KO(kk[ti], pd[1], pd[2]) eq KO(kk[ui], sg[pd[1]], sg[pd[2]]) then
                    nfull +:= 1;
                    if deck eq 1 then ndeck +:= 1; continue; end if;
                    if deck eq -1 then continue; end if;
                    ngenuine +:= 1;
                    pr := ti lt ui select [ti,ui] else [ui,ti];
                    if pr in funneled then continue; end if;
                    Include(~funneled, pr);
                    nf +:= 1;
                    if nf gt MaxFunnel then printf "FUNNEL CAP %o reached\n", MaxFunnel; continue; end if;
                    printf "INSTANCE 2266sigma t=%o u=%o si=%o drop=%o\n", tv, uv, si, drop;
                    EA := 0; EB := 0; okc := true;
                    try EA := E26(tv); EB := E26(uv); catch e okc := false; end try;
                    if not okc then continue; end if;
                    GlueFunnel(EA, EB, Sprintf("l2266s|si=%o|t=%o|u=%o", si, tv, uv), [Integers()|3,3], ~seen, ~nstat);
                end if;
            end for;
        end for;
        printf "sigma=%o drop=%o: %o (t,u) surface points at H=%o\n", sg, drop, nraw, H;
        if Census eq 0 then
            printf "SURFACE si=%o drop=%o rawmatch=%o listed=%o full=%o deckfull=%o genuinefull=%o\n",
                si, drop, nraw, #lines, nfull, ndeck, ngenuine;
            fn := Sprintf("../data/sigma_si%o_d%o_pts.txt", si, drop);
            System(Sprintf("rm -f %o", fn));
            if #lines gt 0 then
                PrintFile(fn, &cat[ l cat "\n" : l in lines ] : Overwrite := true);
            end if;
        end if;
    end for;
end for;

printf "SEARCH_DONE 2266sigma H=%o funneled=%o curves=%o aborts=%o known=%o hits=%o fails=%o %o s\n",
    H, nf, #seen, nstat["abort"], nstat["known"], nstat["hit"], nstat["fail"], Cputime()-t0;
quit;
