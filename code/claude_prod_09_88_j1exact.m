// (8,8) lane, step 1.4: exact-J1 probe at one stage-1 base (m,n).
// Route 3 of notes/claude_prod_09_88.md: local (8,4) filters are provably
// useless on this locus, so we test exact TorsionSubgroup(J1) directly on
// double-stage-1 members -- v = 1/z for MW points (z,w) on
//   E_{s,t}: w^2 = z(z+s^2)(z+A),  A = s^2-t^4+t^2,
// the highest-probability spots for a (8,4)/(8,8) jackpot.
// Prints "J1 UPGRADE" loudly on anything beyond [4,4].
//
// Run (per base, via code/claude_prod_09_88_j1exact_run.sh):
//   magma -b Mnum:=3 Mden:=1 Nnum:=5 Nden:=1 code/claude_prod_09_88_j1exact.m

SetColumns(0);
SetSeed(1);

if not assigned Mnum then Mnum := 3; elif Type(Mnum) eq MonStgElt then Mnum := StringToInteger(Mnum); end if;
if not assigned Mden then Mden := 1; elif Type(Mden) eq MonStgElt then Mden := StringToInteger(Mden); end if;
if not assigned Nnum then Nnum := 5; elif Type(Nnum) eq MonStgElt then Nnum := StringToInteger(Nnum); end if;
if not assigned Nden then Nden := 1; elif Type(Nden) eq MonStgElt then Nden := StringToInteger(Nden); end if;
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned MaxK then MaxK := 2; elif Type(MaxK) eq MonStgElt then MaxK := StringToInteger(MaxK); end if;
if not assigned MaxV then MaxV := 40; elif Type(MaxV) eq MonStgElt then MaxV := StringToInteger(MaxV); end if;
SetMemoryLimit(MemGB*10^9);

load "code/claude_prod_09_88_defs.m";

m := Q_88!Mnum/Mden;
n := Q_88!Nnum/Nden;
s, t := StageOneST(m, n);
A := s^2 - t^4 + t^2;
printf "BASE m=%o n=%o s=%o t=%o\n", m, n, s, t;

// E_{s,t}: w^2 = z(z+s^2)(z+A)
E := EllipticCurve([0, s^2 + A, 0, s^2*A, 0]);
SetClassGroupBounds("GRH");   // descent shortcuts only; torsion results below are unconditional
Em, toEm := MinimalModel(E);
lb, ub := RankBounds(Em);
printf "RANKBOUNDS [%o,%o]\n", lb, ub;

G, mG := MordellWeilGroup(Em);
free := [mG(G.i) : i in [1..Ngens(G)] | Order(G.i) eq 0];
printf "MW rank_used=%o invariants=%o\n", #free, Invariants(G);
if #free eq 0 then
    printf "BASE_RANK0 -- no free part found, nothing to probe\n";
    printf "SEARCH_DONE m=%o n=%o tested=0 upgrades=0\n", m, n;
    quit;
end if;

Tors, mT := TorsionSubgroup(Em);
torpts := [mT(g) : g in Tors];

// v-candidates: z-coords of k.free + tor, k in a small box up to sign.
r := #free;
kvecs := [];
for tup in CartesianPower([-MaxK..MaxK], r) do
    k := [tup[i] : i in [1..r]];
    if &and[ki eq 0 : ki in k] then continue; end if;
    firstnz := [ki : ki in k | ki ne 0][1];
    if firstnz lt 0 then continue; end if;   // z is sign-invariant: quotient by -1
    Append(~kvecs, k);
end for;

vlist := [];
for k in kvecs do
    base := &+[k[i]*free[i] : i in [1..r]];
    for T in torpts do
        pt := base + T;
        ptE := pt @@ toEm;
        if ptE[3] eq 0 or ptE[2] eq 0 then continue; end if;   // infinity / 2-torsion
        z := ptE[1]/ptE[3];
        if z eq 0 then continue; end if;
        Append(~vlist, 1/z);
    end for;
end for;
vlist := Sort(Setseq(Seqset(vlist)));
if #vlist gt MaxV then
    // prefer small height: these minimize downstream model sizes
    hts := [Maximum(Abs(Numerator(v)), Abs(Denominator(v))) : v in vlist];
    ParallelSort(~hts, ~vlist);
    vlist := vlist[1..MaxV];
end if;
printf "V_CANDIDATES %o\n", #vlist;

tested := 0; upgrades := 0;
for v in vlist do
    ok := true;
    try
        h, g1 := Lambda334(s, t, v);
    catch e
        ok := false;
    end try;
    if not ok then continue; end if;
    inv1 := [0];
    try
        inv1 := ExactTorsion(g1);
    catch e
        printf "J1FAIL m=%o n=%o v=%o (torsion computation error)\n", m, n, v;
        continue;
    end try;
    tested +:= 1;
    printf "J1EXACT m=%o n=%o v=%o J1=%o\n", m, n, v, inv1;
    if inv1 ne [4,4] then
        upgrades +:= 1;
        printf "J1 UPGRADE m=%o n=%o v=%o J1=%o\n", m, n, v, inv1;
    end if;
end for;
printf "SEARCH_DONE m=%o n=%o tested=%o upgrades=%o\n", m, n, tested, upgrades;
quit;
