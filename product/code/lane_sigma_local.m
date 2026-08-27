// lane_sigma_local.m — Stage 5a: local solubility (Q_p for p <= PMAX, and R)
// of the FULL 3-condition sigma-surface variety, for each non-identity sigma
// (2026-08-14 plan, section 3.1).  The variety for sigma is the set of (t,u)
// with all three products
//   V_j(t,u) = A(pd[1],pd[2])(t) * A(sg[pd[1]],sg[pd[2]])(u),  pd = PAIRS[j],
// NONZERO squares, where the ordered pair-class VALUES are
//   A(1,2)(x) = (x+3)(x-5),  A(1,3)(x) = 2(x-3),  A(2,3)(x) = -(x-1)(x-9),
// and transposes negate (A(b,a) = -A(a,b)).  A rational pair (t,u) whose
// three products are squares in Q_p certifies Q_p-solubility; the test is
// exact rational arithmetic:
//   p odd : v_p(x) even and unit part a QR mod p (Legendre(a*b, p) = 1)
//   p = 2 : v_2(x) even and unit part = 1 mod 8 ((a*b) mod 8 = 1, a,b odd)
//   R     : x > 0                          (place label p = 0)
// EXPECTATION: the H<=200 full-system joins matched thousands of j-equal DECK
// pairs, which are rational (hence everywhere-local) points on these
// varieties; SOLUBLE everywhere is therefore the expected outcome (ELS), and
// any NOWITNESS/OBSTRUCTED? line is the interesting signal.
// Markers: LOCAL si=.. p=.. SOLUBLE t=.. u=.. | NOWITNESS_H.. | OBSTRUCTED?
//          LOCALSUM si=.. ; LOCAL_DONE
// Usage: cd product/code && magma -b lane_sigma_local.m > ../logs/lane_sigma_local.log
//   optional: HB:=<int> (default 24), PMAX:=<int> (default 50), MemGB:=<int>
SetColumns(0);
if not assigned HB then HB := 24; elif Type(HB) eq MonStgElt then HB := StringToInteger(HB); end if;
if not assigned PMAX then PMAX := 50; elif Type(PMAX) eq MonStgElt then PMAX := StringToInteger(PMAX); end if;
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
SetSeed(1);

RQ := Rationals();
SIGMAS := [ [1,2,3],[2,1,3],[3,2,1],[1,3,2],[2,3,1],[3,1,2] ];
PAIRS  := [ [1,2], [1,3], [2,3] ];
EXCL := {RQ|3,-3,1,5,9};   // all roots of the three class polynomials + poles
HeightOf := func<q | Max(Abs(Numerator(q)), Abs(Denominator(q)))>;

function HeightRats(H)
    S := [];
    for a in [1..H], b in [1..H] do
        if GCD(a,b) eq 1 then Append(~S, a/b); Append(~S, -a/b); end if;
    end for;
    return S;
end function;

function ClassTriple(x)   // [A(1,2), A(1,3), A(2,3)] at x
    return [ (x+3)*(x-5), 2*(x-3), -(x-1)*(x-9) ];
end function;

function AO(c, a, b)      // ordered class VALUE from a precomputed triple
    if a lt b then
        return a eq 1 select (b eq 2 select c[1] else c[2]) else c[3];
    else
        return -(b eq 1 select (a eq 2 select c[1] else c[2]) else c[3]);
    end if;
end function;

function IsSqQp(x, p)     // x nonzero rational; p = 0 means R
    if p eq 0 then return x gt 0; end if;
    v := Valuation(x, p);
    if IsOdd(v) then return false; end if;
    u := x / p^v;
    a := Numerator(u); b := Denominator(u);
    if p eq 2 then return ((a*b) mod 8) eq 1; end if;
    return LegendreSymbol(a*b, p) eq 1;
end function;

// test one (t,u) via precomputed triples ct, cu
function AllSq(ct, cu, sg, p)
    for j in [1..3] do
        pd := PAIRS[j];
        if not IsSqQp(AO(ct, pd[1], pd[2]) * AO(cu, sg[pd[1]], sg[pd[2]]), p) then
            return false;
        end if;
    end for;
    return true;
end function;

vals := [ v : v in HeightRats(HB) | not v in EXCL ];
Sort(~vals, func<a, b | HeightOf(a) - HeightOf(b)>);   // small heights first
NV := #vals;
CT := [ ClassTriple(v) : v in vals ];
printf "LOCALSCAN HB=%o values=%o PMAX=%o\n", HB, NV, PMAX;

places := [0] cat PrimesUpTo(PMAX);
t0 := Cputime();
for si in [2..6] do
    sg := SIGMAS[si];
    nowit := [];
    for p in places do
        found := false;
        // diagonal sweep: max index k grows, so small-height pairs come first
        for k in [1..NV] do
            for ti in [1..k] do
                cands := ti eq k select [ [k,k] ] else [ [ti,k], [k,ti] ];
                for pr in cands do
                    if AllSq(CT[pr[1]], CT[pr[2]], sg, p) then
                        printf "LOCAL si=%o p=%o SOLUBLE t=%o u=%o\n", si, p, vals[pr[1]], vals[pr[2]];
                        found := true; break;
                    end if;
                end for;
                if found then break; end if;
            end for;
            if found then break; end if;
        end for;
        if found then continue; end if;
        printf "LOCAL si=%o p=%o NOWITNESS_H%o\n", si, p, HB;
        // finer pass 1: random pairs of height <= 60
        for trial in [1..400] do
            tv := Random(-60,60)/Random(1,60);
            uv := Random(-60,60)/Random(1,60);
            if tv in EXCL or uv in EXCL then continue; end if;
            if AllSq(ClassTriple(tv), ClassTriple(uv), sg, p) then
                printf "LOCAL si=%o p=%o SOLUBLE t=%o u=%o (random pass)\n", si, p, tv, uv;
                found := true; break;
            end if;
        end for;
        if found then continue; end if;
        // finer pass 2 (p > 0 only): integer residue scan t,u in [0..p^2-1]
        if p gt 0 then
            ntest := 0; partial := false;
            CU := [ <RQ!c0, ClassTriple(RQ!c0)> : c0 in [0..p^2-1] | not (RQ!c0) in EXCL ];
            for tt in CU do
                for uu in CU do
                    ntest +:= 1;
                    if ntest gt 2000000 then partial := true; break; end if;
                    if AllSq(tt[2], uu[2], sg, p) then
                        printf "LOCAL si=%o p=%o SOLUBLE t=%o u=%o (residue scan)\n", si, p, tt[1], uu[1];
                        found := true; break;
                    end if;
                end for;
                if found or partial then break; end if;
            end for;
            if found then continue; end if;
            printf "LOCAL si=%o p=%o OBSTRUCTED?%o\n", si, p, partial select " (partial scan)" else "";
        else
            printf "LOCAL si=%o p=%o OBSTRUCTED? (R)\n", si, p;
        end if;
        Append(~nowit, p);
    end for;
    printf "LOCALSUM si=%o sigma=%o places=%o soluble=%o nowitness=%o\n",
        si, sg, #places, #places - #nowit, nowit;
    System(Sprintf("echo 'sigmalocal si=%o done nowitness=%o' >> ../logs/sigmalocal.progress", si, #nowit));
end for;
printf "LOCAL_DONE %o s\n", Cputime()-t0;
quit;
