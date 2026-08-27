// lane_cong_sieve.m — trace-congruence sieve for (5,*) and (7,*) congruent
// pairs ACROSS the high-torsion families (far beyond the LMFDB conductor
// range of the Frengley/Cremona-Freitas lists).  A pair (E,F) with
// a_p(E) = a_p(F) mod L at >= MINCMP common good primes (no mismatch) is a
// candidate L-congruent pair; non-isogenous candidates with
// gcd(L, #T1*#T2) = 1 and T1 x T2 large feed the analytic gluer
// (analytic_glue.m), which finally accepts exactly the pairs carrying an
// anti-isometry.  Product torsion injects into the glued Jacobian:
//   [2,6]x[2,6] -> [2,2,6,6] (144)   [12]x[12] -> [12,12] (144)
//   [2,8]x[2,6] -> [2,2,2,24] (192)  [2,8]x[2,8] -> [2,2,8,8] (256)
//   [12]x[2,6] -> [2,6,12] (144)     [12]x[2,8] -> [2,4,24] (192)
//   [10]x[10] -> [10,10] (100)       [10]x[12] -> [2,60] (120)
//   [9]x[9] -> [9,9] (81)            [9]x[2,8] -> [2,72] (144)  etc.
// Usage: cd product/code && magma -b lane_cong_sieve.m > ../logs/lane_cong_sieve.log
//   optional: H:=<int> (default 40), Fam:="26" (E26-only pool; default "all")
SetColumns(0);
if not assigned H then H := 40; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned Fam then Fam := "all"; end if;   // "all" or "26"
SetMemoryLimit(6*10^9);
load "split_lab.m";   // KNOWN, RQx

RQ := Rationals();
QT<t> := FunctionField(RQ);

// ---------- families ----------
// tag, builder returning ok, E (with expected torsion checked cheaply)
function E26of(tv)
    if tv in {RQ|3,-3,1,5,9} then return false, 0; end if;
    r1 := (-2*tv+10)/((tv+3)*(tv-3));
    r2 := (-tv^3+7*tv^2-11*tv+5)/(4*(tv+3)*(tv-3)^2);
    r3 := (-2*tv^2+4*tv-2)/((tv+3)^2*(tv-3));
    if r1 eq r2 or r1 eq r3 or r2 eq r3 then return false, 0; end if;
    return true, EllipticCurve((RQx.1-r1)*(RQx.1-r2)*(RQx.1-r3));
end function;
function E28of(tv)
    if tv eq 0 or tv eq -1/4 then return false, 0; end if;
    if 8*tv^2-1 eq 0 then return false, 0; end if;
    r1 := (16*tv^3+12*tv^2+2*tv)/(8*tv^2-1)^2;
    r2 := (32*tv^3+24*tv^2+8*tv+1)/(16*tv^2*(8*tv^2-1));
    r3 := (-32*tv^4-32*tv^3-12*tv^2-2*tv)/((4*tv+1)^2*(8*tv^2-1));
    if r1 eq r2 or r1 eq r3 or r2 eq r3 then return false, 0; end if;
    return true, EllipticCurve((RQx.1-r1)*(RQx.1-r2)*(RQx.1-r3));
end function;
function TateOf(bv, cv)
    if bv eq 0 then return false, 0; end if;
    E := 0;
    try E := EllipticCurve([1-cv, -bv, -bv, 0, 0]); catch e return false, 0; end try;
    return true, E;
end function;
function X1N(N, tv)
    ok := true; bp := RQ!0; cp := RQ!0;
    try
        case N:
        when 7:  bp := tv^3-tv^2; cp := tv^2-tv;
        when 9:  bp := tv^2*(tv-1)*(tv^2-tv+1); cp := tv^2*(tv-1);
        when 10: bp := tv^3*(tv-1)*(2*tv-1)/(tv^2-3*tv+1)^2; cp := -tv*(tv-1)*(2*tv-1)/(tv^2-3*tv+1);
        when 12: bp := tv*(2*tv-1)*(3*tv^2-3*tv+1)*(2*tv^2-2*tv+1)/(tv-1)^4; cp := -tv*(2*tv-1)*(3*tv^2-3*tv+1)/(tv-1)^3;
        end case;
    catch e
        ok := false;
    end try;
    if not ok then return false, 0; end if;
    return TateOf(bp, cp);
end function;

function HeightRats(H)
    S := [];
    for a in [1..H], b in [1..H] do
        if GCD(a,b) eq 1 then Append(~S, a/b); Append(~S, -a/b); end if;
    end for;
    return S;
end function;

PRIMES := [ p : p in [11..230] | IsPrime(p) ];
NP := #PRIMES;

// pool records: <fam, tv, E, torsinvs, j, tracevec5, tracevec7>
pool := [* *];
jseen := {};
FAMS := [* <"26", 0>, <"28", 0>, <"x7", 7>, <"x9", 9>, <"x10", 10>, <"x12", 12> *];
if Fam eq "26" then FAMS := [* <"26", 0> *]; end if;
for fam in FAMS do
    cnt := 0;
    for tv in HeightRats(H) do
        ok := false; E := 0;
        if fam[1] eq "26" then ok, E := E26of(tv);
        elif fam[1] eq "28" then ok, E := E28of(tv);
        else ok, E := X1N(fam[2], tv);
        end if;
        if not ok then continue; end if;
        j := jInvariant(E);
        if j in jseen then continue; end if;
        Include(~jseen, j);
        Ti := Invariants(TorsionSubgroup(E));
        // guard: family degenerations
        n1 := IsEmpty(Ti) select 1 else &*Ti;
        if n1 lt 7 then continue; end if;
        Em := MinimalModel(E);
        NE := Conductor(Em);
        v5 := [ Integers() | ]; v7 := [ Integers() | ]; v3 := [ Integers() | ];
        for p in PRIMES do
            if NE mod p eq 0 then
                Append(~v5, 99); Append(~v7, 99); Append(~v3, 99);
            else
                ap := TraceOfFrobenius(Em, p);
                Append(~v5, ap mod 5); Append(~v7, ap mod 7); Append(~v3, ap mod 3);
            end if;
        end for;
        Append(~pool, < fam[1], tv, Em, Ti, j, v5, v7, v3 >);
        cnt +:= 1;
    end for;
    printf "POOL %o: %o curves\n", fam[1], cnt;
end for;
NPOOL := #pool;
printf "POOL total %o\n", NPOOL;

MINCMP := 14;

// hoisted parallel arrays for the hot loop
jsA  := [ pool[i][5] : i in [1..NPOOL] ];
v5A  := [ pool[i][6] : i in [1..NPOOL] ];
v7A  := [ pool[i][7] : i in [1..NPOOL] ];
v3A  := [ pool[i][8] : i in [1..NPOOL] ];
n1A  := [ IsEmpty(pool[i][4]) select 1 else &*pool[i][4] : i in [1..NPOOL] ];
ok5A := [ n1A[i] mod 5 ne 0 : i in [1..NPOOL] ];
ok7A := [ n1A[i] mod 7 ne 0 : i in [1..NPOOL] ];
ok3A := [ n1A[i] mod 3 ne 0 : i in [1..NPOOL] ];
MINCMP3 := 24;   // 3^-24 ~ 3e-12 per pair

procedure Report(i, k, L, ncmp, ~pool, ~ncand)
    A := pool[i]; B := pool[k];
    iso := IsIsogenous(A[3], B[3]);
    prod := Invariants(AbelianGroup(A[4] cat B[4]));
    nord := &*prod;
    isnew := not prod in KNOWN;
    printf "SIEVECAND L=%o %o(t=%o)%o x %o(t=%o)%o prod=%o order=%o cmp=%o %o%o\n",
        L, A[1], A[2], A[4], B[1], B[2], B[4], prod, nord, ncmp,
        isnew select "NEW" else "known", iso select " ISOGENOUS" else " NONISO";
    ncand +:= 1;
end procedure;

ncand := 0;
t0 := Cputime();
for i in [1..NPOOL] do
    ji := jsA[i]; v5i := v5A[i]; v7i := v7A[i]; v3i := v3A[i];
    o5i := ok5A[i]; o7i := ok7A[i]; o3i := ok3A[i];
    for k in [i+1..NPOOL] do
        if ji eq jsA[k] then continue; end if;
        if o5i and ok5A[k] then
            vb := v5A[k];
            ok := true; ncmp := 0;
            for m in [1..NP] do
                a := v5i[m]; b := vb[m];
                if a eq 99 or b eq 99 then continue; end if;
                if a ne b then ok := false; break; end if;
                ncmp +:= 1;
            end for;
            if ok and ncmp ge MINCMP then Report(i, k, 5, ncmp, ~pool, ~ncand); end if;
        end if;
        if o7i and ok7A[k] then
            vb := v7A[k];
            ok := true; ncmp := 0;
            for m in [1..NP] do
                a := v7i[m]; b := vb[m];
                if a eq 99 or b eq 99 then continue; end if;
                if a ne b then ok := false; break; end if;
                ncmp +:= 1;
            end for;
            if ok and ncmp ge MINCMP then Report(i, k, 7, ncmp, ~pool, ~ncand); end if;
        end if;
        if o3i and ok3A[k] then
            vb := v3A[k];
            ok := true; ncmp := 0;
            for m in [1..NP] do
                a := v3i[m]; b := vb[m];
                if a eq 99 or b eq 99 then continue; end if;
                if a ne b then ok := false; break; end if;
                ncmp +:= 1;
            end for;
            if ok and ncmp ge MINCMP3 then Report(i, k, 3, ncmp, ~pool, ~ncand); end if;
        end if;
    end for;
    if i mod 500 eq 0 then
        printf "PROGRESS sieve %o/%o rows, %o cands, %o s\n", i, NPOOL, ncand, Cputime()-t0;
    end if;
end for;
printf "SIEVE_DONE cands=%o %o s\n", ncand, Cputime()-t0;
quit;
