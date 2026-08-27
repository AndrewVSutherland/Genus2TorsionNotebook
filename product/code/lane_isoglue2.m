// lane_isoglue2.m — STAGE 1 of the isogeny-class harvest, restructured after
// the Mestre bottleneck (HyperellipticCurveFromIgusaClebsch stalls on the
// ~90-digit glued invariants; curve construction now runs as separate
// dedicated jobs).  This lane only SCANS: for the best candidate pairs
// (product order desc, N=5 preferred = 120 Lagrangian M's, conductor asc)
// it runs GlueScan at high precision and prints every Galois-equivariant
// anti-isometry with its rational Igusa-Clebsch invariants; results feed
// the per-hit Mestre queue.
// Usage: cd product/code && magma -b lane_isoglue2.m > ../logs/lane_isoglue2.log
//   optional: H:=<int> (default 25), MaxScan:=<int> (default 6), Prec:=<int> (default 220)
SetColumns(0);
if not assigned H then H := 25; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned MaxScan then MaxScan := 6; elif Type(MaxScan) eq MonStgElt then MaxScan := StringToInteger(MaxScan); end if;
if not assigned Prec then Prec := 220; elif Type(Prec) eq MonStgElt then Prec := StringToInteger(Prec); end if;
SetMemoryLimit(8*10^9);
load "split_lab.m";
load "analytic_glue.m";

RQ := Rationals();

function E26of(tv)
    if tv in {RQ|3,-3,1,5,9} then return false, 0; end if;
    r1 := (-2*tv+10)/((tv+3)*(tv-3));
    r2 := (-tv^3+7*tv^2-11*tv+5)/(4*(tv+3)*(tv-3)^2);
    r3 := (-2*tv^2+4*tv-2)/((tv+3)^2*(tv-3));
    if r1 eq r2 or r1 eq r3 or r2 eq r3 then return false, 0; end if;
    return true, EllipticCurve((RQx.1-r1)*(RQx.1-r2)*(RQx.1-r3));
end function;
function E28of(tv)
    if tv eq 0 or tv eq -1/4 or 8*tv^2-1 eq 0 then return false, 0; end if;
    r1 := (16*tv^3+12*tv^2+2*tv)/(8*tv^2-1)^2;
    r2 := (32*tv^3+24*tv^2+8*tv+1)/(16*tv^2*(8*tv^2-1));
    r3 := (-32*tv^4-32*tv^3-12*tv^2-2*tv)/((4*tv+1)^2*(8*tv^2-1));
    if r1 eq r2 or r1 eq r3 or r2 eq r3 then return false, 0; end if;
    return true, EllipticCurve((RQx.1-r1)*(RQx.1-r2)*(RQx.1-r3));
end function;
function X1N(N, tv)
    ok := true; bp := RQ!0; cp := RQ!0;
    try
        case N:
        when 10: bp := tv^3*(tv-1)*(2*tv-1)/(tv^2-3*tv+1)^2; cp := -tv*(tv-1)*(2*tv-1)/(tv^2-3*tv+1);
        when 12: bp := tv*(2*tv-1)*(3*tv^2-3*tv+1)*(2*tv^2-2*tv+1)/(tv-1)^4; cp := -tv*(2*tv-1)*(3*tv^2-3*tv+1)/(tv-1)^3;
        end case;
    catch e
        ok := false;
    end try;
    if not ok or bp eq 0 then return false, 0; end if;
    E := 0;
    try E := EllipticCurve([1-cp, -bp, -bp, 0, 0]); catch e return false, 0; end try;
    return true, E;
end function;
function HeightRats(H)
    S := [];
    for a in [1..H], b in [1..H] do
        if GCD(a,b) eq 1 then Append(~S, a/b); Append(~S, -a/b); end if;
    end for;
    return S;
end function;

pool := [* *];
jseen := {};
for fam in [* <"26",0>, <"28",0>, <"x10",10>, <"x12",12> *] do
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
        n1 := IsEmpty(Ti) select 1 else &*Ti;
        if n1 lt 10 then continue; end if;
        Append(~pool, < fam[1], tv, MinimalModel(E), Ti, j, n1 >);
    end for;
end for;
printf "POOL total %o\n", #pool;

cands := [* *];
seencls := {};
for i in [1..#pool] do
    E := pool[i][3];
    cls := [* *];
    ok := true;
    try cls := IsogenousCurves(E); catch e ok := false; end try;
    if not ok or #cls eq 0 then continue; end if;
    clskey := Sprintf("%o", Sort([ jInvariant(cls[k]) : k in [1..#cls] ]));
    if clskey in seencls then continue; end if;
    Include(~seencls, clskey);
    Ts := [ Invariants(TorsionSubgroup(cls[k])) : k in [1..#cls] ];
    for k1 in [1..#cls] do
        for k2 in [k1+1..#cls] do
            T1 := Ts[k1]; T2 := Ts[k2];
            n12 := (IsEmpty(T1) select 1 else &*T1) * (IsEmpty(T2) select 1 else &*T2);
            if n12 lt 100 then continue; end if;
            prod := Invariants(AbelianGroup(T1 cat T2));
            if prod in KNOWN then continue; end if;
            for N in [5,7] do
                if GCD(N, n12) ne 1 then continue; end if;
                Append(~cands, < cls[k1], cls[k2], N, 0, T1, T2, prod, Conductor(cls[k1]) >);
            end for;
        end for;
    end for;
end for;
printf "CANDS %o\n", #cands;

// sort: product order desc, then N=5 before N=7, then conductor asc
keyed := [ < &*cands[i][7], cands[i][3] eq 5 select 0 else 1, cands[i][8], i > : i in [1..#cands] ];
Sort(~keyed, func<a,b | a[1] ne b[1] select b[1]-a[1] else (a[2] ne b[2] select a[2]-b[2] else a[3]-b[3]) >);
nscan := 0;
for k in keyed do
    if nscan ge MaxScan then break; end if;
    c := cands[k[4]];
    nscan +:= 1;
    printf "==== SCAN %o: N=%o %o x %o -> %o (cond %o) ====\n", nscan, c[3], c[5], c[6], c[7], c[8];
    hits := GlueScan(c[1], c[2], c[3] : prec := Prec);
    for hh in hits do
        printf "SCANHIT N=%o cond=%o M=%o aE=%o aF=%o\n", c[3], c[8], hh[1],
            aInvariants(c[1]), aInvariants(c[2]);
        printf "SCANINV %o %o %o\n", hh[2], hh[3], hh[4];
    end for;
    printf "SCANSUM %o hits for pair %o\n", #hits, nscan;
end for;
printf "ISOGLUE2_DONE scans=%o\n", nscan;
quit;
