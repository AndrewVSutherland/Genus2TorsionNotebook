// lane_hlp37.m — instantiate the HLP §3.7 "gaining 2-power torsion" families:
//   [2,2,4,8]: glue E28(t) x E28(u) with psi = identity, condition (5):
//              (8t^2-1)(8t^2+8t+1) == (8u^2-1)(8u^2+8u+1) mod squares;
//   [2,2,24]:  glue E26(t) x E28(u), psi = (1 3 2), conditions:
//              2(t-3) = -y^2  and (t+3)(t-5) == (8u^2-1)(8u^2+8u+1) mod squares.
// Universal curves reconstructed from HLP Table 5 (x-coords of 2-torsion):
//   E26(t): roots (-2t+10)/((t+3)(t-3)), (-t^3+7t^2-11t+5)/(4(t+3)(t-3)^2),
//           (-2t^2+4t-2)/((t+3)^2(t-3))
//   E28(u): roots (16u^3+12u^2+2u)/(8u^2-1)^2, (32u^3+24u^2+8u+1)/(16u^2(8u^2-1)),
//           (-32u^4-32u^3-12u^2-2u)/((4u+1)^2(8u^2-1))
// (runtime-verified to have torsion [2,6] resp. [2,8]).
// Usage: magma -b lane_hlp37.m > ../logs/hlp37.log
SetColumns(0);
SetMemoryLimit(4*10^9);
load "split_lab.m";  // run from product/code/

function E26(t)
    r1 := (-2*t+10)/((t+3)*(t-3));
    r2 := (-t^3+7*t^2-11*t+5)/(4*(t+3)*(t-3)^2);
    r3 := (-2*t^2+4*t-2)/((t+3)^2*(t-3));
    f := (RQx.1 - r1)*(RQx.1 - r2)*(RQx.1 - r3);
    return EllipticCurve(RQx!f);
end function;

function E28(u)
    r1 := (16*u^3+12*u^2+2*u)/(8*u^2-1)^2;
    r2 := (32*u^3+24*u^2+8*u+1)/(16*u^2*(8*u^2-1));
    r3 := (-32*u^4-32*u^3-12*u^2-2*u)/((4*u+1)^2*(8*u^2-1));
    f := (RQx.1 - r1)*(RQx.1 - r2)*(RQx.1 - r3);
    return EllipticCurve(RQx!f);
end function;

// runtime verification
ok26 := 0; ok28 := 0;
for tv in [Rationals()| 7, 9/2, -4/3 ] do
    try if Invariants(TorsionSubgroup(E26(tv))) eq [2,6] then ok26 +:= 1; end if; catch e; end try;
    try if Invariants(TorsionSubgroup(E28(tv))) eq [2,8] then ok28 +:= 1; end if; catch e; end try;
end for;
printf "VERIFY E26: %o/3, E28: %o/3\n", ok26, ok28;
error if ok26 lt 2 or ok28 lt 2, "universal curve reconstruction failed";

function HeightRats(H)
    S := [];
    for a in [1..H], b in [1..H] do
        if GCD(a,b) eq 1 then Append(~S, a/b); Append(~S, -a/b); end if;
    end for;
    return S;
end function;

D8 := func<u | (8*u^2-1)*(8*u^2+8*u+1)>;

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

// ---- [2,2,24]: HLP's own point first, then sweep (y,u) ----
printf "== [2,2,24] hunt ==\n";
ys := HeightRats(16); us := HeightRats(16);
n2224 := 0;
for y in [Rationals()| 2/9 ] cat ys do
    tv := 3 - y^2/2;
    lhs := (tv+3)*(tv-5);
    if lhs eq 0 then continue; end if;
    for u in [Rationals()| 1/3 ] cat us do
        if u eq 0 or (8*u^2-1) eq 0 then continue; end if;
        rhs := D8(u);
        if rhs eq 0 then continue; end if;
        if not IsSquare(lhs/rhs) then continue; end if;
        n2224 +:= 1;
        if n2224 gt 40 then break y; end if;
        EA := 0; EB := 0; okc := true;
        try EA := E26(tv); EB := E28(u); catch e okc := false; end try;
        if not okc then continue; end if;
        if jInvariant(EA) eq jInvariant(EB) then continue; end if;
        printf "INSTANCE 2224 y=%o t=%o u=%o\n", y, tv, u;
        GlueFunnel(EA, EB, Sprintf("hlp2224|t=%o|u=%o", tv, u), [Integers()|3], ~seen, ~nstat);
    end for;
end for;

// ---- [2,2,4,8]: sweep (t,u) with (5) ----
printf "== [2,2,4,8] hunt ==\n";
ts := HeightRats(20);
n2248 := 0;
for i in [1..#ts] do
    tv := ts[i];
    if tv eq 0 or (8*tv^2-1) eq 0 then continue; end if;
    dt := D8(tv);
    if dt eq 0 then continue; end if;
    for j in [i+1..#ts] do
        u := ts[j];
        if u eq 0 or (8*u^2-1) eq 0 then continue; end if;
        du := D8(u);
        if du eq 0 then continue; end if;
        if not IsSquare(dt/du) then continue; end if;
        n2248 +:= 1;
        if n2248 gt 40 then break i; end if;
        EA := 0; EB := 0; okc := true;
        try EA := E28(tv); EB := E28(u); catch e okc := false; end try;
        if not okc then continue; end if;
        if jInvariant(EA) eq jInvariant(EB) then continue; end if;
        printf "INSTANCE 2248 t=%o u=%o\n", tv, u;
        GlueFunnel(EA, EB, Sprintf("hlp2248|t=%o|u=%o", tv, u), [Integers()|], ~seen, ~nstat);
    end for;
end for;

printf "SEARCH_DONE hlp37 inst2224 %o inst2248 %o aborts %o known %o hits %o fails %o\n",
    n2224, n2248, nstat["abort"], nstat["known"], nstat["hit"], nstat["fail"];
quit;
