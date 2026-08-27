// lane_2266_fiber2.m — [2,2,6,6] structured cross attack, section-arithmetic
// version (replaces the MordellWeilGroup-per-fiber lane_2266_fiber.m, which
// had unbounded per-fiber cost).  On the fiber over u of HLP's (6)-surface
//   Z^2 = (u+3)(u-5)((u-3)y^2+6)((u-3)y^2-2),  Z = (u+3)(u-5) z,
// we use ONLY the known rational points (y,Z) = (+-1, +-(u+3)(u-5)) and small
// integer combinations of their images on the Jacobian (no descent), giving
// t = 3+(u-3)y^2 on the C1&C2-locus; then test C3 and funnel survivors.
// Usage: cd product/code && magma -b lane_2266_fiber2.m > ../logs/lane2266_fiber2.log
//   optional: HU:=<int> (default 25), NB:=<int> (default 6)
SetColumns(0);
if not assigned HU then HU := 25; elif Type(HU) eq MonStgElt then HU := StringToInteger(HU); end if;
if not assigned NB then NB := 6; elif Type(NB) eq MonStgElt then NB := StringToInteger(NB); end if;
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
load "split_lab.m";  // run from product/code/

RQ := Rationals();
QY<Y> := PolynomialRing(RQ);

function E26(t)
    r1 := (-2*t+10)/((t+3)*(t-3));
    r2 := (-t^3+7*t^2-11*t+5)/(4*(t+3)*(t-3)^2);
    r3 := (-2*t^2+4*t-2)/((t+3)^2*(t-3));
    f := (RQx.1 - r1)*(RQx.1 - r2)*(RQx.1 - r3);
    return EllipticCurve(RQx!f);
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

function HeightRats(H)
    S := [];
    for a in [1..H], b in [1..H] do
        if GCD(a,b) eq 1 then Append(~S, a/b); Append(~S, -a/b); end if;
    end for;
    return S;
end function;

EXCL := {RQ|3,-3,1,5,9};
HeightOf := func<q | Max(Abs(Numerator(q)), Abs(Denominator(q)))>;

survivors := {@ @};
nfib := 0; npts := 0;
t0 := Cputime();
for u in HeightRats(HU) do
    if u in EXCL then continue; end if;
    D0 := (u+3)*(u-5);
    quart := D0*((u-3)*Y^2+6)*((u-3)*Y^2-2);
    if Degree(quart) ne 4 or Discriminant(quart) eq 0 then continue; end if;
    nfib +:= 1;
    C := HyperellipticCurve(quart);
    ok := true; E := 0; mE := 0;
    try E, mE := EllipticCurve(C, C![1, D0]); catch e ok := false; end try;
    if not ok then continue; end if;
    P1 := 0; P2 := 0;
    try P1 := mE(C![-1, D0]); P2 := mE(C![1, -D0]); catch e ok := false; end try;
    if not ok then continue; end if;
    mEi := 0;
    try mEi := Inverse(mE); catch e ok := false; end try;
    if not ok then continue; end if;
    for a in [-NB..NB] do
        for b in [-2..2] do
            P := a*P1 + b*P2;
            if P eq E!0 then continue; end if;
            Q := 0; okq := true;
            try Q := mEi(P); catch e okq := false; end try;
            if not okq or Q[3] eq 0 then continue; end if;
            yv := Q[1]/Q[3];
            tv := 3 + (u-3)*yv^2;
            if tv in EXCL or tv eq u then continue; end if;
            if HeightOf(tv) gt 10^12 then continue; end if;
            npts +:= 1;
            v3 := (tv-1)*(tv-9)*(u-1)*(u-9);
            if v3 eq 0 or not IsSquare(v3) then continue; end if;
            // C1, C2 hold by construction (t = 3+(u-3)y^2 on the surface)
            pr := <tv, u>;
            if pr in survivors then continue; end if;
            Include(~survivors, pr);
            printf "FIBSURV t=%o u=%o y=%o (a,b)=(%o,%o)\n", tv, u, yv, a, b;
        end for;
    end for;
    if nfib mod 200 eq 0 then
        printf "PROGRESS fiberA %o fibers %o pts %o surv %o s\n", nfib, npts, #survivors, Cputime()-t0;
        System(Sprintf("echo 'fiberA %o fibers %o pts %o surv' >> ../logs/fiberA.progress", nfib, npts, #survivors));
    end if;
end for;
printf "FIBER2_SCAN fibers=%o pts=%o survivors=%o %o s\n", nfib, npts, #survivors, Cputime()-t0;

for pr in survivors do
    tv := pr[1]; uv := pr[2];
    EA := 0; EB := 0; okc := true;
    try EA := E26(tv); EB := E26(uv); catch e okc := false; end try;
    if not okc then continue; end if;
    if jInvariant(EA) eq jInvariant(EB) then printf "SKIPISO t=%o u=%o\n", tv, uv; continue; end if;
    printf "INSTANCE 2266fib2 t=%o u=%o\n", tv, uv;
    GlueFunnel(EA, EB, Sprintf("l2266fib2|t=%o|u=%o", tv, uv), [Integers()|3,3], ~seen, ~nstat);
end for;
printf "SEARCH_DONE 2266fiber2 survivors=%o curves=%o aborts=%o known=%o hits=%o fails=%o %o s\n",
    #survivors, #seen, nstat["abort"], nstat["known"], nstat["hit"], nstat["fail"], Cputime()-t0;
quit;
