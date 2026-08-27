// lane_2266_fiber.m — [2,2,6,6] structured cross attack (plan section 1,
// step 4).  The sigma=id two-gain system on E26(t) x E26(u) is
//   C1: (t-3)(u-3) sq,  C2: (t+3)(t-5)(u+3)(u-5) sq,  C3: (t-1)(t-9)(u-1)(u-9) sq.
// HLP's (6)-surface parametrizes C1 & C2:  t = 3 + (u-3)y^2 with
//   Z^2 = (u+3)(u-5)*((u-3)y^2+6)*((u-3)y^2-2),  Z = (u+3)(u-5)z,
// a genus-1 fibration over the u-line with rational points (y,Z) =
// (+-1, +-(u+3)(u-5)).  For each small u we enumerate the fiber's Mordell-Weil
// points and test C3; a survivor (t,u) gives full rational J[2] plus odd part
// [3,3] on the glued Jacobian: J(Q)_tors >= [2,2,6,6] (order 144).
// Usage: cd product/code && magma -b lane_2266_fiber.m > ../logs/lane2266_fiber.log
//   optional: HU:=<int> (default 12) height bound for u; NB:=<int> (default 10)
SetColumns(0);
if not assigned HU then HU := 12; elif Type(HU) eq MonStgElt then HU := StringToInteger(HU); end if;
if not assigned NB then NB := 10; elif Type(NB) eq MonStgElt then NB := StringToInteger(NB); end if;
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
SetClassGroupBounds("GRH");   // speeds the per-fiber descents
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

survivors := {@ @};   // <t,u>
nfib := 0; nfibok := 0; npts := 0;
t0 := Cputime();
for u in HeightRats(HU) do
    if u in EXCL then continue; end if;
    quart := (u+3)*(u-5)*((u-3)*Y^2+6)*((u-3)*Y^2-2);
    if Degree(quart) ne 4 or Discriminant(quart) eq 0 then continue; end if;
    nfib +:= 1;
    C := HyperellipticCurve(quart);
    pt := C![1, (u+3)*(u-5)];
    E := 0; mE := 0; ok := true;
    try E, mE := EllipticCurve(C, pt); catch e ok := false; end try;
    if not ok then continue; end if;
    G := 0; mw := 0;
    try G, mw := MordellWeilGroup(E); catch e ok := false; end try;
    if not ok then printf "MWFAIL u=%o\n", u; continue; end if;
    nfibok +:= 1;
    mEi := 0;
    try mEi := Inverse(mE); catch e ok := false; end try;
    if not ok then continue; end if;
    ng := Ngens(G);
    rk := #[ i : i in [1..ng] | Order(G.i) eq 0 ];
    ranges := [ Order(G.i) eq 0 select [-NB..NB] else [0..Order(G.i)-1] : i in [1..ng] ];
    CP := ng eq 0 select [ [Integers()|] ] else [ [c[i] : i in [1..ng]] : c in CartesianProduct(ranges) ];
    nloc := 0;
    for tup in CP do
        g := G!0;
        for i in [1..ng] do g +:= tup[i]*G.i; end for;
        Q := mw(g);
        P := 0; okq := true;
        try P := mEi(Q); catch e okq := false; end try;
        if not okq then continue; end if;
        if P[3] eq 0 then continue; end if;
        yv := P[1]/P[3];
        tv := 3 + (u-3)*yv^2;
        if tv in EXCL or tv eq u then continue; end if;
        if HeightOf(tv) gt 10^24 then continue; end if;
        nloc +:= 1;
        // C1, C2 hold by construction; test C3
        v3 := (tv-1)*(tv-9)*(u-1)*(u-9);
        if v3 eq 0 then continue; end if;
        if not IsSquare(v3) then continue; end if;
        pr := <tv, u>;
        if pr in survivors then continue; end if;
        Include(~survivors, pr);
        printf "FIBSURV t=%o u=%o y=%o rk=%o\n", tv, u, yv, rk;
    end for;
    npts +:= nloc;
    if nfib mod 10 eq 0 then
        printf "PROGRESS fiber %o fibers (%o ok), %o pts, %o survivors, %o s\n",
            nfib, nfibok, npts, #survivors, Cputime()-t0;
        System(Sprintf("echo 'fiberA-mw %o fibers %o ok %o pts %o surv %o s' >> ../logs/fiberAmw.progress", nfib, nfibok, npts, #survivors, Cputime()-t0));
    end if;
end for;
printf "FIBER_SCAN done: %o fibers (%o with MW), %o points tested, %o survivors, %o s\n",
    nfib, nfibok, npts, #survivors, Cputime()-t0;

// ---------- funnel survivors ----------
for pr in survivors do
    tv := pr[1]; uv := pr[2];
    if HeightOf(tv) gt 10^12 then printf "SKIPBIG t=%o\n", tv; continue; end if;
    EA := 0; EB := 0; okc := true;
    try EA := E26(tv); EB := E26(uv); catch e okc := false; end try;
    if not okc then continue; end if;
    if jInvariant(EA) eq jInvariant(EB) then printf "SKIPISO t=%o u=%o\n", tv, uv; continue; end if;
    printf "INSTANCE 2266fib t=%o u=%o\n", tv, uv;
    GlueFunnel(EA, EB, Sprintf("l2266fib|t=%o|u=%o", tv, uv), [Integers()|3,3], ~seen, ~nstat);
end for;
printf "SEARCH_DONE 2266fiber survivors=%o curves=%o aborts=%o known=%o hits=%o fails=%o %o s\n",
    #survivors, #seen, nstat["abort"], nstat["known"], nstat["hit"], nstat["fail"], Cputime()-t0;
quit;
