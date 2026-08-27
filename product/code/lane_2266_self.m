// lane_2266_self.m — [2,2,6,6] via SELF-gluing E26(t) ~ E26(t) along a
// 3-cycle psi (valid: Aut(E)=+-1 acts trivially on E[2], so psi != id is
// never induced by an isomorphism; sigma = id is the degenerate product and
// the transposition psi's carry a -1 obstruction).
// Conditions (the two 3-cycle sigma-systems evaluated on the diagonal t=u,
// from the validated lane_2266_derive.m tables; 2 independent):
//   C1: -2(t-3)(t+3)(t-5) = square
//   C2: -(t+3)(t-5)(t-1)(t-9) = square
//   (C3: 2(t-3)(t-1)(t-9) = C1*C2 mod squares; any two imply the third)
// A survivor t has FULL rational J[2] on the self-glued Jacobian, odd part
// [3,3]: J(Q)_tors >= [2,2,6,6], order 144.
// Strategy: the simultaneous locus is a genus-3 bidouble cover; enumerate
// Mordell-Weil points of the genus-1 subcovers and cross-test.
// Usage: cd product/code && magma -b lane_2266_self.m > ../logs/lane2266_self.log
SetColumns(0);
SetMemoryLimit(4*10^9);
load "split_lab.m";  // run from product/code/

RQ := Rationals();
QY<Y> := PolynomialRing(RQ);

c1p := -2*(Y-3)*(Y+3)*(Y-5);          // cubic
c2p := -(Y+3)*(Y-5)*(Y-1)*(Y-9);      // quartic
c3p := 2*(Y-3)*(Y-1)*(Y-9);           // cubic

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

EXCL := {RQ|3,-3,1,5,9};
HeightOf := func<q | Max(Abs(Numerator(q)), Abs(Denominator(q)))>;

survivors := {@ RQ | @};

// try a t-value: needs C1 and C2 (test both; enumeration guarantees one)
procedure TryT(tv, src, ~survivors)
    if tv in EXCL then return; end if;
    v1 := Evaluate(c1p, tv); v2 := Evaluate(c2p, tv);
    if v1 eq 0 or v2 eq 0 then return; end if;
    if IsSquare(v1) and IsSquare(v2) then
        if tv in survivors then return; end if;
        Include(~survivors, tv);
        printf "SELFSURV t=%o height=%o src=%o\n", tv, HeightOf(tv), src;
    end if;
end procedure;

// ---------- subcover analysis ----------
// enumerate MW points of a genus-1 curve given by y^2 = f (deg 3 or 4)
procedure EnumCover(fp, name, NB, ~survivors)
    C := HyperellipticCurve(fp);
    basepts := Points(C : Bound := 10000);
    printf "== %o: y^2 = %o ==\n", name, fp;
    printf "%o small points: %o\n", name, #basepts;
    if #basepts eq 0 then
        printf "%o: NO small rational points (possible torsor obstruction)\n", name;
        return;
    end if;
    // scan the small points directly
    for P in basepts do
        if P[3] ne 0 then TryT(P[1]/P[3], name cat "-small", ~survivors); end if;
    end for;
    E, mE := EllipticCurve(C, Rep(basepts));
    Em, mmin := MinimalModel(E);
    printf "%o: E = %o\n", name, aInvariants(Em);
    G, mw := MordellWeilGroup(Em);
    printf "%o: MW group invariants %o\n", name, Invariants(G);
    mEi := Inverse(mE);
    mmini := Inverse(mmin);
    ng := Ngens(G);
    ranges := [ Order(G.i) eq 0 select [-NB..NB] else [0..Order(G.i)-1] : i in [1..ng] ];
    cnt := 0;
    CP := ng eq 0 select [ [Integers()|] ] else [ [c[i] : i in [1..ng]] : c in CartesianProduct(ranges) ];
    for tup in CP do
        g := G!0;
        for i in [1..ng] do g +:= tup[i]*G.i; end for;
        Q := mw(g);
        Qe := 0; okq := true;
        try Qe := mmini(Q); catch e okq := false; end try;
        if not okq then continue; end if;
        P := 0;
        try P := mEi(Qe); catch e okq := false; end try;
        if not okq then continue; end if;
        if P[3] eq 0 then continue; end if;
        tv := P[1]/P[3];
        if HeightOf(tv) gt 10^40 then continue; end if;
        cnt +:= 1;
        TryT(tv, name cat "-mw", ~survivors);
    end for;
    printf "%o: enumerated %o affine MW points\n", name, cnt;
end procedure;

// G.i with Order 0 test needs G in scope; wrap generically (Magma: Order of
// free generator is 0).  NB box 14 per free generator.
EnumCover(c1p, "C1", 14, ~survivors);
EnumCover(c3p, "C3", 14, ~survivors);
EnumCover(c2p, "C2", 14, ~survivors);

// ---------- also a direct small-height scan (belt and braces) ----------
function HeightRats(H)
    S := [];
    for a in [1..H], b in [1..H] do
        if GCD(a,b) eq 1 then Append(~S, a/b); Append(~S, -a/b); end if;
    end for;
    return S;
end function;
for tv in HeightRats(200) do
    TryT(tv, "scan", ~survivors);
end for;
printf "SURVIVORS total %o\n", #survivors;

// ---------- funnel ----------
for tv in survivors do
    if HeightOf(tv) gt 10^12 then
        printf "SKIPBIG t=%o (height too large for auto-funnel)\n", tv;
        continue;
    end if;
    EA := 0; okc := true;
    try EA := E26(tv); catch e okc := false; end try;
    if not okc then continue; end if;
    printf "INSTANCE 2266self t=%o\n", tv;
    GlueFunnel(EA, EA, Sprintf("l2266self|t=%o", tv), [Integers()|3,3], ~seen, ~nstat);
end for;
printf "SEARCH_DONE 2266self survivors=%o curves=%o aborts=%o known=%o hits=%o fails=%o\n",
    #survivors, #seen, nstat["abort"], nstat["known"], nstat["hit"], nstat["fail"];
quit;
