// lane_1212_fiber.m — [12,12] / [10,10] structured cross attack (plan
// section 2, steps 2-4).  For each base t0, the Delta-matched partners u are
// the rational points of the genus-1 fiber  w^2 = D(t0)*D(u)  (D = the
// squarefree class of the 2-division discriminant of the X1(N) curve), which
// has the rational point (u,w) = (t0, D(t0)).  We enumerate each fiber's
// Mordell-Weil points, test the validated beta-match halving condition
// (lane_1212.m: checksums + negative controls + in-sweep [8,8]/[4,8]
// positive/negative funnels all PASS), and funnel survivors.
// Usage: cd product/code && magma -b lane_1212_fiber.m > ../logs/lane1212_fiber.log
//   optional: HT:=<int> (default 13) height bound for t0; NB:=<int> (default 8)
SetColumns(0);
if not assigned HT then HT := 13; elif Type(HT) eq MonStgElt then HT := StringToInteger(HT); end if;
if not assigned NB then NB := 8; elif Type(NB) eq MonStgElt then NB := StringToInteger(NB); end if;
if not assigned MemGB then MemGB := 5; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
SetClassGroupBounds("GRH");   // speeds the per-fiber descents
load "split_lab.m";  // run from product/code/

RQ := Rationals();
QT<t> := FunctionField(RQ);
QY<U> := PolynomialRing(RQ);

BC := AssociativeArray();
BC[10] := [* t^3*(t-1)*(2*t-1)/(t^2-3*t+1)^2, -t*(t-1)*(2*t-1)/(t^2-3*t+1) *];
BC[12] := [* t*(2*t-1)*(3*t^2-3*t+1)*(2*t^2-2*t+1)/(t-1)^4, -t*(2*t-1)*(3*t^2-3*t+1)/(t-1)^3 *];
// squarefree-class discriminant polynomials (from lane_1212_sym.log, x4)
DP := AssociativeArray();
DP[10] := 8*U^3 - 8*U^2 + 1;
DP[12] := 12*U^4 - 24*U^3 + 20*U^2 - 8*U + 1;

function SFrat(x)
    n := Numerator(x)*Denominator(x);
    s := Sign(n); n := Abs(n);
    a := SquarefreeFactorization(n);
    return s*a;
end function;

function TryCurve(bv, cv)
    if bv eq 0 then return false, 0; end if;
    E := 0;
    try
        E := EllipticCurve([1-cv, -bv, -bv, 0, 0]);
    catch e
        return false, 0;
    end try;
    return true, E;
end function;

function DeltaData(N, tv, pair)
    ok := true; bv := 0; cv := 0;
    try
        bv := Evaluate(pair[1], tv); cv := Evaluate(pair[2], tv);
    catch e
        ok := false;
    end try;
    if not ok then return false, 0; end if;
    ok, E := TryCurve(bv, cv);
    if not ok then return false, 0; end if;
    P := 0;
    try P := E![0,0]; catch e return false, 0; end try;
    if Order(P) ne N then return false, 0; end if;
    E2, phi := WeierstrassModel(E);
    f := RQx!(HyperellipticPolynomials(E2));
    rts := Roots(f);
    if #rts eq 3 then return false, 0; end if;
    error if #rts ne 1, "unexpected root count";
    erat := rts[1][1];
    q := f div (RQx.1 - erat);
    q1 := Coefficient(q,1); q0 := Coefficient(q,0);
    D := q1^2 - 4*q0;
    d := SFrat(D);
    error if d eq 1, "square disc with one rational root?";
    okw, w := IsSquare(D/d);
    error if not okw, "D/d not square";
    xg := (phi(P))[1];
    alpha := xg - erat;
    error if alpha eq 0, "x(g1) hits e_rat";
    rec := < tv, E, d, alpha, xg + q1/2, w/2, erat + q1/2, w/2, jInvariant(E), Evaluate(q, erat) >;
    return true, rec;
end function;

function KMul(A1,B1,A2,B2,d)
    return A1*A2 + B1*B2*d, A1*B2 + A2*B1;
end function;

function HalvingMatch(recA, recB)
    d := recA[3];
    error if recB[3] ne d, "d mismatch";
    K := QuadraticField(d);
    A1 := recA[5]; B1 := recA[6]; A2 := recB[5]; B2 := recB[6];
    XA, XB := KMul(A1,B1,A2,B2,d);
    YA, YB := KMul(A1,B1,A2,-B2,d);
    return IsSquare(K!XA + XB*K.1), IsSquare(K!YA + YB*K.1);
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

ODDT := AssociativeArray();
ODDT[10] := [Integers()|5,5]; ODDT[12] := [Integers()|3,3];
HeightOf := func<q | Max(Abs(Numerator(q)), Abs(Denominator(q)))>;

t0 := Cputime();
for N in [12,10] do
    Dpol := DP[N];
    pairBC := BC[N];
    pairseen := {};
    nfib := 0; nmw := 0; ntest := 0; nsurv := 0;
    for tv in HeightRats(HT) do
        okA, recA := DeltaData(N, tv, pairBC);
        if not okA then continue; end if;
        Dt := Evaluate(Dpol, tv);
        if Dt eq 0 then continue; end if;
        quart := Dt*Dpol;      // in QY; w^2 = Dt*D(u), point (tv, Dt)
        if Discriminant(quart) eq 0 then continue; end if;
        nfib +:= 1;
        C := HyperellipticCurve(quart);
        pt := C![tv, Dt];
        E := 0; mE := 0; ok := true;
        try E, mE := EllipticCurve(C, pt); catch e ok := false; end try;
        if not ok then continue; end if;
        G := 0; mw := 0;
        try G, mw := MordellWeilGroup(E); catch e ok := false; end try;
        if not ok then printf "MWFAIL N=%o t=%o\n", N, tv; continue; end if;
        nmw +:= 1;
        mEi := 0;
        try mEi := Inverse(mE); catch e ok := false; end try;
        if not ok then continue; end if;
        ng := Ngens(G);
        ranges := [ Order(G.i) eq 0 select [-NB..NB] else [0..Order(G.i)-1] : i in [1..ng] ];
        CP := ng eq 0 select [ [Integers()|] ] else [ [c[i] : i in [1..ng]] : c in CartesianProduct(ranges) ];
        for tup in CP do
            g := G!0;
            for i in [1..ng] do g +:= tup[i]*G.i; end for;
            Q := mw(g);
            P := 0; okq := true;
            try P := mEi(Q); catch e okq := false; end try;
            if not okq or P[3] eq 0 then continue; end if;
            uv := P[1]/P[3];
            if uv eq tv or HeightOf(uv) gt 10^10 then continue; end if;
            pr := tv lt uv select <tv,uv> else <uv,tv>;
            if pr in pairseen then continue; end if;
            Include(~pairseen, pr);
            okB, recB := DeltaData(N, uv, pairBC);
            if not okB then continue; end if;
            if recA[9] eq recB[9] then continue; end if;   // isomorphic
            if recA[3] ne recB[3] then continue; end if;   // paranoia
            ntest +:= 1;
            if SFrat(recA[4]) ne SFrat(recB[4]) then continue; end if;  // alpha prefilter
            v1, v2 := HalvingMatch(recA, recB);
            if not (v1 or v2) then continue; end if;
            nsurv +:= 1;
            printf "SURVIVOR N=%o t=%o u=%o d=%o V1=%o V2=%o\n", N, tv, uv, recA[3], v1, v2;
            if HeightOf(uv) gt 10^7 then printf "SKIPBIG u=%o\n", uv; continue; end if;
            GlueFunnel(recA[2], recB[2], Sprintf("l1212fib|N=%o|t=%o|u=%o", N, tv, uv), ODDT[N], ~seen, ~nstat);
        end for;
        if nfib mod 10 eq 0 then
            printf "PROGRESS N=%o fibers=%o mw=%o tested=%o surv=%o %o s\n", N, nfib, nmw, ntest, nsurv, Cputime()-t0;
            System(Sprintf("echo 'fiberB-mw N=%o %o fibers %o mw %o tested %o surv %o s' >> ../logs/fiberBmw.progress", N, nfib, nmw, ntest, nsurv, Cputime()-t0));
        end if;
    end for;
    printf "FIBERJOIN N=%o fibers=%o mw-ok=%o beta-tested=%o survivors=%o\n", N, nfib, nmw, ntest, nsurv;
end for;
printf "SEARCH_DONE 1212fiber curves=%o aborts=%o known=%o hits=%o fails=%o %o s\n",
    #seen, nstat["abort"], nstat["known"], nstat["hit"], nstat["fail"], Cputime()-t0;
quit;
