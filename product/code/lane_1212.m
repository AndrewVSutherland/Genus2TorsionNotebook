// lane_1212.m — [12,12] / [10,10] via the g-halving delta-match condition
// (2026-08-13 plan, section 2), with X1(N) pools for N in {8,10,12}.
//
// Theory (derived+validated 2026-08-13): for a class-(b) gluing of two
// X1(N)-curves (one rational 2-torsion point T_rat = (N/2)*g1 each, conjugate
// pair over the same quadratic field K = Q(sqrt(d))), every order-144 (resp.
// order-100) enlargement of the [6,12] (resp. [5,10]) image comes from a
// rational half of the class im((N/2 odd part)...) = im(3g1,3g2)-type, and its
// descent condition collapses (delta(kg)=delta(g) for odd k; delta(T_rat)
// trivial for 4|N) to the SINGLE condition
//     beta_t == beta_u  or  beta_t == conj(beta_u)   in K*/(K*)^2,
// where beta = x(g1) - e_conj on a y^2 = cubic model.  The two variants
// correspond to the two rational gluings psi.  The rational-slot condition
// alpha_t == alpha_u mod squares (alpha = x(g1) - e_rat) is implied, and is
// used as a hash-join prefilter.
//   N=12: survivor pairs give [12,12]  (order 144; [2,6,12] is impossible:
//         A[2](Q) has rank <= 2 for class-(b) gluings).
//   N=10: T_rat = 5*g1 gives delta(T_rat) == delta(g1), so the SAME condition
//         is the 2-rank gain condition; survivors give [10,10] (order 100).
//   N=8:  [8,8] (already realized 2026-08-12) -- used as POSITIVE CONTROL:
//         the x1fam [8,8] hit pairs must satisfy the condition.
// Negative controls: the four exact-[6,12] x1fam pairs must fail both variants.
//
// Checksums per curve: for N in {8,12}: q(e_rat) is a rational square and
// kappa = e_rat - e_conj is a square in K (since T_rat = 2*(N/4)g1 in E(Q));
// for N=10: q(e_rat)*alpha is a square and kappa*beta is a square in K.
//
// DIAGONAL stage (added after the self-gluing insight, same session): a
// class-(b) curve glued to ITSELF via the conjugation psi (T_c <-> T_c-bar;
// valid since Aut(E) = +-1 acts trivially on E[2]) has halving condition
// beta == conj(beta) in K*/(K*)^2, i.e. N(beta) = alpha (mod squares) a square
// in K:  alpha == 1  or  alpha == d  mod rational squares.  ONE condition in
// one variable -> the diagonal is scanned directly.
// Usage: cd product/code && magma -b lane_1212.m > ../logs/lane1212_H120.log
//   optional: H:=<int> (default 120), MaxFunnel:=<int> (default 150)
SetColumns(0);
if not assigned H then H := 120; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned MaxFunnel then MaxFunnel := 150; elif Type(MaxFunnel) eq MonStgElt then MaxFunnel := StringToInteger(MaxFunnel); end if;
SetMemoryLimit(MemGB*10^9);
load "split_lab.m";  // run from product/code/

RQ := Rationals();
QT<t> := FunctionField(RQ);

// Kubert parametrizations (candidate variants, runtime-verified below)
CANDS := AssociativeArray();
CANDS[8]  := [* [* (2*t-1)*(t-1), (2*t-1)*(t-1)/t *] *];
CANDS[10] := [*
  [* t^3*(t-1)*(2*t-1)/(t^2-3*t+1)^2, -t*(t-1)*(2*t-1)/(t^2-3*t+1) *],
  [* t^3*(t-1)*(2*t-1)/(t^2-3*t+1)^2,  t*(t-1)*(2*t-1)/(t^2-3*t+1) *],
  [* -t^3*(t-1)*(2*t-1)/(t^2-3*t+1)^2, -t*(t-1)*(2*t-1)/(t^2-3*t+1) *]
*];
CANDS[12] := [*
  [* -t*(2*t-1)*(3*t^2-3*t+1)*(2*t^2-2*t+1)/(t-1)^4, -t*(2*t-1)*(3*t^2-3*t+1)/(t-1)^3 *],
  [*  t*(2*t-1)*(3*t^2-3*t+1)*(2*t^2-2*t+1)/(t-1)^4,  t*(2*t-1)*(3*t^2-3*t+1)/(t-1)^3 *],
  [*  t*(2*t-1)*(3*t^2-3*t+1)*(2*t^2-2*t+1)/(t-1)^4, -t*(2*t-1)*(3*t^2-3*t+1)/(t-1)^3 *],
  [* -t*(2*t-1)*(3*t^2-3*t+1)*(2*t^2-2*t+1)/(t-1)^4,  t*(2*t-1)*(3*t^2-3*t+1)/(t-1)^3 *]
*];

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

function VerifyFamily(N, pair)
    nok := 0;
    for tv in [RQ| 5, 7/2, -3/4 ] do
        ok := true; bv := 0; cv := 0;
        try
            bv := Evaluate(pair[1], tv); cv := Evaluate(pair[2], tv);
        catch e
            ok := false;
        end try;
        if not ok then continue; end if;
        ok, E := TryCurve(bv, cv);
        if not ok then continue; end if;
        P := E![0,0];
        if Order(P) eq N then nok +:= 1; end if;
    end for;
    return nok ge 2;
end function;

FAM := AssociativeArray();
for N in [8,10,12] do
    for pair in CANDS[N] do
        if VerifyFamily(N, pair) then
            FAM[N] := pair;
            printf "FAMILY N=%o verified\n", N;
            break;
        end if;
    end for;
    error if not IsDefined(FAM, N), "family verification failed";
end for;

function SFrat(x)
    n := Numerator(x)*Denominator(x);
    s := Sign(n); n := Abs(n);
    a := SquarefreeFactorization(n);
    return s*a;
end function;

// per-curve delta data on a y^2 = monic cubic model with exactly one rational
// 2-torsion point.  returns ok, rec
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
    if #rts eq 3 then return false, 0; end if;   // class (a) specialization
    error if #rts ne 1, "unexpected root count";
    erat := rts[1][1];
    q := f div (RQx.1 - erat);
    q1 := Coefficient(q,1); q0 := Coefficient(q,0);
    D := q1^2 - 4*q0;
    d := SFrat(D);
    error if d eq 1, "square disc with one rational root?";
    okw, w := IsSquare(D/d);
    error if not okw, "D/d not a square?";
    xg := (phi(P))[1];
    alpha := xg - erat;
    error if alpha eq 0, "x(g1) hits e_rat";
    // beta  = x(g1) - e_conj = (xg + q1/2) + (w/2) sqrt(d)
    // kappa = e_rat - e_conj = (erat + q1/2) + (w/2) sqrt(d)
    rec := < tv, E, d, alpha, xg + q1/2, w/2, erat + q1/2, w/2, jInvariant(E), Evaluate(q, erat) >;
    return true, rec;
end function;

// K-square test for z = A + B*sqrt(d), given K
function IsSqK(K, A, B)
    return IsSquare(K!A + B*K.1);
end function;

Kcache := AssociativeArray();
function GetK(d)
    // (associative-array cache is per-call read-only; rebuild is cheap anyway)
    return QuadraticField(d);
end function;

// product (A1+B1 s)(A2+B2 s), s^2 = d
function KMul(A1,B1,A2,B2,d)
    return A1*A2 + B1*B2*d, A1*B2 + A2*B1;
end function;

// ---------- per-curve checksums ----------
procedure Checksums(N, rec)
    d := rec[3]; K := GetK(d);
    qe := rec[10];
    if N mod 4 eq 0 then
        error if not IsSquare(qe), Sprintf("checksum q(erat) N=%o t=%o", N, rec[1]);
        error if not IsSqK(K, rec[7], rec[8]), Sprintf("checksum kappa N=%o t=%o", N, rec[1]);
    else
        error if not IsSquare(qe*rec[4]), Sprintf("checksum q(erat)*alpha N=%o t=%o", N, rec[1]);
        A,B := KMul(rec[7],rec[8],rec[5],rec[6],d);
        error if not IsSqK(K, A, B), Sprintf("checksum kappa*beta N=%o t=%o", N, rec[1]);
    end if;
end procedure;

// ---------- the halving condition for a pair ----------
function HalvingMatch(recA, recB)
    // returns matchV1, matchV2 (psi-signs)
    d := recA[3];
    error if recB[3] ne d, "bucket mismatch";
    K := GetK(d);
    A1 := recA[5]; B1 := recA[6]; A2 := recB[5]; B2 := recB[6];
    XA, XB := KMul(A1,B1,A2,B2,d);
    YA, YB := KMul(A1,B1,A2,-B2,d);
    return IsSqK(K, XA, XB), IsSqK(K, YA, YB);
end function;

// ---------- CONTROLS ----------
// NOTE: the x1fam [8,8]-hit pairs turned out to be class-(a) specializations
// ([2,8]-torsion, "Q"-bucket) -- not usable for the class-(b) machinery.
// Positive validation is done IN-SWEEP instead: N=8 class-(b) pairs that
// beta-match are funneled and must come out exactly [8,8] (a KNOWN group, so
// they print EXACTKNOWN [8,8]); unmatched D-matched N=8 pairs are funneled as
// in-sweep negative controls and must come out strictly smaller.
printf "== controls ==\n";
// negative: x1fam exact-[6,12] pairs (verified class (b))
NEG := [ [RQ|2/17, -3/7], [RQ|1/6, -1/6], [RQ|1/7, -2], [RQ|-3, -3/5] ];
for pr in NEG do
    ok1, rA := DeltaData(12, pr[1], FAM[12]);
    ok2, rB := DeltaData(12, pr[2], FAM[12]);
    error if not (ok1 and ok2), "control pool build failed";
    Checksums(12, rA); Checksums(12, rB);
    error if rA[3] ne rB[3], "control d not identical";
    v1, v2 := HalvingMatch(rA, rB);
    printf "NEGCTRL N=12 (t,u)=(%o,%o) V1=%o V2=%o\n", pr[1], pr[2], v1, v2;
    error if v1 or v2, "NEGATIVE CONTROL FAILED (theory bug!)";
end for;
printf "CONTROLS PASS\n";

// ---------- pool build ----------
function HeightRats(H)
    S := [];
    for a in [1..H], b in [1..H] do
        if GCD(a,b) eq 1 then Append(~S, a/b); Append(~S, -a/b); end if;
    end for;
    return S;
end function;

pool := AssociativeArray();   // N -> list of recs
jseen := {};
for N in [8,10,12] do
    lst := [* *];
    nclassa := 0;
    for tv in HeightRats(H) do
        ok, rec := DeltaData(N, tv, FAM[N]);
        if not ok then continue; end if;
        if rec[9] in jseen then continue; end if;
        Include(~jseen, rec[9]);
        Checksums(N, rec);
        Append(~lst, rec);
    end for;
    pool[N] := lst;
    printf "POOL N=%o: %o curves\n", N, #lst;
end for;

// ---------- bucket join and sweep ----------
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

ODDT := AssociativeArray();
ODDT[8] := [Integers()|]; ODDT[10] := [Integers()|5,5]; ODDT[12] := [Integers()|3,3];

nf := 0; t0 := Cputime();

// ---------- diagonal stage: self-gluing via conjugation psi ----------
// (N=8 included as positive control: diagonal survivors must funnel to [8,8])
for N in [8,10,12] do
    ndiag := 0;
    for rec in pool[N] do
        al := rec[4]; d := rec[3];
        s1 := IsSquare(al); s2 := IsSquare(al*d);
        if not (s1 or s2) then continue; end if;
        ndiag +:= 1;
        printf "DIAGSURV N=%o t=%o d=%o alsq=%o aldsq=%o\n", N, rec[1], d, s1, s2;
        if nf ge MaxFunnel then continue; end if;
        nf +:= 1;
        GlueFunnel(rec[2], rec[2], Sprintf("l1212diag|N=%o|t=%o", N, rec[1]), ODDT[N], ~seen, ~nstat);
    end for;
    printf "DIAG N=%o survivors=%o\n", N, ndiag;
end for;

// ---------- cross-pair stage ----------
for N in [8,10,12] do
    lst := pool[N];
    // hash join on d only (alpha-match applied separately so the N=8 negative
    // controls can come from the same buckets)
    B := AssociativeArray();
    for i in [1..#lst] do
        key := lst[i][3];
        if IsDefined(B, key) then B[key] := B[key] cat [i]; else B[key] := [i]; end if;
    end for;
    npairs := 0; nsurv := 0; nnegf := 0;
    for key in Keys(B) do
        idxs := B[key];
        if #idxs lt 2 then continue; end if;
        for ii in [1..#idxs] do
            for jj in [ii+1..#idxs] do
                rA := lst[idxs[ii]]; rB := lst[idxs[jj]];
                if rA[9] eq rB[9] then continue; end if;   // same j
                npairs +:= 1;
                match := false;
                if SFrat(rA[4]) eq SFrat(rB[4]) then      // alpha prefilter
                    v1, v2 := HalvingMatch(rA, rB);
                    match := v1 or v2;
                    if match then
                        nsurv +:= 1;
                        printf "SURVIVOR N=%o t=%o u=%o d=%o V1=%o V2=%o\n", N, rA[1], rB[1], rA[3], v1, v2;
                        if nf lt MaxFunnel then
                            nf +:= 1;
                            GlueFunnel(rA[2], rB[2], Sprintf("l1212|N=%o|t=%o|u=%o", N, rA[1], rB[1]), ODDT[N], ~seen, ~nstat);
                        end if;
                    end if;
                end if;
                // in-sweep negative controls: a few unmatched N=8 pairs
                if N eq 8 and not match and nnegf lt 3 and nf lt MaxFunnel then
                    nnegf +:= 1; nf +:= 1;
                    printf "NEGFUNNEL N=8 t=%o u=%o d=%o\n", rA[1], rB[1], rA[3];
                    GlueFunnel(rA[2], rB[2], Sprintf("l1212negctrl|t=%o|u=%o", rA[1], rB[1]), ODDT[N], ~seen, ~nstat);
                end if;
            end for;
        end for;
    end for;
    printf "JOIN N=%o d-matched pairs=%o survivors=%o\n", N, npairs, nsurv;
end for;
printf "SEARCH_DONE 1212 H=%o funneled=%o curves=%o aborts=%o known=%o hits=%o fails=%o %o s\n",
    H, nf, #seen, nstat["abort"], nstat["known"], nstat["hit"], nstat["fail"], Cputime()-t0;
quit;
