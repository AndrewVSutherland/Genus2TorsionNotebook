
//////////////////////////////////////////////////////////////////////
//  Per-R local-solvability scan for the P_R-halving fiber product in
//  M_1(8,4), on the Pell-parametrized C1 locus.
//
//  From agent_m18_416_sstage_symbolic/crossproducts: at fixed R, the
//  halving condition at m is
//    S_A: V_A = 2*alpha_A(m) +- 2*h_A(m)*sqrt(G(m)) contains a square
//    S_B: same with B data,
//  and G = square is implied by either (exactly (V-2a)^2 = 4Gh^2).
//
//  For each fiber R (taken from the C1&C2 survivor list, regenerated
//  here) and each prime p in plist (+ infinity), we test whether some
//  m in Q_p satisfies S_A and S_B simultaneously, by enumerating
//  representatives m = a + b*p (units+first digit), m = c/p, m = c*p
//  (mod p^2 classes; mod 16/32 for p=2) and doing exact p-adic square
//  tests.  A fiber with NO obstructed prime is Everywhere Locally
//  Solvable (up to the tested places) -- a candidate for a global point
//  and Mordell-Weil analysis.  A uniform obstruction pattern across
//  fibers would point at a nonexistence proof.
//
//  Usage: magma -b hR:=30 hM:=30 maxfibers:=60 agent_m18_416_fiber_local_scan.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

if not assigned hR then hR := 30;
elif Type(hR) eq MonStgElt then hR := StringToInteger(hR); end if;
if not assigned hM then hM := 30;
elif Type(hM) eq MonStgElt then hM := StringToInteger(hM); end if;
if not assigned maxfibers then maxfibers := 60;
elif Type(maxfibers) eq MonStgElt then maxfibers := StringToInteger(maxfibers); end if;

plist := [2,3,5,7,11,13,17,19,23];

// ---- symbolic fiber data over Q(R,m) ----
K2<R,m> := RationalFunctionField(Q, 2);
PX<x> := PolynomialRing(K2);
K := -2*R*(R^2-1);
w := (m^2 + K)/(m^2 - K);
W := w^2;
t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
c4 := R + 2 + 4*t;
Apol := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
Bpol := c4*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);
XR := -c4*R;
At := PX![K2!co : co in Coefficients(c4^2*Evaluate(Apol, x/c4))];
Bt := PX![K2!co : co in Coefficients(c4*Evaluate(Bpol, x/c4))];
G := 2*(R^2-1)*(R*(2*R+1) - W*(R+2));
alphaA := XR + Coefficient(At,1)/2;
alphaB := XR + Coefficient(Bt,1)/2;
dA := Discriminant(At);
dB := Discriminant(Bt);
okA, hA := IsSquare(Evaluate(At, XR)/G); assert okA;
okB, hB := IsSquare(Evaluate(Bt, XR)/G); assert okB;

// numeric evaluation helpers
funcs := [G, alphaA, hA, alphaB, hB, dA, dB];
fnum := [Numerator(f) : f in funcs];
fden := [Denominator(f) : f in funcs];

function EvalAll(Rv, mv)
    vals := [];
    for i in [1..#funcs] do
        dd := Evaluate(fden[i], [Rv, mv]);
        if dd eq 0 then return false, []; end if;
        Append(~vals, Evaluate(fnum[i], [Rv, mv])/dd);
    end for;
    return true, vals;   // [G, aA, hA, aB, hB, dA, dB]
end function;

// degenerate-m guard: denominators of w, curve nondegeneracy handled by
// EvalAll failure + G=0 checks.

// p-adic joint test at a specific rational m
function PointPasses(Rv, mv, Qp)
    ok, v := EvalAll(Rv, mv);
    if not ok then return false; end if;
    Gv := v[1]; aA := v[2]; hAv := v[3]; aB := v[4]; hBv := v[5];
    if Gv eq 0 or hAv eq 0 or hBv eq 0 then return false; end if;
    okG, sg := IsSquare(Qp!Gv);
    if not okG then return false; end if;
    passA := false;
    for sgn in [1,-1] do
        VA := Qp!(2*aA) + sgn*(Qp!(2*hAv))*sg;
        if VA ne 0 and IsSquare(VA) then passA := true; break; end if;
    end for;
    if not passA then return false; end if;
    for sgn in [1,-1] do
        VB := Qp!(2*aB) + sgn*(Qp!(2*hBv))*sg;
        if VB ne 0 and IsSquare(VB) then return true; end if;
    end for;
    return false;
end function;

// local solvability of fiber R at p: try representative m classes
function FiberSolvableAt(Rv, p)
    Qp := pAdicField(p, 60);
    reps := [];
    if p eq 2 then
        for a in [1..31 by 2] do Append(~reps, Q!a); end for;          // units mod 32
        for a in [1..15] do Append(~reps, Q!(2*a)); end for;           // v=1..
        for a in [1..31 by 2] do Append(~reps, Q!a/2); end for;        // v=-1
        for a in [1..31 by 2] do Append(~reps, Q!a/4); end for;        // v=-2
        for a in [1..15 by 2] do Append(~reps, Q!(4*a)); end for;      // v=2
    else
        for a in [1..p-1] do for b in [0..p-1] do
            Append(~reps, Q!(a + b*p));
        end for; end for;
        for c in [1..p-1] do Append(~reps, Q!c/p); end for;
        for c in [1..p-1] do Append(~reps, Q!(c*p)); end for;
    end if;
    for mv in reps do
        if PointPasses(Rv, mv, Qp) then return true, mv; end if;
    end for;
    return false, 0;
end function;

// real place: sign conditions on a sample grid
function FiberSolvableAtInfinity(Rv)
    for num in [-400..400] do
        mv := Q!num/7;
        if mv eq 0 then continue; end if;
        ok, v := EvalAll(Rv, mv);
        if not ok then continue; end if;
        Gv := v[1]; aA := v[2]; hAv := v[3]; aB := v[4]; hBv := v[5];
        if Gv le 0 or hAv eq 0 or hBv eq 0 then continue; end if;
        sg := Sqrt(RealField(30)!Gv);
        passA := false;
        for sgn in [1,-1] do
            if RealField(30)!(2*aA) + sgn*RealField(30)!(2*hAv)*sg gt 0 then passA := true; break; end if;
        end for;
        if not passA then continue; end if;
        for sgn in [1,-1] do
            if RealField(30)!(2*aB) + sgn*RealField(30)!(2*hBv)*sg gt 0 then return true; end if;
        end for;
    end for;
    return false;
end function;

// ---- regenerate C1&C2 survivors to get the fiber list ----
function IsSquareQ(qv)
    qv := Q!qv;
    if qv le 0 then return false; end if;
    return IsSquare(Numerator(qv)) and IsSquare(Denominator(qv));
end function;
function RationalParametersOfHeight(Bd)
    vals := [];
    for den in [1..Bd] do for num in [-Bd..Bd] do
        if GCD(num, den) eq 1 then Append(~vals, Q!num/den); end if;
    end for; end for;
    return Sort(Setseq(Seqset(vals)));
end function;

printf "collecting C1&C2 survivors at hR=%o hM=%o ...\n", hR, hM;
Rparams := RationalParametersOfHeight(hR);
Mparams := RationalParametersOfHeight(hM);
fiberR := [];   // distinct R with a C2 point, with an example m
fiberSeen := {};
nsurv := 0;
for Rv in Rparams do
    if Rv in {Q!0, Q!1, Q!-1} then continue; end if;
    Kv := -2*Rv*(Rv^2 - 1);
    for mv in Mparams do
        den := mv^2 - Kv;
        if den eq 0 then continue; end if;
        wv := (mv^2 + Kv)/den;
        if wv in {Q!0, Q!1, Q!-1} then continue; end if;
        Wv := wv^2;
        Gv := 2*(Rv^2-1)*(Rv*(2*Rv+1) - Wv*(Rv+2));
        if not IsSquareQ(Gv) then continue; end if;
        nsurv +:= 1;
        if Rv notin fiberSeen then
            Include(~fiberSeen, Rv);
            Append(~fiberR, <Rv, mv>);
        end if;
    end for;
end for;
printf "survivors=%o distinct fibers=%o (scanning up to %o)\n",
    nsurv, #fiberR, maxfibers;

// ---- the scan ----
nELS := 0; nObs := 0;
obsCount := AssociativeArray();
for i in [1..Min(#fiberR, maxfibers)] do
    Rv := fiberR[i][1]; mex := fiberR[i][2];
    obstructed := [];
    if not FiberSolvableAtInfinity(Rv) then Append(~obstructed, 0); end if;   // 0 = infinity
    for p in plist do
        ok, mwit := FiberSolvableAt(Rv, p);
        if not ok then Append(~obstructed, p); end if;
    end for;
    if #obstructed eq 0 then
        nELS +:= 1;
        printf "FIBER R=%o (example m=%o): ELS up to p<=%o  <== global candidate\n",
            Rv, mex, plist[#plist];
    else
        nObs +:= 1;
        printf "FIBER R=%o: obstructed at %o\n", Rv,
            [pv eq 0 select "oo" else Sprint(pv) : pv in obstructed];
        for pv in obstructed do
            key := pv eq 0 select "oo" else Sprint(pv);
            if IsDefined(obsCount,key) then obsCount[key] +:= 1;
            else obsCount[key] := 1; end if;
        end for;
    end if;
end for;
printf "\nSUMMARY: fibers scanned=%o ELS=%o obstructed=%o\n",
    Min(#fiberR,maxfibers), nELS, nObs;
print "obstruction place histogram:";
for key in Sort([k : k in Keys(obsCount)]) do
    printf "  %o : %o\n", key, obsCount[key];
end for;
print "DONE";
quit;
