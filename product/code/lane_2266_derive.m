// lane_2266_derive.m — symbolic Prop-12 slot machinery for the [2,2,6,6] and
// [2,2,2,24] routes (2026-08-13 plan, section 1).
//   * iota tables (2-descent images of the 2-torsion points) of the HLP
//     universal curves E26 = E^t_{2,6}, E28 = E^u_{2,8}, mod squares over Q(t);
//   * the sigma in S3 condition systems "a half of the graph class
//     (T_i, T'_sigma(i)) descends" -- per sigma the FULL two-gain system
//     collapses to 3 unordered-pair class-matching conditions;
//   * anchors: iota(T1) of E26 = HLP's printed value; sigma=id T1-system on
//     E26xE26 = HLP system (6); E28xE28 sigma=id = single condition (5);
//     E28 has a trivial iota row (it has rational 4-torsion); the [2,2,24]
//     witness satisfies a full T_i-subsystem of some E26xE28 sigma;
//   * extraction of the u-free subsystems of the E26xE28 sigmas => t-side
//     genus-1 quartics for the [2,2,2,24] route + a point search on each.
// Usage: cd product/code && magma -b lane_2266_derive.m > ../logs/derive2266.log
SetColumns(0);
SetMemoryLimit(3*10^9);

RQ := Rationals();
QT<t> := FunctionField(RQ);
Pt<T> := PolynomialRing(RQ);   // class representatives live here (T = the parameter)

// ---------- signed squarefree part of a nonzero rational ----------
function SFint(x)
    n := Numerator(x)*Denominator(x);
    s := Sign(n); n := Abs(n);
    a := SquarefreeFactorization(n);
    return s*a;
end function;

// ---------- mod-square class of a nonzero element of Q(t) ----------
// returns the squarefree representative in Pt: signed squarefree integer
// times the product of the odd-multiplicity monic irreducible factors.
// NOTE Factorization drops the unit -- reconstructed explicitly (magma-lab #1).
function SqClassFF(q)
    error if q eq 0, "SqClassFF: zero input";
    N := Numerator(q); D := Denominator(q);
    P := Pt!N * Pt!D;
    fac := Factorization(P);
    if #fac eq 0 then
        unit := P;
    else
        unit := P div &*[ f[1]^f[2] : f in fac ];
    end if;
    error if Degree(unit) ne 0, "unit not constant";
    c := RQ!Coefficient(unit, 0);
    cls := Pt!SFint(c);
    for f in fac do if IsOdd(f[2]) then cls *:= f[1]; end if; end for;
    return cls;
end function;

// q1 == q2 mod squares in Q(t)?  (q1, q2 nonzero elements of QT)
function SameClassFF(q1, q2)
    return SqClassFF(q1*q2) eq Pt!1;
end function;

// ---------- the universal curves (roots of the 2-division polynomial) ----------
// (from HLP Table 5 x-coordinates, runtime-verified in lane_hlp37.m)
r26 := [ (-2*t+10)/((t+3)*(t-3)),
         (-t^3+7*t^2-11*t+5)/(4*(t+3)*(t-3)^2),
         (-2*t^2+4*t-2)/((t+3)^2*(t-3)) ];
r28 := [ (16*t^3+12*t^2+2*t)/(8*t^2-1)^2,
         (32*t^3+24*t^2+8*t+1)/(16*t^2*(8*t^2-1)),
         (-32*t^4-32*t^3-12*t^2-2*t)/((4*t+1)^2*(8*t^2-1)) ];

// ---------- iota tables: row i = class vector of delta(T_i), slot k ----------
function IotaTable(rs)
    M := [];
    for i in [1..3] do
        row := [];
        for k in [1..3] do
            q := i eq k select &*[ rs[i]-rs[j] : j in [1..3] | j ne i ] else rs[i]-rs[k];
            Append(~row, SqClassFF(q));
        end for;
        Append(~M, row);
    end for;
    return M;
end function;

M26 := IotaTable(r26);
M28 := IotaTable(r28);
printf "== iota table E26 (classes mod squares, rows = T_i, cols = slots) ==\n";
for i in [1..3] do printf "  iota(T%o) = %o\n", i, M26[i]; end for;
printf "== iota table E28 ==\n";
for i in [1..3] do printf "  iota(T%o) = %o\n", i, M28[i]; end for;

// ---------- ANCHOR 1: iota(T1) of E26 equals HLP's printed value ----------
HLP1 := [ 2*(t-3)*(t+3)*(t-5), (t+3)*(t-5), 2*(t-3) ];
ok1 := forall{ k : k in [1..3] | SameClassFF(Evaluate(M26[1][k], t), HLP1[k]) };
printf "ANCHOR1 iota(T1) E26 vs HLP printed: %o\n", ok1 select "PASS" else "FAIL";
error if not ok1, "anchor 1 failed";

// ---------- ANCHOR 2: E28 has a trivial iota row (rational 4-torsion) ----------
trivrows := [ i : i in [1..3] | forall{ k : k in [1..3] | k eq i or M28[i][k] eq Pt!1 } ];
printf "ANCHOR2 E28 trivial iota rows: %o\n", trivrows;
error if #trivrows ne 1, "expected exactly one trivial row for E28";

// ---------- sigma condition systems ----------
// A half of the graph class (T_i, T'_sigma(i)) descends iff
//   iota_t(T_i)_k == iota_u(T'_sigma(i))_sigma(k) mod squares for all k.
// Any two of the three T_i-systems imply the third; the union over i reduces
// to 3 unordered-pair matching conditions:
//   class(rA_i - rA_k)(t) * class(rB_sigma(i) - rB_sigma(k))(u) is a square,
//   for {i,k} in {{1,2},{1,3},{2,3}}  (transposed orders flip both signs
//   simultaneously, so the i<k orientation is consistent).
SIGMAS := [ [1,2,3],[2,1,3],[3,2,1],[1,3,2],[2,3,1],[3,1,2] ];
PAIRS  := [ [1,2],[1,3],[2,3] ];

function SigmaSystem(rsA, rsB, sg)
    // sequence of <LA, LB> in Pt x Pt: condition is LA(t)*LB(u) square
    return [ < SqClassFF(rsA[p[1]]-rsA[p[2]]), SqClassFF(rsB[sg[p[1]]]-rsB[sg[p[2]]]) > : p in PAIRS ];
end function;

printf "== E26 x E26 sigma systems (condition: LA(t)*LB(u) square) ==\n";
for si in [1..6] do
    sys := SigmaSystem(r26, r26, SIGMAS[si]);
    printf "sigma=%o:\n", SIGMAS[si];
    for c in sys do printf "   LA = %o   |   LB = %o\n", c[1], c[2]; end for;
end for;

// ---------- ANCHOR 3: sigma=id T1-system on E26xE26 = HLP system (6) ----------
sysid := SigmaSystem(r26, r26, [1,2,3]);
ok3 := SameClassFF(Evaluate(sysid[1][1], t), (t+3)*(t-5)) and SameClassFF(Evaluate(sysid[1][2], t), (t+3)*(t-5))
   and SameClassFF(Evaluate(sysid[2][1], t), 2*(t-3))     and SameClassFF(Evaluate(sysid[2][2], t), 2*(t-3));
printf "ANCHOR3 sigma=id T1-system reproduces HLP (6): %o\n", ok3 select "PASS" else "FAIL";
error if not ok3, "anchor 3 failed";

// ---------- ANCHOR 4: E28xE28 sigma=id collapses to the single condition (5) ----------
D8v := (8*t^2-1)*(8*t^2+8*t+1);
sys88 := SigmaSystem(r28, r28, [1,2,3]);
nauto := 0; nD8 := 0;
for c in sys88 do
    if Degree(c[1]) eq 0 and Degree(c[2]) eq 0 and IsSquare(RQ!(Coefficient(c[1],0)*Coefficient(c[2],0))) then
        nauto +:= 1;
    elif SameClassFF(Evaluate(c[1], t), D8v) and SameClassFF(Evaluate(c[2], t), D8v) then
        nD8 +:= 1;
    end if;
end for;
printf "ANCHOR4 E28xE28 sigma=id: %o auto conditions, %o = condition (5): %o\n",
    nauto, nD8, (nauto eq 2 and nD8 eq 1) select "PASS" else "FAIL";
error if nauto ne 2 or nD8 ne 1, "anchor 4 failed";

// ---------- ANCHOR 5: the [2,2,24] witness satisfies a full T_i-subsystem ----------
// witness from lane_hlp37/verify_witnesses2: y=2/9 -> t0 = 3 - y^2/2, u0 = 1/3.
t0 := RQ!3 - (RQ!2/9)^2/2;  u0 := RQ!1/3;
printf "== E26 x E28 sigma systems at the [2,2,24] witness (t,u)=(%o,%o) ==\n", t0, u0;
function EvalCond(c, tv, uv)  // c = <LA,LB>; true iff LA(tv)*LB(uv) is a nonzero rational square
    va := Evaluate(c[1], tv); vb := Evaluate(c[2], uv);
    if va eq 0 or vb eq 0 then return false; end if;
    return IsSquare(va*vb);
end function;
found5 := false;
for si in [1..6] do
    sys := SigmaSystem(r26, r28, SIGMAS[si]);
    holds := [ EvalCond(c, t0, u0) : c in sys ];
    printf "sigma=%o holds=%o\n", SIGMAS[si], holds;
    for i in [1..3] do
        idx := [ j : j in [1..3] | i in PAIRS[j] ];
        if forall{ j : j in idx | holds[j] } then
            printf "   -> full T%o-subsystem HOLDS at witness (sigma=%o)\n", i, SIGMAS[si];
            found5 := true;
        end if;
    end for;
end for;
printf "ANCHOR5 [2,2,24] witness satisfies a full T_i-subsystem: %o\n", found5 select "PASS" else "FAIL";
error if not found5, "anchor 5 failed";

// ---------- [2,2,2,24] route: u-free subsystems of E26xE28 sigmas ----------
// E28's two trivial-class pairs make exactly 2 of the 3 conditions u-free for
// every sigma.  t-side system:  cB1*LA1(t) and cB2*LA2(t) both squares.
// Parametrize one condition, substitute into the other -> genus-1 quartic;
// point-search.  Reported t* feed the u-side D8-matching sweep in lane_2266.m.
printf "== [2,2,2,24] t-side systems and quartics ==\n";
tstars := [];   // tuples <si, tv>
EXCL26 := {RQ|3,-3,1,5,9};
for si in [1..6] do
    sys := SigmaSystem(r26, r28, SIGMAS[si]);
    ufree := [ j : j in [1..3] | Degree(sys[j][2]) eq 0 ];
    if #ufree ne 2 then printf "sigma=%o: %o u-free conditions (skip)\n", SIGMAS[si], #ufree; continue; end if;
    rep := [ j : j in [1..3] | not j in ufree ][1];
    cs  := [ RQ!Coefficient(sys[j][2], 0) : j in ufree ];
    LAs := [ sys[j][1] : j in ufree ];
    printf "sigma=%o: t-side: %o*(%o) and %o*(%o) both squares; u-cond: (%o)(t) * (%o)(u) square\n",
        SIGMAS[si], cs[1], LAs[1], cs[2], LAs[2], sys[rep][1], sys[rep][2];
    if Degree(LAs[1]) gt Degree(LAs[2]) then
        tmp := LAs[1]; LAs[1] := LAs[2]; LAs[2] := tmp;
        tmpc := cs[1]; cs[1] := cs[2]; cs[2] := tmpc;
    end if;
    // build the substitution t(m) parametrizing condition 1, in QT (t plays m)
    tsub := QT!0; okp := true;
    if Degree(LAs[1]) eq 1 then
        // c1*(a t + b) = m^2  ->  t = (m^2/c1 - b)/a
        a := Coefficient(LAs[1],1); b := Coefficient(LAs[1],0);
        tsub := (t^2/cs[1] - b)/a;
    else
        // conic c1*LA1(t) = z^2 through a rational root of LA1
        rts := Roots(LAs[1]);
        if #rts eq 0 then
            printf "   (deg-2 base condition has no rational root -- skipped)\n";
            okp := false;
        else
            rt := rts[1][1];
            lc := Coefficient(LAs[1],2);
            othr := #rts eq 2 select rts[2][1] else -Coefficient(LAs[1],1)/lc - rt;
            e := rt - othr;      // LA1 = lc*(t-rt)*(t-othr); t = rt+s: lc*s*(s+e)
            // z = m*s: s = c1*lc*e/(m^2 - c1*lc)
            tsub := rt + cs[1]*lc*e/(t^2 - cs[1]*lc);
        end if;
    end if;
    if not okp then continue; end if;
    val := cs[2]*Evaluate(LAs[2], tsub);      // in QT; need: square
    quart := SqClassFF(val);                  // squarefree representative, deg <= 4
    printf "   quartic in m (need square): %o\n", quart;
    if Degree(quart) le 2 then
        printf "   (degree <= 2 -- degenerate, note for manual follow-up)\n";
        continue;
    end if;
    C := HyperellipticCurve(quart);
    pts := Points(C : Bound := 3000);
    printf "   points (bound 3000): %o\n", #pts;
    for P in pts do
        if P[3] eq 0 then continue; end if;
        mv := P[1]/P[3];
        okv := true; tv := RQ!0;
        try tv := Evaluate(tsub, mv); catch e okv := false; end try;
        if not okv or tv in EXCL26 then continue; end if;
        // belt and braces: verify both u-free conditions at tv via the class reps
        v1 := Evaluate(LAs[1], tv); v2 := Evaluate(LAs[2], tv);
        if v1 eq 0 or v2 eq 0 then continue; end if;
        if not (IsSquare(cs[1]*v1) and IsSquare(cs[2]*v2)) then
            printf "   WARN t*=%o fails direct check (mv=%o)\n", tv, mv;
            continue;
        end if;
        Append(~tstars, <si, tv>);
        printf "   TSTAR sigma=%o t*=%o (m=%o)\n", SIGMAS[si], tv, mv;
    end for;
end for;
printf "TSTARS total %o\n", #tstars;
printf "DERIVE_DONE 2266\n";
quit;
