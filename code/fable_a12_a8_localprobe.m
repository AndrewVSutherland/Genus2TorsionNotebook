//////////////////////////////////////////////////////////////////////
// fable_a12_a8_localprobe.m  (2026-07-18, Fable session)
//
// Finite-field good-open viability probes on TWO charts:
//
//  (A) A(12) chart (p,z,r) [agent_a12_224_descent_setup.m formulas]:
//      count points whose J(F_q) CONTAINS [3,12].  The M(2,12)+3 lane
//      is dead mod 5 (halving wall); question: is the A(12) chart also
//      dead mod 5?  Canary: containment of [12] should be common.
//
//  (B) A(8) chart (r,p,t): count points whose J(F_q) contains [8,8].
//      Canary: containment of [8].
//
// Local logic: a global curve in the chart with the target torsion
// reduces, at every good prime q, to a chart point over F_q whose
// J(F_q) contains the target.  hits=0 across the full chart at a good
// prime is a chart-level local obstruction.
//
// Run:  magma -b code/fable_a12_a8_localprobe.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(2*10^9);

// A ascending-invariant containment test: target ta embeds in group
// with invariants iv (both divisor chains) iff aligned tails divide.
function ContainsInv(iv, ta)
    k := #iv; m := #ta;
    if k lt m then return false; end if;
    return &and[ IsDivisibleBy(iv[k-m+j], ta[j]) : j in [1..m] ];
end function;

// ---------- (A) A(12) chart over F_q ----------
for q in [5,7,11] do
    K := GF(q);
    Pol<x> := PolynomialRing(K);
    tested := 0; good := 0; c12 := 0; c312 := 0;
    ex := [];
    for pv0, zv0, rv0 in K do
        pv := pv0; zv := zv0; rv := rv0;
        if pv eq 0 or zv eq 0 then continue; end if;
        s := (zv^2 - 4*pv^2 + 1)/(2*zv);
        if s^2 eq 1 then continue; end if;
        t := (zv^2 + 4*pv^2 - 1)^2/(8*pv^2*zv);
        mu := ((s^2 - 1)*(2*pv*rv + 1) - pv^2*(2*s*t - 4))/(4*pv^3);
        lam := (4 - mu^2)*pv^2/(s^2 - 1);
        if lam eq 0 then continue; end if;
        tested +:= 1;
        T1 := pv*x + rv;
        R := (T1^2 + x - 1)/lam;
        F := R*x^2 + 4*(R + x - 1)*(R - 1);
        f := R*F;
        if Degree(f) ne 6 or not IsSeparable(f) then continue; end if;
        okc := true;
        try
            C := HyperellipticCurve(f);
            G := AbelianGroup(Jacobian(C));
        catch ee okc := false; end try;
        if not okc then continue; end if;
        good +:= 1;
        iv := Invariants(G);
        if ContainsInv(iv, [12]) then c12 +:= 1; end if;
        if ContainsInv(iv, [3,12]) then
            c312 +:= 1;
            if #ex lt 4 then Append(~ex, <pv,zv,rv,iv>); end if;
        end if;
    end for;
    printf "A12PLUS3 q=%o tested=%o good=%o canary12=%o contains312=%o\n", q, tested, good, c12, c312;
    for e in ex do printf "A12HIT q=%o p=%o z=%o r=%o inv=%o\n", q, e[1], e[2], e[3], e[4]; end for;
end for;

// ---------- (B) A(8) chart over F_q ----------
for q in [7,11,13] do
    K := GF(q);
    Pol<x> := PolynomialRing(K);
    tested := 0; good := 0; c8 := 0; c88 := 0;
    ex := [];
    for rv0, pv0, tv0 in K do
        rv := rv0; pv := pv0; tv := tv0;
        if rv eq 0 or tv eq 0 then continue; end if;
        e := tv^2 - 2*pv*tv/rv;
        d := e + 2*pv - rv^2;
        lam := rv/tv;
        aa := rv^2 - lam;
        bb := 2*rv*pv - 2*lam*(pv + rv*tv) + 2*rv*lam;
        cc := pv^2 + 2*pv*rv^2 - rv^4 - rv^3*tv - rv*pv^2/tv
              - lam*(rv^2 + e) + 2*lam*(rv*pv + rv^2*tv - 3*pv*tv + rv*tv^2);
        tested +:= 1;
        Qq := x^2 + d;
        qq := aa*x^2 + bb*x + cc;
        f := qq*(Qq^2 + qq);
        if Degree(f) ne 6 or not IsSeparable(f) then continue; end if;
        okc := true;
        try
            C := HyperellipticCurve(f);
            G := AbelianGroup(Jacobian(C));
        catch ee okc := false; end try;
        if not okc then continue; end if;
        good +:= 1;
        iv := Invariants(G);
        if ContainsInv(iv, [8]) then c8 +:= 1; end if;
        if ContainsInv(iv, [8,8]) then
            c88 +:= 1;
            if #ex lt 4 then Append(~ex, <rv,pv,tv,iv>); end if;
        end if;
    end for;
    printf "A8X88 q=%o tested=%o good=%o canary8=%o contains88=%o\n", q, tested, good, c8, c88;
    for e in ex do printf "A8HIT q=%o r=%o p=%o t=%o inv=%o\n", q, e[1], e[2], e[3], e[4]; end for;
end for;

printf "PROBE_DONE\n";
quit;
