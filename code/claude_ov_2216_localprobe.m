//////////////////////////////////////////////////////////////////////
//  claude_ov_2216_localprobe.m   -- lane OV/2216
//
//  Is the FULL halving condition (norm surface + the eps/field condition)
//  even soluble mod p?  If some prime p admits NO nondegenerate solution
//  over F_p, that is a candidate everywhere-local obstruction and would
//  close chart 1 outright, rather than by box exhaustion.
//
//  Uses the closed form validated in results/claude_ov_2216_closedform.log:
//      g12 = square, g13 = square, u*v*(u+v+1) = t^2,
//      lam0*(P + 2*s*t) = square  OR  lam0*(P - 2*s*t) = square,
//      lam0 = (u-1)(v-1)(u+v+2)  (the same for all four condition sets).
//
//  For each prime and each condition set we count nondegenerate (u,v) in
//  F_p^2 satisfying:
//      A) the three rational conditions only            [the norm surface]
//      B) all four                                      [the true locus]
//  A prime with countB = 0 but countA > 0 is the interesting case.
//
//  CAVEAT (stated in the log too): an F_p count is a probe, not a proof.
//  It tests the unit stratum only; a genuine Q_p-insolubility claim needs
//  the non-unit strata and Hensel bookkeeping as well.
//
//  Usage:
//    code/claude_magma_slot.sh -b PMAX:=150 code/claude_ov_2216_localprobe.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(4*10^9);

if not assigned PMAX then PMAX := 150;
elif Type(PMAX) eq MonStgElt then PMAX := StringToInteger(PMAX); end if;

FF<u,v> := RationalFunctionField(Rationals(), 2);
R := PolynomialRing(Rationals(), 2);
PX<X> := PolynomialRing(FF);

B := u^2*v - u^2 + u*v^2 - u - v^2 - v - 2;
C := u^2 + u*v + v^2 + u + v + 1;
qtm := X^2 - B*X + C;
f   := ((1-u)*X + 1)*((1-v)*X + 1)*((u+v+2)*X + 1)*(-qtm);
r   := [ 1/(u-1), 1/(v-1), -1/(u+v+2) ];
lin := [ X - r[i] : i in [1..3] ];

function SqFree(gin)
    g := FF!gin;
    if g eq 0 then return R!0; end if;
    h := R!Numerator(g) * R!Denominator(g);
    fa := Factorization(h);  unit := h;
    for t in fa do unit := unit div t[1]^t[2]; end for;
    out := R!1;
    for t in fa do if IsOdd(t[2]) then out *:= t[1]; end if; end for;
    cst := Rationals()!unit;
    n := Numerator(cst); d := Denominator(cst);
    m := n*d; s := Sign(m); m := AbsoluteValue(m); rr := 1;
    for pe in Factorization(m) do if IsOdd(pe[2]) then rr *:= pe[1]; end if; end for;
    return (s*rr)*out;
end function;
function SqrtFF(gin)
    g := FF!gin;
    function HalfPoly(h)
        fa := Factorization(h); unit := h;
        for t in fa do unit := unit div t[1]^t[2]; end for;
        out := R!1;
        for t in fa do assert IsEven(t[2]); out *:= t[1]^(t[2] div 2); end for;
        ok, sq := IsSquare(Rationals()!unit); assert ok;
        return sq*out;
    end function;
    return FF!HalfPoly(R!Numerator(g)) / FF!HalfPoly(R!Denominator(g));
end function;
function ModQ(g) return g mod qtm; end function;
function NormK(g)
    h := ModQ(g); a := Coefficient(h,0); b := Coefficient(h,1);
    return a^2 + a*b*B + b^2*C;
end function;
function DeltaCop(uD)
    return [* Evaluate(uD,r[1]), Evaluate(uD,r[2]), Evaluate(uD,r[3]), ModQ(uD) *];
end function;
function DeltaTwo(uS)
    cof := f div uS; assert uS*cof eq f;
    out := [* *];
    for i in [1..3] do
        if Evaluate(uS,r[i]) eq 0 then Append(~out, -Evaluate(cof,r[i]));
        else Append(~out, Evaluate(uS,r[i])); end if;
    end for;
    vS := ModQ(uS);
    if vS eq 0 then Append(~out, -ModQ(cof)); else Append(~out, vS); end if;
    return out;
end function;
function DMul(a,b) return [* a[1]*b[1], a[2]*b[2], a[3]*b[3], ModQ(a[4]*b[4]) *]; end function;

dQ := DeltaCop(X+1);
uSrep  := [ PX!1, lin[1]*lin[2], lin[1]*lin[3], lin[2]*lin[3] ];
setlab := [ "SET1(surface0)", "SET2", "SET3", "SET4(surface1)" ];

// numerators/denominators as bivariate polynomials, for fast mod-p evaluation
NU := [];  DE := [];
for j in [1..4] do
    dp := (j eq 1) select dQ else DMul(dQ, DeltaTwo(uSrep[j]));
    lam0 := FF!SqFree(dp[1]*dp[2]*dp[3]*NormK(dp[4]));
    h := ModQ(dp[4]);
    Pv := 2*Coefficient(h,0) + Coefficient(h,1)*B;
    s  := SqrtFF(NormK(dp[4]) / (u*v*(u+v+1)));
    g12 := FF!SqFree(dp[1]*dp[2]);
    g13 := FF!SqFree(dp[1]*dp[3]);
    Append(~NU, [ R!Numerator(g12), R!Numerator(g13), R!Numerator(lam0),
                  R!Numerator(Pv), R!Numerator(s) ]);
    Append(~DE, [ R!Denominator(g12), R!Denominator(g13), R!Denominator(lam0),
                  R!Denominator(Pv), R!Denominator(s) ]);
end for;

printf "claude_ov_2216_localprobe PMAX=%o\n", PMAX;
print "# countA = nondegenerate (u,v) in F_p^2 on the norm surface (3 rational conditions)";
print "# countB = those that ALSO satisfy the eps (quadratic-field) condition";
print "# CAVEAT: an F_p count probes the unit stratum only; it is evidence, not a proof.";

zeroB := [];
p := 5;
while p le PMAX do
    Fp := GF(p);
    RP := PolynomialRing(Fp, 2);
    // quadratic-residue table
    QR := [ false : i in [1..p] ];
    for x in [1..p-1] do QR[(x*x mod p) + 1] := true; end for;
    SQRT := [ -1 : i in [1..p] ];
    for x in [0..p-1] do SQRT[(x*x mod p) + 1] := x; end for;

    for j in [1..4] do
        nu := [ RP!g : g in NU[j] ];
        de := [ RP!g : g in DE[j] ];
        cA := 0; cB := 0;
        for ui in [0..p-1] do for vi in [0..p-1] do
            uu := Fp!ui; vv := Fp!vi;
            ww := -(uu+vv+1);
            // nondegeneracy: {1,u,v,w} distinct and none zero
            if uu eq 0 or vv eq 0 or ww eq 0 then continue; end if;
            if uu eq 1 or vv eq 1 or ww eq 1 then continue; end if;
            if uu eq vv or uu eq ww or vv eq ww then continue; end if;
            bad := false;
            for k in [1..5] do
                if Evaluate(de[k], [uu,vv]) eq 0 then bad := true; break; end if;
            end for;
            if bad then continue; end if;
            val := [ Evaluate(nu[k],[uu,vv]) / Evaluate(de[k],[uu,vv]) : k in [1..5] ];
            g12 := val[1]; g13 := val[2]; lam0 := val[3]; Pv := val[4]; sv := val[5];
            c3 := uu*vv*(uu+vv+1);
            if c3 eq 0 or g12 eq 0 or g13 eq 0 then continue; end if;
            if not QR[Integers()!c3 + 1] then continue; end if;
            if not QR[Integers()!g12 + 1] then continue; end if;
            if not QR[Integers()!g13 + 1] then continue; end if;
            cA +:= 1;
            t := Fp!SQRT[Integers()!c3 + 1];
            e1 := lam0*(Pv + 2*sv*t);
            e2 := lam0*(Pv - 2*sv*t);
            ok1 := (e1 eq 0) or QR[Integers()!e1 + 1];
            ok2 := (e2 eq 0) or QR[Integers()!e2 + 1];
            if ok1 or ok2 then cB +:= 1; end if;
        end for; end for;
        printf "p=%-4o %-15o countA=%-7o countB=%-7o ratio=%o\n",
               p, setlab[j], cA, cB, cA gt 0 select RealField(4)!(cB/cA) else 0;
        if cA gt 0 and cB eq 0 then Append(~zeroB, <p, j>); end if;
    end for;
    p := NextPrime(p);
end while;

printf "\nPRIMES_WITH_NO_FULL_SOLUTION %o\n", zeroB;
printf "LOCAL_OBSTRUCTION_CANDIDATE %o\n", #zeroB gt 0;
print "SEARCH_DONE";
quit;
