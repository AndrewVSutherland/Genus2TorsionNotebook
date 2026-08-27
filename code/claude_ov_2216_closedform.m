//////////////////////////////////////////////////////////////////////
//  claude_ov_2216_closedform.m   -- lane OV/2216
//
//  Turn the halving criterion on M_1(8,2,2) into FOUR explicit polynomial
//  square conditions in (u,v), with NO x-T module and no number field --
//  i.e. into something a C sieve or a PARI script can impose directly.
//
//  Setup.  delta'(Q+T) = (d1,d2,d3,dK), dK = a + b*theta in K,
//  theta^2 = B*theta - C, and lambda0 = squarefree(d1 d2 d3 N(dK)).
//  Halving <=> lambda0*dj square in Q (j=1,2,3) and xi = lambda0*dK
//  square in K.  Writing theta = (B + sqrt(DK))/2, DK = B^2 - 4C:
//
//      xi = alpha + beta*sqrt(DK),  alpha = lambda0*(a + b*B/2),
//                                    beta  = lambda0*b/2,
//      N(xi) = lambda0^2 * N(dK),   and  N(dK) = u*v*(u+v+1) * s^2
//
//  for an explicit s in Q(u,v) (the norm class is u*v*(u+v+1) for ALL
//  eight twists -- results/claude_ov_2216_s3check.log).  So on the locus
//  u*v*(u+v+1) = t^2 the square root of N(xi) is c = lambda0*s*t and
//
//      eps_+- = (alpha +- c)/2 ~ lambda0 * ( P +- 2*s*t ),  P := 2a + b*B,
//
//  modulo squares.  Halving <=> eps_+ or eps_- is a rational SQUARE.
//
//  This script prints lambda0, P and s for every condition set, and then
//  VALIDATES the closed form against code/claude_ov_2216_delta.m + Magma
//  IsDivisibleBy at every point of a candidate file.
//
//  Usage:
//    code/claude_magma_slot.sh -b candidate_file:=data/claude_ov_2216_all15_uv.txt \
//        code/claude_ov_2216_closedform.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(4*10^9);

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
    fa := Factorization(h);
    unit := h;
    for t in fa do unit := unit div t[1]^t[2]; end for;
    out := R!1;
    for t in fa do if IsOdd(t[2]) then out *:= t[1]; end if; end for;
    cst := Rationals()!unit;
    n := Numerator(cst); d := Denominator(cst);
    m := n*d; s := Sign(m); m := AbsoluteValue(m); rr := 1;
    for pe in Factorization(m) do if IsOdd(pe[2]) then rr *:= pe[1]; end if; end for;
    return (s*rr)*out;
end function;

// exact square root in Q(u,v) (assumes the argument IS a square)
function SqrtFF(gin)
    g := FF!gin;
    nu := R!Numerator(g); de := R!Denominator(g);
    function HalfPoly(h)
        fa := Factorization(h);
        unit := h;
        for t in fa do unit := unit div t[1]^t[2]; end for;
        out := R!1;
        for t in fa do
            assert IsEven(t[2]);
            out *:= t[1]^(t[2] div 2);
        end for;
        ok, sq := IsSquare(Rationals()!unit);
        assert ok;
        return sq*out;
    end function;
    return FF!HalfPoly(nu) / FF!HalfPoly(de);
end function;

function ModQ(g) return g mod qtm; end function;
function NormK(g)
    h := ModQ(g);
    a := Coefficient(h,0); b := Coefficient(h,1);
    return a^2 + a*b*B + b^2*C;
end function;
function DeltaCop(uD)
    return [* Evaluate(uD,r[1]), Evaluate(uD,r[2]), Evaluate(uD,r[3]), ModQ(uD) *];
end function;
function DeltaTwo(uS)
    cof := f div uS;  assert uS*cof eq f;
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
setlab := [ "SET1 = surface 0  (T = 0, B_quad)",
            "SET2             (T = A_r1r2, C_r3)",
            "SET3             (T = A_r1r3, C_r2)",
            "SET4 = surface 1  (T = A_r2r3, C_r1)" ];
uSrep  := [ PX!1, lin[1]*lin[2], lin[1]*lin[3], lin[2]*lin[3] ];

LAM := []; PP_ := []; SS := []; G12 := []; G13 := [];
print "######## closed-form halving criterion, per condition set ########";
print "  conditions:  lam0*d1*d2 = square,  lam0*d1*d3 = square  (printed as g12, g13),";
print "               u*v*(u+v+1) = t^2,  and  lam0*(P +- 2*s*t) = square.";
for j in [1..4] do
    dp := (j eq 1) select dQ else DMul(dQ, DeltaTwo(uSrep[j]));
    lam0 := SqFree(dp[1]*dp[2]*dp[3]*NormK(dp[4]));
    h := ModQ(dp[4]);
    a := Coefficient(h,0); b := Coefficient(h,1);
    Pv := 2*a + b*B;
    Nk := NormK(dp[4]);
    s := SqrtFF(Nk / (u*v*(u+v+1)));
    Append(~LAM, lam0); Append(~PP_, Pv); Append(~SS, s);
    Append(~G12, SqFree(dp[1]*dp[2])); Append(~G13, SqFree(dp[1]*dp[3]));
    printf "\n---- %o\n", setlab[j];
    printf "   g12   = %o\n", Factorization(SqFree(dp[1]*dp[2]));
    printf "   g13   = %o\n", Factorization(SqFree(dp[1]*dp[3]));
    printf "   lam0  = %o\n", Factorization(lam0);
    printf "   P     = %o\n", Pv;
    printf "   s     = %o\n", s;
end for;

// =====================================================================
print "";
print "######## numerical validation against the delta module + IsDivisibleBy ########";
load "code/claude_ov_2216_delta.m";
Qq := Rationals();
PZ<Y> := PP;

if not assigned candidate_file then
    candidate_file := "data/claude_ov_2216_all15_uv.txt";
end if;
function ParseQ(s)
    if "/" in s then
        pr := Split(s, "/");
        return Qq!StringToInteger(pr[1]) / Qq!StringToInteger(pr[2]);
    end if;
    return Qq!StringToInteger(s);
end function;

lines := Split(Read(candidate_file), "\n");
npt := 0; nchk := 0; nagree := 0; nonsurf := 0;
for raw in lines do
    s0 := raw;
    if #s0 gt 0 and s0[#s0] eq "\r" then s0 := s0[1..#s0-1]; end if;
    if #s0 eq 0 or s0[1] eq "#" then continue; end if;
    prt := [t : t in Split(s0, " ") | #t gt 0];
    if #prt lt 2 then continue; end if;
    uu := ParseQ(prt[1]); vv := ParseQ(prt[2]);
    npt +:= 1;

    Bn := uu^2*vv - uu^2 + uu*vv^2 - uu - vv^2 - vv - 2;
    Cn := uu^2 + uu*vv + vv^2 + uu + vv + 1;
    qn := Y^2 - Bn*Y + Cn;
    fn := ((1-uu)*Y + 1)*((1-vv)*Y + 1)*((uu+vv+2)*Y + 1)*(-qn);
    L := 1; for i in [0..5] do L := LCM(L, Denominator(Coefficient(fn,i))); end for;
    fI := PZ!(L^2*fn);
    rn := [ 1/(uu-1), 1/(vv-1), -1/(uu+vv+2) ];
    KK<th> := NumberField(qn);
    Cv := HyperellipticCurve(fI);
    J := Jacobian(Cv);
    DQ := J![Y + 1, Qq!(L*uu*vv*(uu+vv+1))];
    linn := [ Y - rn[i] : i in [1..3] ];
    JT := [ J!0, J![linn[1]*linn[2],0], J![linn[1]*linn[3],0], J![linn[2]*linn[3],0] ];

    okt, tt := IsSquare(uu*vv*(uu+vv+1));
    assert okt;                       // condition (3) holds at every candidate

    for j in [1..4] do
        lamn := Evaluate(LAM[j], [uu,vv]);
        Pn   := Evaluate(Numerator(FF!PP_[j]), [uu,vv]) / Evaluate(Denominator(FF!PP_[j]), [uu,vv]);
        sn   := Evaluate(Numerator(FF!SS[j]),  [uu,vv]) / Evaluate(Denominator(FF!SS[j]),  [uu,vv]);
        g12n := Evaluate(G12[j], [uu,vv]);
        g13n := Evaluate(G13[j], [uu,vv]);
        ratOK := IsSquare(g12n) and IsSquare(g13n);
        cf := IsSquare(lamn*(Pn + 2*sn*tt)) or IsSquare(lamn*(Pn - 2*sn*tt));
        closed := ratOK and cf;
        magdiv := IsDivisibleBy(DQ + JT[j], 2);
        nchk +:= 1;
        if closed eq magdiv then nagree +:= 1; end if;
        if ratOK then
            nonsurf +:= 1;
            printf "  u=%o v=%o set=%o ONSURFACE eps+=%o eps-=%o field=%o magma=%o\n",
                uu, vv, j, CO_SqClass(lamn*(Pn + 2*sn*tt)), CO_SqClass(lamn*(Pn - 2*sn*tt)),
                cf, magdiv;
        end if;
        assert closed eq magdiv;
    end for;
end for;
printf "\npoints=%o  (point,set) checks=%o  closedform_agrees_with_magma=%o  on-surface=%o\n",
       npt, nchk, nagree, nonsurf;
assert nagree eq nchk;
print "CLOSED_FORM_VALIDATED";
print "SEARCH_DONE";
quit;
