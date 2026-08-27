
//////////////////////////////////////////////////////////////////////
//  R = -8 candidate certificate via the rank-0 V_B quotient: VALIDATION.
//
//  agent_m18_416_R8_SA_product_quotients.m found the quotient
//     y^2 = 3969*l^3 + 3654*l^2 - 4095*l - 3528        (*)
//  (squarefree model of y_B^2 = V_B(l)) to have rank bounds [0,0] and
//  torsion [2,2].  If V_B is the correct S_B-condition function on the
//  lambda-line (covering BOTH g-signs via the conic involution l -> l'),
//  then every S_B-passing rational point of the fiber R=-8 maps to a
//  rational point of (*), so S_B can hold only at finitely many lambda.
//
//  NOTE: agent_m18_416_R8_mwsieve_attempt.md has a DIFFERENT (rank-1)
//  E_Y for the S_B quotient; one normalization must be wrong.  This
//  script validates V_B against ground truth:
//   1. exact S_B test (number-field square) at sample w(lambda) vs the
//      prediction [V_B(l) or V_B(l') is a square];
//   2. enumerates ALL rational points of the rank-0 curve (*) via the
//      elliptic model, maps them to lambda;
//   3. runs the full exact second-stage test at each such lambda.
//  If 1 validates and 3 all fail ==> R=-8 fiber is EMPTY (certificate).
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

Rv := Q!-8;
Kv := Q!1008;
P<xq> := PolynomialRing(Q);

// lambda data
function XofL(l) return (-2016*l - 1824)/(l^2 - 1); end function;
function SofL(l) return (-1008*l^2 - 1824*l - 1008)/(l^2 - 1); end function;
// conjugate lambda': other solution of X(T) = X(l)
function ConjL(l)
    // (-2016T - 1824)(l^2-1) = (-2016 l - 1824)(T^2 - 1)
    PT<T> := PolynomialRing(Q);
    pol := (-2016*T - 1824)*(l^2-1) - (-2016*l - 1824)*(T^2 - 1);
    rts := [r[1] : r in Roots(pol) | r[1] ne l];
    if #rts eq 0 then return l; end if;   // fixed point
    return rts[1];
end function;

// family data at w
function FamilyData(wv)
    tv := (2*Rv^2 + (1-wv^2)*Rv - 2*wv^2)/(4*(wv^2-1));
    c4v := Rv + 2 + 4*tv;
    Av := xq^2 + (Rv^3 + 4*Rv^2*tv + Rv - 8*Rv*tv + 4*tv)*xq + Rv^4;
    Bv := c4v*xq^2 + (Rv^2 + 4*Rv + 1 + 8*tv)*xq + (2*Rv^2 + Rv + 4*tv);
    return Av, Bv, c4v;
end function;

function IsSquareQ(qv)
    qv := Q!qv;
    if qv le 0 then return false; end if;
    return IsSquare(Numerator(qv)) and IsSquare(Denominator(qv));
end function;

// exact S_B test at w (component square test in Q[T]/Bt)
function ExactSB(wv)
    Av, Bv, c4v := FamilyData(wv);
    XR := -c4v*Rv;
    Btq := P![Q!co : co in Coefficients(c4v*Evaluate(Bv, xq/c4v))];
    dsc := Discriminant(Btq);
    if dsc eq 0 then return false, "degenerate"; end if;
    if IsSquareQ(dsc) then
        for rt in Roots(Btq) do
            val := XR - rt[1];
            if val eq 0 then return false, "zero"; end if;
            if not IsSquareQ(val) then return false, "split-fail"; end if;
        end for;
        return true, "split-pass";
    end if;
    Kf<th> := NumberField(Btq);
    return IsSquare(XR - th), "quad";
end function;

// exact S_A test at w
function ExactSA(wv)
    Av, Bv, c4v := FamilyData(wv);
    XR := -c4v*Rv;
    Atq := P![Q!co : co in Coefficients(c4v^2*Evaluate(Av, xq/c4v))];
    dsc := Discriminant(Atq);
    if dsc eq 0 then return false, "degenerate-discA=0"; end if;
    if IsSquareQ(dsc) then
        for rt in Roots(Atq) do
            val := XR - rt[1];
            if val eq 0 then return false, "zero"; end if;
            if not IsSquareQ(val) then return false, "split-fail"; end if;
        end for;
        return true, "split-pass";
    end if;
    Kf<th> := NumberField(Atq);
    return IsSquare(XR - th), "quad";
end function;

// V_B(l) exactly (rebuild from symbolic pieces at R=-8)
QX<Xv> := RationalFunctionField(Q);
wX := (Xv + Kv)/(Xv - Kv);
tX := (2*Rv^2 + (1-wX^2)*Rv - 2*wX^2)/(4*(wX^2-1));
c4X := Rv + 2 + 4*tX;
PXx<xx> := PolynomialRing(QX);
BX := c4X*xx^2 + (Rv^2 + 4*Rv + 1 + 8*tX)*xx + (2*Rv^2 + Rv + 4*tX);
XRX := -c4X*Rv;
BtX := PXx![QX!co : co in Coefficients(c4X*Evaluate(BX, xx/c4X))];
GX := 2*(Rv^2-1)*(Rv*(2*Rv+1) - wX^2*(Rv+2));
alphaBX := XRX + Coefficient(BtX,1)/2;
okB, hBX := IsSquare(Evaluate(BtX, XRX)/GX); assert okB;

function VBofL(l)
    Xl := XofL(l); Sl := SofL(l);
    g := 126*Sl/(Xl - Kv);
    return 2*Evaluate(alphaBX, Xl) + 2*Evaluate(hBX, Xl)*g;
end function;

// ---- 1. validation on sample lambdas ----
print "== validation: exact S_B vs V_B-quotient prediction ==";
nchk := 0; nagree := 0;
for lnum in [-30..30] do
    for lden in [1,2,3,7,9] do
        l := Q!lnum/lden;
        if l in {Q!1, Q!-1} then continue; end if;
        Xl := XofL(l);
        if Xl eq 0 or Xl eq Kv then continue; end if;
        wv := (Xl + Kv)/(Xl - Kv);
        if wv in {Q!0,Q!1,Q!-1} then continue; end if;
        lp := ConjL(l);
        okB1, _ := ExactSB(wv);
        v1 := VBofL(l);
        pred := false;
        if v1 ne 0 and IsSquareQ(v1) then pred := true; end if;
        if not pred and lp ne l then
            v2 := VBofL(lp);
            if v2 ne 0 and IsSquareQ(v2) then pred := true; end if;
        end if;
        nchk +:= 1;
        if okB1 eq pred then nagree +:= 1;
        else
            printf "MISMATCH l=%o w=%o exact=%o pred=%o\n", l, wv, okB1, pred;
        end if;
    end for;
end for;
printf "validation: %o/%o agree\n", nagree, nchk;

// ---- 2. all rational points of the rank-0 quotient ----
print "== rational points of y^2 = 3969 l^3 + 3654 l^2 - 4095 l - 3528 ==";
F := 3969*xq^3 + 3654*xq^2 - 4095*xq - 3528;
C := HyperellipticCurve(F);
pts := Points(C : Bound := 100000);
printf "small points: %o\n", [<pt[1],pt[2],pt[3]> : pt in pts];
// rank 0 + torsion [2,2]: rational points = preimages of torsion; the
// cubic has roots exactly at 2-torsion x-coords.
rts := Roots(F);
printf "cubic roots (2-torsion lambdas, y=0): %o\n", [r[1] : r in rts];

// ---- 3. full second-stage test at each candidate lambda ----
print "== full test at each rank-0-allowed lambda ==";
cands := [r[1] : r in rts] cat [pt[1]/pt[3] : pt in pts | pt[3] ne 0];
cands := Sort(Setseq(Seqset(cands)));
for l in cands do
    if l in {Q!1,Q!-1} then printf "l=%o: boundary (l^2=1)\n", l; continue; end if;
    Xl := XofL(l);
    if Xl eq 0 or Xl eq Kv then printf "l=%o: X boundary\n", l; continue; end if;
    wv := (Xl + Kv)/(Xl - Kv);
    okA, whyA := ExactSA(wv);
    okB, whyB := ExactSB(wv);
    okX, mrt := IsSquareQ(Xl);
    printf "l=%o : X=%o m-rational=%o | S_A=%o (%o) | S_B=%o (%o)%o\n",
        l, Xl, okX, okA, whyA, okB, whyB,
        (okA and okB and okX) select "   <<<< FULL PASS ?!" else "";
end for;
print "If no FULL PASS above and validation agrees, the R=-8 fiber is EMPTY.";
print "DONE";
quit;
