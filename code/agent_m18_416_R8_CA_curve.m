
//////////////////////////////////////////////////////////////////////
//  R = -8: the CORRECT S_A-condition curve C_A and its structure.
//
//  C_A : y^4 - 4*alpha_A(X(l))*y^2 + d_A(X(l)) = 0  over the lambda-line.
//  Rational points of C_A with y != 0 and bb != 0 biject with rational
//  (l, m, y): m^2 = X(l), y^2 = V_A = 2*alpha_A + m*bb
//  via m = (y^2 - 2*alpha_A)/bb.  So C_A(Q) finite+listed would reduce
//  the S_A condition on the fiber R=-8 to finitely many lambda.
//
//  z = y^2 maps C_A -> E_m (rank 1), so Jac(C_A) ~ E_m x Prym.
//  This script:
//   1. builds C_A as an exact plane curve, computes its genus;
//   2. searches small rational points;
//   3. attempts to split the Prym: quotients of C_A by visible
//      involutions (y -> -y composed with the conjugations of the
//      lambda-line), producing candidate lower-genus quotients whose
//      Jacobian pieces can be rank-bounded;
//   4. prints exact equations for external Chabauty work.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

Rv := Q!-8;
Kv := Q!1008;

// symbolic data in X
QX<Xv> := RationalFunctionField(Q);
wX := (Xv + Kv)/(Xv - Kv);
tX := (2*Rv^2 + (1-wX^2)*Rv - 2*wX^2)/(4*(wX^2-1));
c4X := Rv + 2 + 4*tX;
PXx<xx> := PolynomialRing(QX);
AX := xx^2 + (Rv^3 + 4*Rv^2*tX + Rv - 8*Rv*tX + 4*tX)*xx + Rv^4;
XRX := -c4X*Rv;
AtX := PXx![QX!co : co in Coefficients(c4X^2*Evaluate(AX, xx/c4X))];
GX := 2*(Rv^2-1)*(Rv*(2*Rv+1) - wX^2*(Rv+2));
alphaAX := XRX + Coefficient(AtX,1)/2;
dAX := Discriminant(AtX);
okA, hAmX := IsSquare(Xv*Evaluate(AtX, XRX)/GX); assert okA;

// lambda-line data
Ql<l> := RationalFunctionField(Q);
Xl := (-2016*l - 1824)/(l^2 - 1);
Sl := (-1008*l^2 - 1824*l - 1008)/(l^2 - 1);
assert (126*Sl)^2 eq Evaluate(Numerator(GX*(Xv-Kv)^2), Xl)/Evaluate(Denominator(GX*(Xv-Kv)^2), Xl);

aL := 2*Evaluate(alphaAX, Xl);          // 2*alpha_A
dL := Evaluate(dAX, Xl);                // d_A
gl := 126*Sl/(Xl - Kv);
bbL := 2*gl*Evaluate(hAmX, Xl)/Xl;      // bb, with V_A = aL/1... V_A = aL + m*bbL? careful:
// V_A = 2*alpha_A + 2*h_A*g = aL + m*bbL/... define exactly:
// h_A = hAm/m, so 2*h_A*g = 2*g*hAm/m = m*(2*g*hAm/X). So bb = 2*g*hAm/X. OK:
assert aL^2 - Xl*bbL^2 eq dL;   // wait: V_A*V_A' = (aL + m bb)(aL - m bb) = aL^2 - X bb^2
print "identity (2 alpha_A)^2 - X*bb^2 = d_A verified on lambda-line";

// ---- C_A as a plane curve ----
// y^4 - 2*aL*y^2 + dL = 0   (note V_A + V_A' = 2*aL... with aL = 2 alpha:
// roots of z^2 - 2*aL*z + (aL^2 - X bb^2)?? No: V_A + V_A' = 2*aL, product aL^2 - X bb^2 = dL.
// So z^2 - 2*aL*z + dL = 0 with z = y^2:
//   C_A : y^4 - 2*aL(l)*y^2 + dL(l) = 0.
A2<lu, yu> := AffineSpace(Q, 2);
num_a := Numerator(aL); den_a := Denominator(aL);
num_d := Numerator(dL); den_d := Denominator(dL);
// common denominator D = lcm(den_a, den_d); multiply through by D:
Dl := LCM(den_a, den_d);
ca := Dl div den_a; cd := Dl div den_d;
// C_A cleared: D*y^4 - 2*num_a*ca*y^2 + num_d*cd = 0
polyCA := Evaluate(Dl, lu)*yu^4 - 2*Evaluate(num_a*ca, lu)*yu^2 + Evaluate(num_d*cd, lu);
printf "C_A cleared plane model:\n%o = 0\n\n", polyCA;
CA := Curve(A2, polyCA);
try
    gCA := Genus(ProjectiveClosure(CA));
    printf "genus(C_A) = %o\n", gCA;
catch e
    printf "genus computation failed: %o\n", e`Object;
end try;

// ---- small points ----
print "small rational points (lambda height <= 60, y from quadratic):";
Pz<zz> := PolynomialRing(Q);
found := [];
for dn in [1..60] do for nu in [-60..60] do
    if GCD(nu,dn) ne 1 then continue; end if;
    lv := Q!nu/dn;
    if lv in {Q!1, Q!-1} then continue; end if;
    if Evaluate(den_a, lv) eq 0 or Evaluate(den_d, lv) eq 0 then continue; end if;
    av := Evaluate(num_a, lv)/Evaluate(den_a, lv);
    dv := Evaluate(num_d, lv)/Evaluate(den_d, lv);
    // z^2 - 2 av z + dv = 0
    disc := av^2 - dv;
    if disc lt 0 then continue; end if;
    oks, sd := IsSquare(disc);
    if not oks then continue; end if;
    for z in [av + sd, av - sd] do
        if z le 0 then continue; end if;
        oky, yv := IsSquare(z);
        if oky then
            Append(~found, <lv, yv>);
            printf "  (l, y) = (%o, %o)   [z=%o]\n", lv, yv, z;
        end if;
    end for;
end for; end for;
printf "total small points found: %o\n", #found;

// ---- quotient structure ----
// z = y^2: quotient is E_m (rank 1).  The Prym: quotient of C_A by the
// composite involution iota: (l, y) -> (l', y') where l' = conj(l)
// (deck of X: same X, S -> -S) -- under which aL is invariant? aL is a
// function of X only => invariant; dL too.  bb -> ? S -> -S flips bb.
// V_A(l', m) = aL - m bb: so (l,y) -> (l', y) maps C_A to itself
// (the defining quartic in y has X-only coefficients!).  Quotients by
// tau: (l,y) -> (l', y) and tau': (l,y) -> (l', -y):
// Note the quartic C_A depends only on X(l)!  So C_A = pullback of
//    C_A^X : y^4 - 2*aX(X)*y^2 + dX(X) = 0
// along the double cover lambda-line -> X-line.  Rational points of C_A
// map to rational points of C_A^X: analyze C_A^X first!
num_aX := Numerator(2*alphaAX*2); den_aX := Denominator(2*alphaAX*2);
// aX = 2*(2 alphaA)?? careful: aL = 2*alpha_A evaluated at Xl.  In X-coords:
aXfun := 2*alphaAX;   // = aL as function of X
dXfun := dAX;
DX := LCM(Denominator(aXfun), Denominator(dXfun));
caX := DX div Denominator(aXfun); cdX := DX div Denominator(dXfun);
A2X<xu, yu2> := AffineSpace(Q, 2);
polyCAX := Evaluate(DX, xu)*yu2^4 - 2*Evaluate(Numerator(aXfun)*caX, xu)*yu2^2
           + Evaluate(Numerator(dXfun)*cdX, xu);
printf "\nC_A^X (X-line model, C_A is its pullback along X = X(l)):\n%o = 0\n\n", polyCAX;
CAX := Curve(A2X, polyCAX);
try
    gX := Genus(ProjectiveClosure(CAX));
    printf "genus(C_A^X) = %o\n", gX;
catch e
    printf "genus(C_A^X) failed: %o\n", e`Object;
end try;
// z = y^2 quotient of C_A^X: z^2 - 2 aX z + dX = 0: conic-like in (X,z):
// (z - aX)^2 = aX^2 - dX = X*bbX^2-ish: u^2 = X * (stuff)^2?? check:
delta := aXfun^2 - dXfun;
printf "aX^2 - dX factorization (should be X * square):\n num %o\n den %o\n",
    Factorization(Numerator(delta)), Factorization(Denominator(delta));

// ---- B-side X-line curve C_B^X ----
BX := c4X*xx^2 + (Rv^2 + 4*Rv + 1 + 8*tX)*xx + (2*Rv^2 + Rv + 4*tX);
BtX := PXx![QX!co : co in Coefficients(c4X*Evaluate(BX, xx/c4X))];
alphaBX := XRX + Coefficient(BtX,1)/2;
dBX := Discriminant(BtX);
aBfun := 2*alphaBX;
DXB := LCM(Denominator(aBfun), Denominator(dBX));
caB := DXB div Denominator(aBfun); cdB := DXB div Denominator(dBX);
polyCBX := Evaluate(DXB, xu)*yu2^4 - 2*Evaluate(Numerator(aBfun)*caB, xu)*yu2^2
           + Evaluate(Numerator(dBX)*cdB, xu);
printf "\nC_B^X (X-line model):\n%o = 0\n\n", polyCBX;
CBX := Curve(A2X, polyCBX);
try
    gBX := Genus(ProjectiveClosure(CBX));
    printf "genus(C_B^X) = %o\n", gBX;
catch e
    printf "genus(C_B^X) failed: %o\n", e`Object;
end try;
print "DONE";
quit;
