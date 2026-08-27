
//////////////////////////////////////////////////////////////////////
//  R = -8: quotient analysis of the S_A cover and the full S_A & S_B
//  fiber product, hunting a rank-zero killer quotient.
//
//  Background: the corrected S_B V_4 method (agent_m18_416_SB_v4_scan_els.m)
//  killed R=-25/4 and R=-29/8 but not R=-8 (all three S_B quotients have
//  rank 1).  The S_A side and the mixed characters of the fiber product
//  were not yet examined.  Any rank-0 elliptic quotient (or provably
//  point-finite higher-genus quotient) forces the lambda values into a
//  finite checkable set.
//
//  Objects (R=-8, K=1008; lambda-parametrization of E_R from
//  agent_m18_416_R8_mwsieve_attempt.md):
//    X(l)  = m^2 = (-2016*l - 1824)/(l^2 - 1)
//    S(l)  =        (-1008*l^2 - 1824*l - 1008)/(l^2 - 1),  S^2 = Gamma(X)
//    g = S/(X - K)  (so g^2 = G);
//    V_A(l) = 2*alpha_A(X) + 2*h_A(X)*S/(X-K)     (S_A condition: V_A = sq)
//    V_B(l) = 2*alpha_B(X) + 2*h_B(X)*S/(X-K)
//  Full halving needs m, y_A, y_B all rational:
//    m^2 = X(l),  y_A^2 = V_A(l),  y_B^2 = V_B(l)
//  => quotient curves for each nontrivial character:
//    u^2 = X, y_A^2 = V_A, y_B^2 = V_B, (y_A y_B)^2 = V_A V_B,
//    (m y_A)^2 = X V_A, (m y_B)^2 = X V_B, (m y_A y_B)^2 = X V_A V_B.
//  This script builds each as y^2 = squarefree poly(l), reports genus,
//  and for genus-1 quotients with a visible point computes rank bounds.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

Rv := Q!-8;
Kv := -2*Rv*(Rv^2-1);   // 1008
printf "R=%o K=%o\n", Rv, Kv;

// symbolic fiber data specialized at R=-8, as functions of X=m^2
QX<Xv> := RationalFunctionField(Q);
// w = (X+K)/(X-K); reconstruct alpha_A, h_A, alpha_B, h_B, G in X
wv := (Xv + Kv)/(Xv - Kv);
tv := (2*Rv^2 + (1-wv^2)*Rv - 2*wv^2)/(4*(wv^2-1));
c4v := Rv + 2 + 4*tv;
PXx<xx> := PolynomialRing(QX);
Av := xx^2 + (Rv^3 + 4*Rv^2*tv + Rv - 8*Rv*tv + 4*tv)*xx + Rv^4;
Bv := c4v*xx^2 + (Rv^2 + 4*Rv + 1 + 8*tv)*xx + (2*Rv^2 + Rv + 4*tv);
XRv := -c4v*Rv;
Atv := PXx![QX!co : co in Coefficients(c4v^2*Evaluate(Av, xx/c4v))];
Btv := PXx![QX!co : co in Coefficients(c4v*Evaluate(Bv, xx/c4v))];
Gv := 2*(Rv^2-1)*(Rv*(2*Rv+1) - wv^2*(Rv+2));
alphaA := XRv + Coefficient(Atv,1)/2;
alphaB := XRv + Coefficient(Btv,1)/2;
// h_A is ODD in m (h_A = hAm/m with hAm even), h_B is EVEN:
okA, hAm := IsSquare(Xv*Evaluate(Atv, XRv)/Gv); assert okA;   // hAm = m*h_A
okB, hB := IsSquare(Evaluate(Btv, XRv)/Gv); assert okB;

// Gamma(X) = G*(X-K)^2 should be a polynomial (quadratic in X)
Gam := Gv*(Xv - Kv)^2;
GamNum := Numerator(Gam); GamDen := Denominator(Gam);
assert Degree(GamDen) eq 0;
GamPoly := GamNum/LeadingCoefficient(GamDen);
printf "Gamma(X) = %o\n", GamPoly;

// lambda parametrization
Ql<l> := RationalFunctionField(Q);
Xl := (-2016*l - 1824)/(l^2 - 1);
Sl := (-1008*l^2 - 1824*l - 1008)/(l^2 - 1);
// verify S^2 = Gamma(X(l))
lhs := (126*Sl)^2;   // reduced coordinate: Gamma = (126*S)^2, 15876 = 126^2
rhs := Evaluate(QX!GamPoly/1, Xl);
assert lhs eq rhs;
print "lambda parametrization verified: (126*S(l))^2 = Gamma(X(l))";

// build the lambda-rational data.
// V_B(l) = 2*alpha_B + 2*h_B*S/(X-K) is even in m => lambda-rational.
// V_A = a + m*bb with a = 2*alpha_A, bb = 2*S*hAm/((X-K)*X): NOT
// lambda-rational, but its m-conjugate norm is:
//   (y_A * y_A')^2 = a^2 - X*bb^2 = d_A(X)   (exact identity).
function SubstX(fX, Xl)
    return Evaluate(fX, Xl);
end function;
gl := 126*Sl/(Xl - Kv);   // g = sqrt(Gamma)/(X-K)
VB := 2*SubstX(alphaB, Xl) + 2*SubstX(hB, Xl)*gl;
al := 2*SubstX(alphaA, Xl);
bbl := 2*gl*SubstX(hAm, Xl)/Xl;   // bb = 2*g*hAm/X (g includes the 126)
dAl := al^2 - Xl*bbl^2;
// verify against the direct discriminant
dAdirect := SubstX(Discriminant(Atv), Xl);
assert dAl eq dAdirect;
print "A-norm identity verified: a^2 - X*bb^2 = d_A(X(l))";
printf "V_B(l) = %o\n\n", VB;
printf "d_A(l) = %o\n\n", dAl;

// characters (all necessary conditions on the lambda-line)
chars := [
    <"X", Xl>,
    <"V_B", VB>,
    <"X*V_B", Xl*VB>,
    <"d_A", dAl>,
    <"X*d_A", Xl*dAl>,
    <"d_A*V_B", dAl*VB>,
    <"X*d_A*V_B", Xl*dAl*VB>
];

Pl<T> := PolynomialRing(Q);

function SquarefreeModel(fl)   // y^2 = fl(l): return squarefree poly F(T)
    num := Numerator(fl); den := Denominator(fl);
    F := Pl!num * Pl!den;      // y^2 f_den^2 = num*den
    // strip square factors, KEEPING the squarefree part of the leading
    // unit (Factorization returns monic factors x unit -- dropping the
    // unit changes the curve!)
    lcF := LeadingCoefficient(F);
    sf := Pl!1;
    for pr in Factorization(F) do
        if IsOdd(pr[2]) then sf *:= pr[1]; end if;
    end for;
    // unit = lcF / lc(monic sf-part times even parts) -- recover exactly:
    even := F/(lcF*sf);   // monic square rational function? compute unit properly
    // F = lcF * sf * (even monic part); the model y^2 = F/den^2 == lcF*sf mod squares
    csq := Q!lcF;
    // reduce csq modulo squares
    cnum := Numerator(csq)*Denominator(csq);
    csf := Sign(cnum);
    for pe in Factorization(Z!Abs(cnum)) do
        if IsOdd(pe[2]) then csf *:= pe[1]; end if;
    end for;
    sf := csf*sf;
    dn := LCM([Denominator(co) : co in Coefficients(sf)]);
    return Pl!(dn^2*sf);
end function;

for ch in chars do
    nm := ch[1]; fl := ch[2];
    F := SquarefreeModel(fl);
    dg := Degree(F);
    printf "== character %o: y^2 = deg-%o squarefree %o\n", nm, dg, F;
    if dg le 2 then
        printf "   genus 0 (conic/rational)\n";
        continue;
    end if;
    if dg in {3,4} then
        // genus 1: find a rational point to convert
        C := HyperellipticCurve(F);
        printf "   genus 1 candidate; searching small points ...\n";
        pts := Points(C : Bound := 3000);
        if #pts eq 0 then
            printf "   no small point found; possibly ELS-failing quotient -- check local solvability!\n";
            // quick local check at small primes
            for pp in [2,3,5,7,11,13] do
                Qp := pAdicField(pp, 30);
                found := false;
                for aa in [0..pp^2-1] do
                    lv := Q!aa;
                    vv := Evaluate(F, lv);
                    if vv eq 0 then found := true; break; end if;
                    if IsSquare(Qp!vv) then found := true; break; end if;
                end for;
                if not found then
                    printf "   LOCALLY INSOLUBLE AT p=%o (quotient torsor obstruction!)\n", pp;
                end if;
            end for;
            continue;
        end if;
        pt := pts[1];
        E, mp := EllipticCurve(C, pt);
        Emin, _ := MinimalModel(E);
        r1, r2 := RankBounds(Emin);
        Tt, _ := TorsionSubgroup(Emin);
        printf "   elliptic: %o  rank bounds [%o,%o] torsion %o%o\n",
            aInvariants(Emin), r1, r2, Invariants(Tt),
            (r2 eq 0 select "   <== RANK ZERO KILLER CANDIDATE" else "");
    else
        printf "   genus >= 2 (deg %o): point-finite by Faltings; Chabauty target\n", dg;
        // report small points for the record
        C := HyperellipticCurve(F);
        pts := Points(C : Bound := 2000);
        printf "   small points (bound 2000): %o\n",
            [<pt[1],pt[2],pt[3]> : pt in pts];
    end if;
end for;
print "DONE";
quit;
