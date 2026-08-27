//////////////////////////////////////////////////////////////////////
// Corrected S_B V_4 quotient scan for ELS fibers in the M_1(8,4)
// P_R-halving problem.
//
// For a fixed R, the C2/ELS fiber is
//
//     S^2 = Gamma(m),     Gamma even quartic in m.
//
// The corrected S_B pullback is
//
//     Y^2 = 2*alpha_B*m^2
//           + sqrt(content(Gamma_raw))*2*h_B*m^2/(m^2-K)*S.
//
// Setting X=m^2 and parameterizing the conic S^2=Gamma(X) gives a V_4
// cover over P^1_lambda.  The three elliptic quotients correspond to
// square classes X(lambda), Y2(lambda), and X(lambda)*Y2(lambda).
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

if not assigned Case then Case := "all"; end if;
if not assigned SearchHeight then SearchHeight := 0; end if;
if Type(SearchHeight) eq MonStgElt then SearchHeight := StringToInteger(SearchHeight); end if;

cases := [
    <Q!-8, Q!-28, "R=-8">,
    <Q!-25/4, Q!-25/2, "R=-25/4">,
    <Q!-29/8, Q!-29/4, "R=-29/8">
];
if Case ne "all" then
    if Case eq "-8" then
        cases := [cases[1]];
    elif Case eq "-25/4" then
        cases := [cases[2]];
    elif Case eq "-29/8" then
        cases := [cases[3]];
    else
        error "unknown Case";
    end if;
end if;

P<m> := PolynomialRing(Q);
K2<R,mm> := RationalFunctionField(Q, 2);
PX<x> := PolynomialRing(K2);

Ksym := -2*R*(R^2-1);
wsym := (mm^2 + Ksym)/(mm^2 - Ksym);
tsym := (2*R^2 + (1-wsym^2)*R - 2*wsym^2)/(4*(wsym^2-1));
c4 := R + 2 + 4*tsym;
Bpol := c4*x^2 + (R^2 + 4*R + 1 + 8*tsym)*x + (2*R^2 + R + 4*tsym);
XR := -c4*R;
Bt := PX![K2!co : co in Coefficients(c4*Evaluate(Bpol, x/c4))];
alphaBsym := XR + Coefficient(Bt,1)/2;
dBsym := Discriminant(Bt);
NBsym := Evaluate(Bt, XR);
Gsym := 2*(R^2-1)*(R*(2*R+1) - wsym^2*(R+2));
okBsym, hBsym := IsSquare(NBsym/Gsym);
assert okBsym;

function PolyFromRat(rat)
    rat := P!Numerator(rat) / Q!Denominator(rat);
    return P!rat;
end function;

function SquareQuotient(poly)
    PXq<Xq> := PolynomialRing(Q);
    out := PXq!0;
    for i in [0..Degree(poly)] do
        c := Coefficient(poly, i);
        if c ne 0 then
            assert i mod 2 eq 0;
            out +:= c*Xq^(i div 2);
        end if;
    end for;
    return out;
end function;

function SquarefreeInteger(n)
    n := Z!n;
    if n eq 0 then return Z!0; end if;
    sgn := n lt 0 select -1 else 1;
    n := Abs(n);
    out := Z!1;
    for pe in Factorization(n) do
        if pe[2] mod 2 eq 1 then out *:= pe[1]; end if;
    end for;
    return sgn*out;
end function;

function SquarefreeRational(c)
    c := Q!c;
    if c eq 0 then return Q!0; end if;
    return Q!SquarefreeInteger(Numerator(c))*Q!SquarefreeInteger(Denominator(c));
end function;

function SquarefreeCoreRat(rat, L)
    num := L!Numerator(rat);
    den := L!Denominator(rat);
    f := num*den;
    cont := RationalGCD(Coefficients(f));
    f /:= cont;
    // Magma's polynomial factorization reports monic factors, so the
    // primitive leading coefficient is part of the squareclass too.
    core := L!SquarefreeRational(cont*LeadingCoefficient(f));
    for pe in Factorization(f) do
        if pe[2] mod 2 eq 1 then
            core *:= pe[1];
        end if;
    end for;
    return core;
end function;

function IsSquareQ(q)
    q := Q!q;
    if q lt 0 then return false, Q!0; end if;
    okn, sn := IsSquare(Numerator(q));
    okd, sd := IsSquare(Denominator(q));
    if okn and okd then return true, Q!sn/sd; end if;
    return false, Q!0;
end function;

procedure AnalyzeQuotient(label, fcore)
    printf "  quotient %o squarefree f(lambda) = %o\n", label, fcore;
    printf "    degree=%o factorization=%o\n", Degree(fcore), Factorization(fcore);
    if Degree(fcore) lt 3 or Degree(fcore) gt 4 then
        printf "    not genus-one degree 3/4; skipping elliptic model\n";
        return;
    end if;
    C := HyperellipticCurve(fcore);
    pts := [];
    ptsInf := PointsAtInfinity(C);
    if #ptsInf gt 0 then
        Append(~pts, ptsInf[1]);
    else
        for rt in Roots(fcore) do
            Append(~pts, C![rt[1], Q!0, Q!1]);
            break;
        end for;
    end if;
    if #pts eq 0 then
        printf "    no obvious rational base point for elliptic conversion\n";
        return;
    end if;
    try
        Eraw, phi := EllipticCurve(C, pts[1]);
        Emin, mp := MinimalModel(Eraw);
        printf "    minimal model=%o\n", Emin;
        printf "    rank bounds=%o\n", RankBounds(Emin);
        printf "    torsion=%o\n", Invariants(TorsionSubgroup(Emin));
    catch e
        printf "    elliptic conversion/rank failed: %o\n", e`Object;
    end try;
end procedure;

procedure AnalyzeCase(Rv, mex, label)
    printf "\n================ %o ================\n", label;
    Kv := -2*Rv*(Rv^2-1);
    wv := (m^2 + Kv)/(m^2 - Kv);
    G := 2*(Rv^2-1)*(Rv*(2*Rv+1) - wv^2*(Rv+2));
    rawGamma := Numerator(G*(m^2-Kv)^2);
    gammaContent := RationalGCD(Coefficients(rawGamma));
    Gamma := rawGamma/gammaContent;
    okContent, sqrtContent := IsSquare(gammaContent);
    printf "R=%o K=%o mex=%o\n", Rv, Kv, mex;
    printf "Gamma=%o\n", Gamma;
    printf "raw content=%o sqrtContent=%o ok=%o\n", gammaContent, sqrtContent, okContent;
    if not okContent then
        printf "  content is not rational square; skipping corrected model\n";
        return;
    end if;

    alphaB := Evaluate(Numerator(alphaBsym), [Rv,m]) / Evaluate(Denominator(alphaBsym), [Rv,m]);
    dB := Evaluate(Numerator(dBsym), [Rv,m]) / Evaluate(Denominator(dBsym), [Rv,m]);
    hB := Evaluate(Numerator(hBsym), [Rv,m]) / Evaluate(Denominator(hBsym), [Rv,m]);
    alphaP := PolyFromRat(alphaB*m^2);
    dP := PolyFromRat(dB*m^2);
    coefP := PolyFromRat(sqrtContent*2*m^2*hB/(m^2-Kv));
    AP := 2*alphaP;
    printf "alphaP=%o\n", alphaP;
    printf "coefP(S)=%o\n", coefP;
    printf "dP=%o\n", dP;

    gammaX := SquareQuotient(Gamma);
    AX := SquareQuotient(AP);
    BX := SquareQuotient(coefP);
    dX := SquareQuotient(dP);
    printf "gammaX=%o\n", gammaX;
    printf "Y2 = (%o) + (%o)*S\n", AX, BX;

    okSconst, Sconst := IsSquareQ(Evaluate(gammaX, Q!0));
    if okSconst then
        X0 := Q!0;
        S0 := Sconst;
        printf "base point: boundary X0=0 S0=%o\n", S0;
    else
        X0 := mex^2;
        okS0, S0 := IsSquareQ(Evaluate(Gamma, mex));
        printf "base point: example X0=%o S0 square=%o S0=%o\n", X0, okS0, S0;
        if not okS0 then
            printf "  no rational base point; skipping\n";
            return;
        end if;
    end if;

    L<lam> := PolynomialRing(Q);
    RL<X> := PolynomialRing(L);
    gammaRL := Evaluate(gammaX, X);
    Fline := (L!S0 + lam*(X - L!X0))^2 - gammaRL;
    qline := ExactQuotient(Fline, X - L!X0);
    cX := Coefficient(qline, 1);
    c0 := Coefficient(qline, 0);
    Xrat := -c0/cX;
    Srat := L!S0 + lam*(Xrat - L!X0);
    Y2rat := Evaluate(AX, Xrat) + Evaluate(BX, Xrat)*Srat;
    printf "X(lambda)=%o\n", Xrat;
    printf "S(lambda)=%o\n", Srat;
    printf "Y2(lambda)=%o\n", Y2rat;

    fX := SquarefreeCoreRat(Xrat, L);
    fY := SquarefreeCoreRat(Y2rat, L);
    fXY := SquarefreeCoreRat(Xrat*Y2rat, L);
    AnalyzeQuotient("E_m", fX);
    AnalyzeQuotient("E_Y", fY);
    AnalyzeQuotient("E_mY", fXY);

    // If a rank-zero quotient has only 2-torsion, the finite lambda
    // candidates are among the rational roots of its squarefree cubic.
    for pair in [<"E_m", fX>, <"E_Y", fY>, <"E_mY", fXY>] do
        f := pair[2];
        if Degree(f) in {3,4} then
            C := HyperellipticCurve(f);
            pts := PointsAtInfinity(C);
            if #pts eq 0 and #Roots(f) gt 0 then
                pts := [C![Roots(f)[1][1], Q!0, Q!1]];
            end if;
            if #pts gt 0 then
                try
                    Eraw, phi := EllipticCurve(C, pts[1]);
                    Emin, mp := MinimalModel(Eraw);
                    rb := RankBounds(Emin);
                    tors := Invariants(TorsionSubgroup(Emin));
                    if rb eq 0 and tors eq [2,2] and #Roots(f) gt 0 then
                        printf "  candidate killer quotient %o: roots %o plus infinity/boundary\n", pair[1], [rt[1] : rt in Roots(f)];
                        for rt in Roots(f) do
                            lv := rt[1];
                            if Evaluate(Denominator(Xrat), lv) eq 0 then
                                printf "    lambda=%o: denominator boundary\n", lv;
                                continue;
                            end if;
                            Xv := Evaluate(Xrat, lv);
                            Sv := Evaluate(Srat, lv);
                            Y2v := Evaluate(Y2rat, lv);
                            okm, mv := IsSquareQ(Xv);
                            okY, Yv := IsSquareQ(Y2v);
                            wval := (Xv + Kv)/(Xv - Kv);
                            printf "    lambda=%o: X=m^2=%o square=%o", lv, Xv, okm;
                            if okm then printf " m=+/- %o", mv; end if;
                            printf ", S=%o, Y2=%o square=%o", Sv, Y2v, okY;
                            if okY then printf " Y=+/- %o", Yv; end if;
                            printf ", w=%o", wval;
                            if okm then
                                printf ", dP(+m)=%o", Evaluate(dP, mv);
                            end if;
                            printf "\n";
                        end for;
                    end if;
                catch e
                    printf "  quotient %o candidate check skipped: %o\n", pair[1], e`Object;
                end try;
            end if;
        end if;
    end for;

    if SearchHeight gt 0 then
        printf "  lambda search height <= %o for nondegenerate S_B points\n", SearchHeight;
        seen := {};
        hits := 0;
        for den in [1..SearchHeight] do
            for num in [-SearchHeight..SearchHeight] do
                if GCD(num, den) ne 1 then continue; end if;
                lv := Q!num/den;
                if lv in seen then continue; end if;
                Include(~seen, lv);
                if Evaluate(Denominator(Xrat), lv) eq 0 or Evaluate(Denominator(Y2rat), lv) eq 0 then
                    continue;
                end if;
                Xv := Evaluate(Xrat, lv);
                Y2v := Evaluate(Y2rat, lv);
                okm, mv := IsSquareQ(Xv);
                okY, Yv := IsSquareQ(Y2v);
                if okm and okY and Xv ne 0 then
                    dval := Evaluate(dP, mv);
                    if dval ne 0 then
                        Sv := Evaluate(Srat, lv);
                        wval := (Xv + Kv)/(Xv - Kv);
                        hits +:= 1;
                        printf "    HIT lambda=%o: m=+/- %o S=%o Y=+/- %o w=%o dP=%o\n",
                            lv, mv, Sv, Yv, wval, dval;
                    end if;
                end if;
            end for;
        end for;
        printf "  lambda search nondegenerate hits=%o\n", hits;
    end if;
end procedure;

for c in cases do
    AnalyzeCase(c[1], c[2], c[3]);
end for;

quit;
