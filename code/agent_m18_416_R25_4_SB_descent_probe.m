//////////////////////////////////////////////////////////////////////
// Probe the S_B two-cover on the ELS fiber R = -25/4.
//
// The newest note in agent_m18_416_p7_blowup_notes.md suggests
// certifying the obstruction on this fiber by doing a genuine descent on
// the small S_B cover
//
//     y^4 - 4 alpha_B(m) y^2 + d_B(m) = 0
//
// pulled back to the elliptic fiber E_R : g^2 = G(R,m).  This script
// builds that cover in a few equivalent forms and reports what Magma can
// recognize directly.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

Rv := Q!-25/4;
P<m> := PolynomialRing(Q);

K := -2*Rv*(Rv^2 - 1);
w := (m^2 + K)/(m^2 - K);
W := w^2;
G := 2*(Rv^2-1)*(Rv*(2*Rv+1) - W*(Rv+2));

// The quartic model used by the earlier ELS scripts.  The raw numerator
// has a rational-square content; the reduced coordinate S satisfies
// raw_sqrt = sqrtContent*S.
rawGamma := Numerator(G*(m^2-K)^2);
gammaContent := RationalGCD(Coefficients(rawGamma));
Gamma := rawGamma/gammaContent;
okContent, sqrtContent := IsSquare(gammaContent);
assert okContent;

// Rebuild alpha_B,d_B from the same symbolic construction to avoid
// transcription errors.
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

alphaB := Evaluate(Numerator(alphaBsym), [Rv,m]) / Evaluate(Denominator(alphaBsym), [Rv,m]);
dB := Evaluate(Numerator(dBsym), [Rv,m]) / Evaluate(Denominator(dBsym), [Rv,m]);
hB := Evaluate(Numerator(hBsym), [Rv,m]) / Evaluate(Denominator(hBsym), [Rv,m]);

printf "R = %o\n", Rv;
printf "K = %o\n", K;
printf "Gamma = %o\n", Gamma;
printf "raw Gamma content=%o sqrtContent=%o\n", gammaContent, sqrtContent;
printf "alphaB = %o\n", alphaB;
printf "dB = %o\n", dB;
printf "hB = %o\n", hB;
printf "identity alphaB^2 - dB/4 = G*hB^2: %o\n", alphaB^2 - dB/4 eq G*hB^2;
if alphaB^2 - dB/4 ne G*hB^2 then
    printf "identity difference = %o\n", alphaB^2 - dB/4 - G*hB^2;
end if;

function PolyFromRat(rat)
    rat := P!Numerator(rat) / Q!Denominator(rat);
    return P!rat;
end function;

// Plane quartic-ish S_B equation after clearing denominators in m.
// Use Y = m*y, so denominators m^2 are removed.
P2<m2,Y> := PolynomialRing(Q, 2);
alphaP := PolyFromRat(alphaB*m^2);
dP := PolyFromRat(dB*m^2);
alpha2 := Evaluate(alphaP, m2);
d2 := Evaluate(dP, m2);
SBplane := Y^4 - 4*alpha2*Y^2 + d2*m2^2;
SBplane /:= RationalGCD(Coefficients(SBplane));
printf "SB plane equation in (m,Y=m*y):\n%o\n", SBplane;
printf "total degree %o, degree_m %o, degree_Y %o\n",
    TotalDegree(SBplane), Degree(SBplane,1), Degree(SBplane,2);
printf "factorization: %o\n", Factorization(SBplane);

// Quotient by the visible involutions m -> -m and Y -> -Y:
// X = m^2, Z = Y^2.  Rational points on the original S_B curve must
// be rational points on this quotient with both X and Z rational squares.
// This is not the full obstruction: its discriminant is the ELS fiber
// again, because the square root used to solve for Z is the fiber
// coordinate up to a rational factor.
P2q<X,Zq> := PolynomialRing(Q, 2);

PX<Xu> := PolynomialRing(Q);
function SquareQuotientUnivariate(f)
    out := PX!0;
    for i in [0..Degree(f)] do
        c := Coefficient(f, i);
        if c ne 0 then
            assert i mod 2 eq 0;
            out +:= c*Xu^(i div 2);
        end if;
    end for;
    return out;
end function;

alphaXu := SquareQuotientUnivariate(alphaP);
dXu := SquareQuotientUnivariate(dP);
GammaX := SquareQuotientUnivariate(Gamma);
HX := SquareQuotientUnivariate(PolyFromRat(hB*m^2));
alphaX := Evaluate(alphaXu, X);
dX := Evaluate(dXu, X);
SBquot := Zq^2 - 4*alphaX*Zq + X*dX;
SBquot /:= RationalGCD(Coefficients(SBquot));
printf "quotient equation F(X,Z)=0, X=m^2, Z=Y^2:\n%o\n", SBquot;
printf "quotient bidegrees X=%o Z=%o\n", Degree(SBquot,1), Degree(SBquot,2);
discZ := 16*alphaXu^2 - 4*Xu*dXu;
discZ /:= RationalGCD(Coefficients(discZ));
printf "disc_Z quotient = %o\n", discZ;
printf "disc_Z degree %o factorization %o\n", Degree(discZ), Factorization(discZ);
printf "disc_Z/GammaX is square: ";
okDisc, sqDisc := IsSquare(discZ/GammaX);
printf "%o", okDisc;
if okDisc then
    printf ", sqrt=%o", sqDisc;
end if;
printf "\n";

try
    Cq := HyperellipticCurve(discZ);
    printf "quotient hyperelliptic genus=%o\n", Genus(Cq);
    rootsX := Roots(discZ);
    printf "disc roots over Q: %o\n", rootsX;
    try
        ptsq := Points(Cq : Bound := 1000);
        printf "quotient hyperelliptic points Bound 1000: %o\n", ptsq;
    catch e
        printf "Points quotient failed: %o\n", e`Object;
    end try;
    if Genus(Cq) eq 1 then
        try
            ptsq := Points(Cq : Bound := 100);
            if #ptsq gt 0 then
                Eq, phiq := EllipticCurve(Cq, ptsq[1]);
                Eqmin, mpq := MinimalModel(Eq);
                printf "quotient elliptic model: %o\n", Eqmin;
                printf "quotient rank bounds: %o %o\n", RankBounds(Eqmin);
                printf "quotient torsion: %o\n", Invariants(TorsionSubgroup(Eqmin));
            end if;
        catch e
            printf "quotient elliptic conversion failed: %o\n", e`Object;
        end try;
    end if;
catch e
    printf "quotient hyperelliptic construction failed: %o\n", e`Object;
end try;

// Affine pulled-back cover over E_R.  Let gnum = g*(m^2-K), so
// gnum^2 = Gamma and
// (m*y)^2 = m^2*(2 alphaB + 2 hB*g).
// The right side is polynomial in m and gnum.
A3<m3,gnum,Y3> := AffineSpace(Q, 3);
R3 := CoordinateRing(A3);
mm3 := R3.1; gg3 := R3.2; yy3 := R3.3;

Gamma3 := Evaluate(Gamma, mm3);
alpha_m2_3 := Evaluate(PolyFromRat(alphaB*m^2), mm3);
// Coefficient of the reduced quartic coordinate S in y^2*m^2.
coef_g_3 := Evaluate(PolyFromRat(sqrtContent*2*m^2*hB/(m^2-K)), mm3);
eq1 := gg3^2 - Gamma3;
eq2 := yy3^2 - 2*alpha_m2_3 - coef_g_3*gg3;
denoms := [Denominator(c) : c in Coefficients(eq2)];
eq2c := LCM(denoms)*eq2;
eq1c := eq1;
printf "pullback eq1: %o\n", eq1c;
printf "pullback eq2 cleared: %o\n", eq2c;

Caff := Curve(A3, [eq1c, eq2c]);
printf "affine curve defined; dimension=%o\n", Dimension(Caff);
try
    Cp := ProjectiveClosure(Caff);
    printf "projective closure degree=%o dimension=%o\n", Degree(Cp), Dimension(Cp);
    try
        printf "genus(projective closure)=%o\n", Genus(Cp);
    catch e
        printf "Genus(projective closure) failed: %o\n", e`Object;
    end try;
    try
        G1 := GenusOneModel(Cp);
        printf "GenusOneModel(Cp) succeeded, degree=%o\n", Degree(G1);
        printf "G1 equation(s): %o\n", Equations(Curve(G1));
        printf "Jacobian(G1) = %o\n", MinimalModel(Jacobian(G1));
        for pp in [2,3,5,7,11,13,17,19,23,29,31,37,41] do
            try
                ok, pt := IsLocallySolvable(G1, pp);
                printf "local G1 p=%o: %o", pp, ok;
                if ok then printf " pt=%o", pt; end if;
                printf "\n";
            catch e
                printf "local G1 p=%o failed: %o\n", pp, e`Object;
            end try;
        end for;
        try
            pts := Points(Curve(G1) : Bound := 1000);
            printf "G1 rational points Bound 1000: %o\n", pts;
        catch e
            printf "Points(G1) failed: %o\n", e`Object;
        end try;
    catch e
        printf "GenusOneModel(Cp) failed: %o\n", e`Object;
    end try;
catch e
    printf "ProjectiveClosure/geometry failed: %o\n", e`Object;
end try;

// Also ask Magma about the elliptic fiber itself.
Cfib := HyperellipticCurve(Gamma);
okpt0, sqrtpt0 := IsSquare(Evaluate(Gamma, Q!-25/2));
assert okpt0;
pt0 := Cfib![Q!-25/2, sqrtpt0, Q!1];
E, phi := EllipticCurve(Cfib, pt0);
Emin, mp := MinimalModel(E);
printf "elliptic fiber minimal model: %o\n", Emin;
printf "rank bounds: %o\n", RankBounds(Emin);
MW, toE := MordellWeilGroup(Emin);
invs := Invariants(MW);
printf "MW invariants: %o\n", invs;
freeGens := {};
for i in [1..#invs] do
    Pi := toE(MW.i);
    printf "MW generator %o invariant=%o point=%o\n", i, invs[i], Pi;
    if invs[i] eq 0 then
        Include(~freeGens, Pi);
    end if;
end for;

printf "Starting ordinary elliptic TwoDescent after removing torsion and free generators\n";
try
    SetVerbose("TwoDescent", 1);
    covers, maps := TwoDescent(Emin : RemoveTorsion := true, RemoveGens := freeGens);
    printf "TwoDescent residual num_covers=%o\n", #covers;
    for i in [1..#covers] do
        printf "residual cover %o: %o\n", i, covers[i];
        try
            pts := Points(covers[i] : Bound := 1000);
            printf "residual cover %o Points Bound 1000: %o\n", i, pts;
        catch e
            printf "residual cover %o point search failed: %o\n", i, e`Object;
        end try;
    end for;
catch e
    printf "TwoDescent residual failed: %o\n", e`Object;
end try;

try
    printf "MordellWeilShaInformation: %o\n", MordellWeilShaInformation(Emin);
catch e
    printf "MordellWeilShaInformation failed: %o\n", e`Object;
end try;
printf "done\n";
quit;
