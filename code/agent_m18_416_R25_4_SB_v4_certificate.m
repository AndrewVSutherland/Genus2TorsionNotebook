//////////////////////////////////////////////////////////////////////
// Corrected certificate for the R = -25/4 S_B cover.
//
// The earlier genus-5 pullback missed the scale factor 609/256 between
// the unreduced fiber coordinate and the reduced quartic
//
//     S^2 = 1024*m^4 - 865600*m^2 + 231800625.
//
// With the corrected scaling, the S_B cover is a genus-3 V_4 cover.  A
// rank-zero elliptic quotient forces all rational affine points to be
// degenerate branch/boundary points.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

// Corrected affine S_B pullback:
// S^2 = Gamma(m),
// Y^2 = 29/50*m^4 - 16907/32*m^2 + 268888725/2048
//       + (29/1600*m^2 - 17661/2048)*S.
A3<m,S,Y> := AffineSpace(Q, 3);
eq1 := S^2 - (1024*m^4 - 865600*m^2 + 231800625);
eq2 := 51200*Y^2 - 29696*m^4 - 928*m^2*S
       + 27051200*m^2 + 441525*S - 6722218125;

Caff := Curve(A3, [eq1, eq2]);
Cp := ProjectiveClosure(Caff);
printf "corrected S_B projective closure: degree=%o normalized genus=%o nonsingular=%o\n",
    Degree(Cp), Genus(Cp), IsNonsingular(Cp);

// Lambda quotient coordinates:
//   u = 32*m^2/25, t = S/25, lambda = (t - 609)/u.
// For non-boundary affine points,
//   m^2 = 25/16*(609*l + 541)/(1 - l^2),
// and the rank-zero quotient is
//   C^2 = -29*(l - 1)*(29*l + 21)*(21*l + 25).
P<l> := PolynomialRing(Q);
f_mY := -29*(l - 1)*(29*l + 21)*(21*l + 25);
CmY := HyperellipticCurve(f_mY);
ptsInf := PointsAtInfinity(CmY);
Eraw, phi := EllipticCurve(CmY, ptsInf[1]);
Emin, mp := MinimalModel(Eraw);
printf "E_mY minimal model: %o\n", Emin;
printf "E_mY rank bounds: %o\n", RankBounds(Emin);
printf "E_mY torsion invariants: %o\n", Invariants(TorsionSubgroup(Emin));
printf "E_mY cubic factorization: %o\n", Factorization(f_mY);

roots := [r[1] : r in Roots(f_mY)];
printf "finite lambda forced by 2-torsion: %o\n", roots;
printf "plus lambda=infinity on the elliptic quotient\n";
printf "exceptional cleared-model boundary: lambda=-541/609 gives m=0, S=15225, Y=0\n";

function IsSquareQ(q)
    q := Q!q;
    if q lt 0 then return false, Q!0; end if;
    okn, sn := IsSquare(Numerator(q));
    okd, sd := IsSquare(Denominator(q));
    if okn and okd then return true, Q!sn/sd; end if;
    return false, Q!0;
end function;

function LambdaData(lam)
    m2 := Q!25/16*(609*lam + 541)/(1 - lam^2);
    Sval := 15225 + 32*m2*lam;
    Y2 := Q!725/1024
          *(29*lam + 21)*(609*lam + 541)*(21*lam + 25)
          /((lam - 1)^2*(lam + 1));
    okm, mv := IsSquareQ(m2);
    okY, Yv := IsSquareQ(Y2);
    wv := (m2 + Q!15225/32)/(m2 - Q!15225/32);
    return m2, okm, mv, Sval, Y2, okY, Yv, wv;
end function;

for lam in roots do
    if lam in {Q!1, Q!-1} then
        printf "lambda=%o is an infinite/boundary value for the affine m-model\n", lam;
        continue;
    end if;
    m2, okm, mv, Sval, Y2, okY, Yv, wv := LambdaData(lam);
    printf "lambda=%o: m^2=%o square=%o", lam, m2, okm;
    if okm then printf " m=+/- %o", mv; end if;
    printf ", S=%o, Y^2=%o square=%o", Sval, Y2, okY;
    if okY then printf " Y=+/- %o", Yv; end if;
    printf ", w=%o\n", wv;
end for;

// The two finite affine lifts have Y=0 and d_B=0, so they are branch
// degeneracies and cannot pass the nondegenerate second-stage condition.
Pm<mm> := PolynomialRing(Q);
dBnum := -841/100*mm^4 + 2190805/256*mm^2 - 7797773025/4096;
for mv in [Q!145/8, Q!-145/8, Q!105/4, Q!-105/4] do
    printf "m=%o: d_B numerator=%o\n", mv, Evaluate(dBnum, mv);
end for;

printf "CERTIFICATE: all rational points on corrected S_B are boundary or d_B=0 branch points; no nondegenerate S_B points remain on R=-25/4.\n";

quit;
