// claude_z96_family_setup.m — shared setup: the reconstructed Elkies [32] family
// as a one-parameter family of quintics over Q(u), via parametrization of the
// genus-0 component C32(z,r)=0 and the recorded scale formula p=(F2G0-G2F0)/(G2F1-F2G1).
// Loaded by the Z/96 (32 x 3-contact) and Z/160 (32 x second-5-contact) deciders.

function BuildFamily()
  Q := Rationals();
  Pzr<z, r> := PolynomialRing(Q, 2);
  C32 := 3*z^4*r^4 + 9*z^4*r^3 - 16*z^3*r^4 + 10*z^4*r^2 - 56*z^3*r^3 + 32*z^2*r^4
       + 5*z^4*r - 72*z^3*r^2 + 144*z^2*r^3 - 64*z*r^4 + z^4 - 40*z^3*r + 208*z^2*r^2
       - 224*z*r^3 + 80*r^4 - 8*z^3 + 120*z^2*r - 288*z*r^2 + 160*r^3 + 24*z^2
       - 160*z*r + 160*r^2 - 32*z + 80*r + 16;
  // parametrize the genus-0 plane curve from the printed point (11/6, 1/15)
  A2 := AffineSpace(Pzr);
  Cc := Curve(A2, C32);
  Pc := ProjectiveClosure(Cc);
  ptp := Pc![11/6, 1/15, 1];
  para := Parametrization(Pc, ptp);
  dp := DefiningPolynomials(para);
  Ku<u> := FunctionField(Q);
  zu := Evaluate(dp[1], [u, 1]) / Evaluate(dp[3], [u, 1]);
  ru := Evaluate(dp[2], [u, 1]) / Evaluate(dp[3], [u, 1]);
  // scale p from the recorded common-root formula: build F,G as P-quadratics
  PzrP<zz, rr, PP> := PolynomialRing(Q, 3);
  Fbig := -zz^8*PP^2*rr^4 + 8*zz^7*PP^2*rr^5 - 2*zz^8*PP^2*rr^3 + 24*zz^7*PP^2*rr^4 - 48*zz^6*PP^2*rr^5
    - zz^8*PP^2*rr^2 + 24*zz^7*PP^2*rr^3 - 120*zz^6*PP^2*rr^4 + 96*zz^5*PP^2*rr^5 + 8*zz^7*PP^2*rr^2
    - 104*zz^6*PP^2*rr^3 + 8*zz^6*PP*rr^4 + 272*zz^5*PP^2*rr^4 - 32*zz^5*PP*rr^5 - 64*zz^4*PP^2*rr^5
    - 32*zz^6*PP^2*rr^2 + 24*zz^6*PP*rr^3 + 256*zz^5*PP^2*rr^3 - 128*zz^5*PP*rr^4 - 336*zz^4*PP^2*rr^4
    + 24*zz^6*PP*rr^2 + 80*zz^5*PP^2*rr^2 - 192*zz^5*PP*rr^3 - 384*zz^4*PP^2*rr^3 + 64*zz^4*PP*rr^4
    + 192*zz^3*PP^2*rr^4 + 128*zz^3*PP*rr^5 + 8*zz^6*PP*rr - 128*zz^5*PP*rr^2 - 128*zz^4*PP^2*rr^2
    + 192*zz^4*PP*rr^3 + 320*zz^3*PP^2*rr^3 - 16*zz^4*rr^4 + 320*zz^3*PP*rr^4 - 64*zz^2*PP^2*rr^4
    - 32*zz^5*PP*rr + 192*zz^4*PP*rr^2 + 128*zz^3*PP^2*rr^2 - 64*zz^4*rr^3 + 192*zz^3*PP*rr^3
    - 128*zz^2*PP^2*rr^3 - 384*zz^2*PP*rr^4 + 64*zz^4*PP*rr - 96*zz^4*rr^2 - 64*zz^3*PP*rr^2
    - 64*zz^2*PP^2*rr^2 - 640*zz^2*PP*rr^3 + 128*zz^2*rr^4 + 256*zz*PP*rr^4 - 64*zz^4*rr
    - 64*zz^3*PP*rr - 256*zz^2*PP*rr^2 + 384*zz^2*rr^3 + 512*zz*PP*rr^3 - 16*zz^4 + 384*zz^2*rr^2
    + 256*zz*PP*rr^2 - 256*rr^4 + 128*zz^2*rr - 512*rr^3 - 256*rr^2;
  Gbig := zz^7*PP^2*rr^4 + 3*zz^7*PP^2*rr^3 - 12*zz^6*PP^2*rr^4 + 3*zz^7*PP^2*rr^2 - 30*zz^6*PP^2*rr^3
    + 36*zz^5*PP^2*rr^4 + zz^7*PP^2*rr - 26*zz^6*PP^2*rr^2 + 100*zz^5*PP^2*rr^3 - 4*zz^5*PP*rr^4
    - 24*zz^4*PP^2*rr^4 - 8*zz^6*PP^2*rr + 92*zz^5*PP^2*rr^2 - 16*zz^5*PP*rr^3 - 144*zz^4*PP^2*rr^3
    + 28*zz^5*PP^2*rr - 24*zz^5*PP*rr^2 - 168*zz^4*PP^2*rr^2 + 16*zz^4*PP*rr^3 + 96*zz^3*PP^2*rr^3
    + 48*zz^3*PP*rr^4 - 16*zz^5*PP*rr - 56*zz^4*PP^2*rr + 48*zz^4*PP*rr^2 + 160*zz^3*PP^2*rr^2
    + 128*zz^3*PP*rr^3 - 32*zz^2*PP^2*rr^3 - 4*zz^5*PP + 48*zz^4*PP*rr + 64*zz^3*PP^2*rr
    + 96*zz^3*PP*rr^2 - 64*zz^2*PP^2*rr^2 - 192*zz^2*PP*rr^3 + 16*zz^4*PP - 32*zz^2*PP^2*rr
    - 320*zz^2*PP*rr^2 + 32*zz^2*rr^3 + 128*zz*PP*rr^3 - 16*zz^3*PP - 128*zz^2*PP*rr + 96*zz^2*rr^2
    + 256*zz*PP*rr^2 + 96*zz^2*rr + 128*zz*PP*rr - 128*rr^3 + 32*zz^2 - 256*rr^2 - 128*rr;
  // coefficients as polynomials in PP
  UP<T> := PolynomialRing(FieldOfFractions(Pzr));
  FP := &+[ UP | Evaluate(Coefficient(Fbig, 3, i), [Pzr.1, Pzr.2, 0]) * T^i : i in [0..2] ];
  GP := &+[ UP | Evaluate(Coefficient(Gbig, 3, i), [Pzr.1, Pzr.2, 0]) * T^i : i in [0..2] ];
  F2 := Coefficient(FP, 2); F1 := Coefficient(FP, 1); F0 := Coefficient(FP, 0);
  G2 := Coefficient(GP, 2); G1 := Coefficient(GP, 1); G0 := Coefficient(GP, 0);
  pnum := F2*G0 - G2*F0;  pden := G2*F1 - F2*G1;
  // sanity at printed point
  pval := Evaluate(Numerator(pnum), [11/6, 1/15]) / Evaluate(Numerator(pden), [11/6, 1/15])
          * Evaluate(Denominator(pden), [11/6, 1/15]) / Evaluate(Denominator(pnum), [11/6, 1/15]);
  assert pval eq -1440/11;
  // specialize along the parametrization
  pu := Evaluate(Numerator(pnum), [zu, ru]) / Evaluate(Numerator(pden), [zu, ru])
        * Evaluate(Denominator(pden), [zu, ru]) / Evaluate(Denominator(pnum), [zu, ru]);
  au := zu * pu;
  cu := (-au^4 + 4*au^3*pu - 8*au^2*pu^2 + 8*au*pu^3 + 8*au^2*pu - 16*pu^3)/(4*au^2*pu);
  bu := au + cu - 1 - pu;
  Px<x> := PolynomialRing(Ku);
  hu := au*x^3 + bu*x^2 + cu*x + 1;
  fu := hu^2 - au^2*x^5*(x+1);
  assert Degree(fu) eq 5;
  // clear denominators: ftil = D^2 * fu integral in u
  D := LCM([Denominator(Coefficient(fu, i)) : i in [0..5]]);
  ftil := D^2 * fu;
  PxZ<xx> := PolynomialRing(Ku);
  return Ku, PxZ, ftil, fu, zu, ru;
end function;
