//////////////////////////////////////////////////////////////////////
// Stage 2: hyperelliptic model + Mordell-Weil rank bounds for the
// genus-2 S3-quotient E4 of the signed genus-10 relative-3 cover.
//
// Input: data/contact6_m612_E4_mp4Q.txt — the exact monic quartic
// mp4Q over Q(e) written by code/contact6_m612_relative3_E4_exact.m
// (Sprintf "%m" format).  E4 = Q(e)[Y]/(mp4Q), certified genus 2.
//
// Steps:
//  1. rebuild the plane model G(Y,e)=0 over Q (clear denominators);
//  2. convert to a hyperelliptic model y^2 = f(x) over Q
//     (IsHyperelliptic, with a function-field Riemann-Roch fallback);
//  3. J = Jacobian: RankBounds, TorsionSubgroup, small-height points;
//  4. if rank bound <= 1: attempt Chabauty.
// Every rational point of the [6,12] P8 relative-3 cover maps to a
// rational point of this curve, so its point set gates the route.
//
// Usage: magma -b code/contact6_m612_relative3_E4_rank.m
//////////////////////////////////////////////////////////////////////

SetColumns(0); SetSeed(1);
SetMemoryLimit(6*10^9);
Q := Rationals();
K<e> := FunctionField(Q); Kz<zz> := PolynomialRing(K);

// load mp4Q
mp4Q := eval Read("data/contact6_m612_E4_mp4Q.txt");
mp4Q := Kz!mp4Q;
assert Degree(mp4Q) eq 4 and IsMonic(mp4Q);
printf "loaded mp4Q, degree %o\n", Degree(mp4Q);

// plane model over Q: clear denominators of the coefficients
R2<Y,E> := PolynomialRing(Q, 2);
D := LCM([Denominator(Coefficient(mp4Q,i)) : i in [0..4]]);
G := R2!0;
for i in [0..4] do
    c := Coefficient(mp4Q, i) * K!D;   // polynomial in e now
    cn := Numerator(c);
    G +:= (&+[R2 | Coefficient(cn,h)*E^h : h in [0..Degree(cn)]]) * Y^i;
end for;
printf "plane model: total degree %o, terms %o\n", TotalDegree(G), #Terms(G);
A2 := AffineSpace(R2);
C := Curve(A2, G);
gC := Genus(C);
printf "plane model genus = %o\n", gC;
assert gC eq 2;

// hyperelliptic model
okH := false;
try
    bH, H, mpH := IsHyperelliptic(C);
    if bH then
        okH := true;
        printf "IsHyperelliptic: %o\n", H;
    end if;
catch ee
    printf "IsHyperelliptic failed: %o\n", ee`Object;
end try;

if not okH then
    // function-field Riemann-Roch fallback: degree-2 map from canonical class
    print "fallback: function-field construction";
    FC<fx,fy> := FunctionField(C);
    FF := AlgorithmicFunctionField(FC);
    KD := CanonicalDivisor(FF);
    B := Basis(KD);                      // L(K) is 2-dim for genus 2
    printf "dim L(K) = %o\n", #B;
    assert #B eq 2;
    xx := B[2]/B[1];                     // degree-2 map to P^1
    mpy := MinimalPolynomial(xx);
    // find y : [FF : Q(xx)] = 2; use any element not in Q(xx)
    // canonical approach: FF as degree-2 extension of Q(xx)
    Fx<t> := FunctionField(Q);
    // represent FF over Q(xx): minimal polynomial of a generator
    gen := FF.1;
    mpg := MinimalPolynomial(gen, sub< FF | xx >);
    print "relative minpoly degree:", Degree(mpg);
end if;

if okH then
    // simplify + Jacobian arithmetic
    H2, mpsimp := SimplifiedModel(H);
    printf "simplified: %o\n", H2;
    fH, hH := HyperellipticPolynomials(H2);
    printf "f = %o\n", fH;
    J := Jacobian(H2);
    tt := Cputime();
    rlo, rhi := RankBounds(J);
    printf "RANK BOUNDS: %o <= rank <= %o  (%.1os)\n", rlo, rhi, Cputime(tt);
    T, mpT := TorsionSubgroup(J);
    printf "torsion subgroup: %o\n", Invariants(T);
    pts := Points(H2 : Bound := 5000);
    printf "small points (height 5000): %o\n", pts;
    if rhi eq 0 then
        print "RANK 0: rational points = torsion pullback; Chabauty0 applies";
    elif rhi eq 1 then
        print "RANK <= 1: Chabauty applicable once a generator is known";
    end if;
end if;
print "E4_RANK_STAGE_DONE";
quit;
