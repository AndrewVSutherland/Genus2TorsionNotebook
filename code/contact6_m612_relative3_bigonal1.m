//////////////////////////////////////////////////////////////////////
// Bigonal construction, step 1: express Y (the E8-cover function) on
// the Weierstrass model of E4 as Y = (a(x) + b(x)*y)/d(x), exactly.
// Method: explicit map mpH: C(Y,e) -> H(x,y); fit the linear model
// through C(F_p) points for two primes; rational-reconstruct; VERIFY
// exactly modulo the curve equation G(Y,e).
// Then N := a^2 - b^2*fH (norm of Y), its squarefree part (expected
// quadratic: the two branch points of E8->E4), printed for step 2.
// Usage: magma -b code/contact6_m612_relative3_bigonal1.m
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(1); SetMemoryLimit(6*10^9);
Q := Rationals();
K<e> := FunctionField(Q); Kz<zz> := PolynomialRing(K);
mp4Q := Kz ! eval Read("data/contact6_m612_E4_mp4Q.txt");
R2<Y,E> := PolynomialRing(Q, 2);
D := LCM([Denominator(Coefficient(mp4Q,i)) : i in [0..4]]);
G := R2!0;
for i in [0..4] do
    c := Coefficient(mp4Q, i) * K!D; cn := Numerator(c);
    G +:= (&+[R2 | Coefficient(cn,h)*E^h : h in [0..Degree(cn)]]) * Y^i;
end for;
C := Curve(AffineSpace(R2), G);
bH, H, mpH := IsHyperelliptic(C);
assert bH;
fH, hH := HyperellipticPolynomials(H); assert hH eq 0;
printf "H: y^2 = %o\n", fH;
dp := DefiningPolynomials(mpH);
printf "map C->H defining polynomials (%o of them):\n", #dp;
for i in [1..#dp] do printf "  p%o = %o\n", i, dp[i]; end for;

// fit Y = (a(x)+b(x)y)/d(x) through F_p points, two primes, then verify
Px<x> := PolynomialRing(Q);
function FitModP(p, dp, G, dega, degb, degd)
    Fp := GF(p); R2p := PolynomialRing(Fp, 2);
    Gp := R2p ! G;
    dpp := [R2p ! g : g in dp];
    Pp := PolynomialRing(Fp);
    rows := [];
    for e0 in Fp do
        if e0 eq 0 then continue; end if;
        quart := Pp ! [Evaluate(Coefficient(Gp, 1, i), [0, e0]) : i in [0..4]];
        if quart eq 0 then continue; end if;
        for r in Roots(quart) do
            Y0 := r[1];
            // image point on weighted-(1,3,1) model: x=p1/p3, y=p2/p3^3
            v1 := Evaluate(dpp[1], [Y0,e0]); v2 := Evaluate(dpp[2], [Y0,e0]); v3 := Evaluate(dpp[3], [Y0,e0]);
            if v3 eq 0 then continue; end if;
            x0 := v1/v3; y0 := v2/v3^3;
            // row: d-part: -Y0*x0^k (k<=degd) ; a-part: x0^i ; b-part: x0^j*y0
            row := [ -Y0*x0^k : k in [0..degd] ] cat [ x0^i : i in [0..dega] ] cat [ x0^j*y0 : j in [0..degb] ];
            Append(~rows, row);
            if #rows ge 3*(dega+degb+degd+5) then break e0; end if;
        end for;
    end for;
    M := Matrix(Fp, #rows, degd+dega+degb+3, &cat rows);
    NS := Nullspace(Transpose(M));
    return Dimension(NS), NS;
end function;

found := false;
for degs in [<4,2,4>, <5,3,5>, <6,4,6>, <7,5,7>, <8,6,8>] do
    dega := degs[1]; degb := degs[2]; degd := degs[3];
    dim1, NS1 := FitModP(10007, dp, G, dega, degb, degd);
    if dim1 eq 0 then continue; end if;
    printf "degrees a<=%o b<=%o d<=%o : nullspace dim %o (mod 10007)\n", dega, degb, degd, dim1;
    // take a basis vector mod two primes and rationally reconstruct via CRT
    dim2, NS2 := FitModP(10009, dp, G, dega, degb, degd);
    if dim2 ne dim1 then continue; end if;
    v1 := Basis(NS1)[1]; v2 := Basis(NS2)[1];
    // normalize both at the first nonzero coordinate
    i0 := Min([i : i in [1..Ncols(v1)] | v1[i] ne 0]);
    v1 := v1/v1[i0]; v2 := v2/v2[i0];
    N12 := 10007*10009; ZN := Integers(N12);
    co := [];
    okall := true;
    for i in [1..Ncols(v1)] do
        xr := CRT([Integers()!v1[i], Integers()!v2[i]], [10007, 10009]);
        okr, q := RationalReconstruction(ZN!xr);
        if not okr then okall := false; break; end if;
        Append(~co, q);
    end for;
    if not okall then print "  rational reconstruction failed; trying larger degrees"; continue; end if;
    dpoly := Px ! co[1..degd+1];
    apoly := Px ! co[degd+2..degd+dega+2];
    bpoly := Px ! co[degd+dega+3..degd+dega+degb+3];
    if dpoly eq 0 then continue; end if;
    // EXACT verification: d(xi)*Y - a(xi) - b(xi)*eta == 0 mod G, where
    // xi = p1/p3, eta = p2/p3^3 as rational functions on C
    FF := FieldOfFractions(R2);
    xiF := FF!dp[1]/FF!dp[3]; etaF := FF!dp[2]/(FF!dp[3])^3;
    lhs := Evaluate(dpoly, xiF)*FF!Y - Evaluate(apoly, xiF) - Evaluate(bpoly, xiF)*etaF;
    num := Numerator(lhs);
    if NormalForm(num, [G]) eq 0 then
        printf "EXACT VERIFIED: Y = (a + b*y)/d with\n  a = %o\n  b = %o\n  d = %o\n", apoly, bpoly, dpoly;
        // norm of Y
    Nx := (apoly^2 - bpoly^2*(Px!fH))/dpoly^2;
        printf "N = Y*Y^tau = (a^2 - b^2*fH)/d^2 ; numerator factorization:\n";
        printf "  num: %o\n", Factorization(Numerator(Nx));
        printf "  den: %o\n", Factorization(Denominator(Nx));
        PrintFile("data/contact6_m612_bigonal_abd.txt",
            Sprintf("<%m, %m, %m>", apoly, bpoly, dpoly) : Overwrite:=true);
        found := true;
        break;
    else
        print "  exact verification FAILED at these degrees; escalating";
    end if;
end for;
if not found then print "NO FIT FOUND up to degree 8 -- rethink"; end if;
print "BIGONAL1_DONE";
quit;
