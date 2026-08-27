//////////////////////////////////////////////////////////////////////
// Coleman/Chabauty prep: (1) explicit generator of J(C2'_min)(Q);
// (2) isomorphism C2'_min -> C2'_param (s^2 = f2'(t));
// (3) exact fit of e = (ae(x) + be(x)*y)/de(x) on E4's Weierstrass
//     model H (same fit-and-verify method as the Y-fit).
// Usage: magma -b code/contact6_m612_relative3_coleman_prep.m
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(1); SetMemoryLimit(4*10^9);
Q := Rationals(); Px<x> := PolynomialRing(Q);

// (1) generator of J(C2'_min)
fmin := -3*x^6 + 24*x^3 - 75;
Cmin := HyperellipticCurve(fmin);
Jmin := Jacobian(Cmin);
pts := Points(Jmin : Bound := 2000);
printf "Jmin points to 2000: %o\n", #pts;
gens := [P : P in pts | Order(P) eq 0];
printf "infinite-order: %o\n", #gens;
if #gens gt 0 then
    g0 := gens[1];
    printf "GENERATOR (Mumford): a = %o ; b = %o ; d = %o\n", g0[1], g0[2], g0[3];
end if;

// (2) isomorphism to the parametrized model
f2param := -1/192*x^6 + 1/32*x^5 - 5/64*x^4 + 7/16*x^3 - 69/64*x^2 + 33/32*x - 555/64;
Cpar := HyperellipticCurve(f2param);
okI, mpI := IsIsomorphic(Cmin, Cpar);
printf "IsIsomorphic(Cmin, Cpar): %o\n", okI;
if okI then print "map:", mpI; end if;

// (3) exact fit of e on H: y^2 = fH
K<e> := FunctionField(Q); Kz<zz> := PolynomialRing(K);
mp4Q := Kz ! eval Read("data/contact6_m612_E4_mp4Q.txt");
R2<Y,E> := PolynomialRing(Q, 2);
DD := LCM([Denominator(Coefficient(mp4Q,i)) : i in [0..4]]);
G := R2!0;
for i in [0..4] do
    c := Coefficient(mp4Q, i) * K!DD; cn := Numerator(c);
    G +:= (&+[R2 | Coefficient(cn,h)*E^h : h in [0..Degree(cn)]]) * Y^i;
end for;
C := Curve(AffineSpace(R2), G);
bH, H, mpH := IsHyperelliptic(C);
assert bH;
dp := DefiningPolynomials(mpH);
function FitModP(p, dp, G, dega, degb, degd)
    Fp := GF(p); R2p := PolynomialRing(Fp, 2);
    Gp := R2p ! G; dpp := [R2p ! g : g in dp];
    Pp := PolynomialRing(Fp);
    rows := [];
    for e0 in Fp do
        if e0 eq 0 then continue; end if;
        quart := Pp ! [Evaluate(Coefficient(Gp, 1, i), [0, e0]) : i in [0..4]];
        if quart eq 0 then continue; end if;
        for r in Roots(quart) do
            Y0 := r[1];
            v1 := Evaluate(dpp[1], [Y0,e0]); v2 := Evaluate(dpp[2], [Y0,e0]); v3 := Evaluate(dpp[3], [Y0,e0]);
            if v3 eq 0 then continue; end if;
            x0 := v1/v3; y0 := v2/v3^3;
            row := [ -e0*x0^k : k in [0..degd] ] cat [ x0^i : i in [0..dega] ] cat [ x0^j*y0 : j in [0..degb] ];
            Append(~rows, row);
            if #rows ge 3*(dega+degb+degd+5) then break e0; end if;
        end for;
    end for;
    M := Matrix(Fp, #rows, degd+dega+degb+3, &cat rows);
    NS := Nullspace(Transpose(M));
    return Dimension(NS), NS;
end function;
found := false;
for degs in [<8,6,8>, <9,7,9>, <10,8,10>, <12,10,12>] do
    dega := degs[1]; degb := degs[2]; degd := degs[3];
    d1, N1 := FitModP(10007, dp, G, dega, degb, degd);
    if d1 eq 0 then continue; end if;
    d2, N2 := FitModP(10009, dp, G, dega, degb, degd);
    if d2 ne d1 then continue; end if;
    v1 := Basis(N1)[1]; v2 := Basis(N2)[1];
    i0 := Min([i : i in [1..Ncols(v1)] | v1[i] ne 0]);
    v1 := v1/v1[i0]; v2 := v2/v2[i0];
    ZN := Integers(10007*10009);
    co := []; okall := true;
    for i in [1..Ncols(v1)] do
        xr := CRT([Integers()!v1[i], Integers()!v2[i]], [10007, 10009]);
        okr, q := RationalReconstruction(ZN!xr);
        if not okr then okall := false; break; end if;
        Append(~co, q);
    end for;
    if not okall then continue; end if;
    depoly := Px ! co[1..degd+1];
    aepoly := Px ! co[degd+2..degd+dega+2];
    bepoly := Px ! co[degd+dega+3..degd+dega+degb+3];
    if depoly eq 0 then continue; end if;
    FF := FieldOfFractions(R2);
    xiF := FF!dp[1]/FF!dp[3]; etaF := FF!dp[2]/(FF!dp[3])^3;
    lhs := Evaluate(depoly, xiF)*FF!E - Evaluate(aepoly, xiF) - Evaluate(bepoly, xiF)*etaF;
    if NormalForm(Numerator(lhs), [G]) eq 0 then
        printf "EXACT VERIFIED: e = (ae + be*y)/de with\n  ae = %o\n  be = %o\n  de = %o\n",
            aepoly, bepoly, depoly;
        PrintFile("data/contact6_m612_efun_abd.txt",
            Sprintf("<%m, %m, %m>", aepoly, bepoly, depoly) : Overwrite:=true);
        found := true; break;
    end if;
end for;
if not found then print "e-fit NOT FOUND up to degree 7"; end if;
print "COLEMAN_PREP2_DONE";
quit;
