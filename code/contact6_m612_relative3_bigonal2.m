//////////////////////////////////////////////////////////////////////
// Bigonal construction, step 2: the dual tower and the Prym rank.
//   C1' : u^2 = N(x) = -108 (x^2-6x-36)^2 / (x^3 (x+24))
//         ~ conic v^2 = -3 x (x+24), rational point (0,0);
//         parametrize x = -72/(t^2+3), v = -72 t/(t^2+3),
//         u = 6 (x^2-6x-36) v / (x^2 (x+24)).
//   C2' : s^2 = 2*a/d + 2*u  pulled through the parametrization
//         -> s^2 = f2'(t), hyperelliptic; expected genus 2;
//         J(C2') isogenous to Prym(E8/E4)^dual  (Pantazis) => same rank.
// Verification: Prym trace from function fields (mp8/mp4 over GF(p))
// vs trace of J(C2') at p = 5,7,11,13,17.  Then RankBounds(J(C2')).
// Usage: magma -b code/contact6_m612_relative3_bigonal2.m
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(1); SetMemoryLimit(6*10^9);
Q := Rationals();
Kt<t> := FunctionField(Q);
Px<x> := PolynomialRing(Q);

abd := eval Read("data/contact6_m612_bigonal_abd.txt");
a := abd[1]; b := abd[2]; d := abd[3];
fH := (Px.1^6 - 18*Px.1^5 + 162*Px.1^4 - 972*Px.1^3 - 2916*Px.1^2 + 34992*Px.1 + 104976)/4;
// sanity: N = (a^2-b^2*fH)/d^2 = -108 (x^2-6x-36)^2/(x^3(x+24))
Nx := (a^2 - b^2*fH)/d^2;
KX<X> := FunctionField(Q);
NX := Evaluate(Nx, X);
assert NX eq -108*(X^2-6*X-36)^2/(X^3*(X+24));
print "N verified";

// parametrization
xt := -72/(t^2+3);
vt := -72*t/(t^2+3);
assert vt^2 eq -3*xt*(xt+24);
ut := 6*(xt^2-6*xt-36)*vt/(xt^2*(xt+24));
assert ut^2 eq Evaluate(Nx, xt);
print "conic parametrization verified; u^2 = N(x(t)) holds";

// s^2 = 2 a/d + 2 u  at x = x(t)
At := Evaluate(a/d, xt);
S := 2*At + 2*ut;
printf "S(t) = %o\n", S;
num := Numerator(S); den := Denominator(S);
// hyperelliptic form: (s*den)^2 = num*den
f2 := num*den;
// strip square factors for the model
fac := Factorization(f2);
sf := Q!LeadingCoefficient(f2);
sqpart := Parent(f2)!1;
for g in fac do
    if IsOdd(g[2]) then sf := sf*g[1]^1; end if;   // multiply symbolically below
end for;
// rebuild squarefree polynomial exactly: product of odd-multiplicity factors * lead
f2sf := LeadingCoefficient(f2) * &*[Parent(f2) | g[1]^(g[2] mod 2) : g in fac];
printf "f2' (squarefree model rhs) = %o\n", f2sf;
printf "degree %o, factorization %o\n", Degree(f2sf), [<Degree(g[1]),g[2]> : g in Factorization(f2sf)];
C2p := HyperellipticCurve(f2sf);
printf "C2' genus = %o\n", Genus(C2p);

// trace cross-check at small p: Prym trace (E8 minus E4) vs J(C2') trace
print "== L-factor cross-check ==";
K<e> := FunctionField(Q); Kz<zz> := PolynomialRing(K);
mp4Q := Kz ! eval Read("data/contact6_m612_E4_mp4Q.txt");
mp8 := Kz ! [ IsEven(i) select Coefficient(mp4Q, i div 2) else 0 : i in [0..8] ];
for p in [5,7,11,13,17,19] do
    ok := true;
    // function fields over GF(p)
    Kp<ep> := FunctionField(GF(p)); Kpz := PolynomialRing(Kp);
    red := function(c)
        nn := Numerator(c); dd := Denominator(c);
        np := &+[Kp | (GF(p)!Numerator(Coefficient(nn,i)))/(GF(p)!Denominator(Coefficient(nn,i)))*ep^i : i in [0..Degree(nn)]];
        dp2 := &+[Kp | (GF(p)!Numerator(Coefficient(dd,i)))/(GF(p)!Denominator(Coefficient(dd,i)))*ep^i : i in [0..Degree(dd)]];
        return np/dp2;
    end function;
    c4 := [red(Coefficient(mp4Q,i)) : i in [0..4]];
    c8 := [red(Coefficient(mp8,i)) : i in [0..8]];
    if &or[cc eq 0 : cc in [c4[1], c8[1]]] then continue; end if;
    F4p := ext<Kp | Kpz!c4>; F8p := ext<Kp | Kpz!c8>;
    n4 := NumberOfPlacesOfDegreeOneECF(F4p, 1);
    n8 := NumberOfPlacesOfDegreeOneECF(F8p, 1);
    t4 := p + 1 - n4; t8 := p + 1 - n8;
    // C2' over GF(p)
    f2p := PolynomialRing(GF(p)) ! f2sf;
    if not IsSquarefree(f2p) then printf "p=%o: bad reduction for C2', skip\n", p; continue; end if;
    C2pp := HyperellipticCurve(f2p);
    n2 := #Points(C2pp);
    t2 := p + 1 - n2;
    printf "p=%o : trace(E8)-trace(E4) = %o - %o = %o ; trace(J C2') = %o  %o\n",
        p, t8, t4, t8-t4, t2, (t8-t4 eq t2) select "MATCH" else "**MISMATCH**";
end for;

// rank of the Prym
J2 := Jacobian(C2p);
rlo, rhi := RankBounds(J2);
printf "PRYM RANK BOUNDS: %o <= rank <= %o\n", rlo, rhi;
print "BIGONAL2_DONE";
quit;
