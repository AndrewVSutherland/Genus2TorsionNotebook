SetColumns(0); SetMemoryLimit(10*10^9);
Q := Rationals();
K<r,t> := RationalFunctionField(Q, 2); NR := PolynomialRing(Q, 2);
Px<x> := PolynomialRing(K);
pv := r*(r+t)/2; e := t^2 - 2*pv*t/r; lambda := r/t;
A := r^2 - lambda; B := 2*r*pv - 2*lambda*(pv + r*t) + 2*r*lambda;
C := pv^2 + 2*pv*r^2 - r^4 - r^3*t - r*pv^2/t - lambda*(r^2+e) + 2*lambda*(r*pv + r^2*t - 3*pv*t + r*t^2);
qq := A*x^2 + B*x + C; f := qq*(x^4 + qq);
fc := [Coefficient(f,i) : i in [0..6]]; f1:=fc[2];f2:=fc[3];f3:=fc[4];f4:=fc[5];f5:=fc[6];f6:=fc[7];
l := C; j := f1/(2*l); n := (f2 - j^2)/(2*l);
Sm<mm> := PolynomialRing(K); kap := f6 - mm^2; bN := 2*mm*n - f5;
Eq4 := 3*kap*(f4 - n^2 - 2*mm*j) - bN^2; Eq3 := 27*kap^2*(f3 - 2*mm*l - 2*n*j) + bN^3;
R3<rr,tt,mmm> := PolynomialRing(Q, 3);
ClearToR3 := function(poly)
    cs := Coefficients(poly); D := NR!1;
    for c in cs do if c ne 0 then D := LCM(D, NR!Denominator(c)); end if; end for;
    res := R3!0;
    for i in [1..#cs] do if cs[i] eq 0 then continue; end if;
        Nn := NR!Numerator(cs[i]*(K!D)); mons := Monomials(Nn); cof := Coefficients(Nn);
        for u in [1..#mons] do ex := Exponents(mons[u]); res +:= (R3!cof[u])*rr^ex[1]*tt^ex[2]*mmm^(i-1); end for;
    end for; return res; end function;
P4 := ClearToR3(Eq4); P3 := ClearToR3(Eq3);
g0 := GCD(P4,P3); if TotalDegree(g0) ge 1 then P4:=P4 div g0; P3:=P3 div g0; end if;

// eliminate t -> (r,m) model
Rt := Resultant(P4, P3, tt);
R2<Rr2,Mm2> := PolynomialRing(Q, 2);
DropT := function(poly) res := R2!0; mons := Monomials(poly); cof := Coefficients(poly);
    for u in [1..#mons] do ex := Exponents(mons[u]); res +:= (R2!cof[u])*Rr2^ex[1]*Mm2^ex[3]; end for; return res; end function;
Plane_rm := DropT(Rt);
printf "(r,m) resultant total deg = %o\n", TotalDegree(Plane_rm);
Grm := 0;
for g in Factorization(Plane_rm) do
    printf "  [mult %o deg %o] B(1/3,2/9)? %o\n", g[2], TotalDegree(g[1]), Evaluate(g[1],[1/3,2/9]) eq 0;
    if Evaluate(g[1],[1/3,2/9]) eq 0 and TotalDegree(g[1]) ge 3 then Grm := g[1]; end if;
end for;
if Grm cmpne 0 then
    printf "genuine (r,m) component deg = %o\n", TotalDegree(Grm);
    Crm := Curve(AffineSpace(R2), Grm); PCrm := ProjectiveClosure(Crm);
    printf "genus(r,m model)=%o\n", Genus(PCrm);
    pB := PCrm ! [1/3, 2/9];
    try
        E, mp := EllipticCurve(PCrm, pB); E := MinimalModel(E);
        printf "E=%o cond=%o tors=%o\n", E, Conductor(E), Invariants(TorsionSubgroup(E));
        rk, pf := Rank(E);
        printf "RANK = %o (proven=%o)\n", rk, pf;
        if rk gt 0 then printf "gens=%o\n", [P:P in Generators(E)]; end if;
    catch ee printf "EllipticCurve(r,m) failed: %o\n", ee`Object; end try;
end if;
print "DONE"; quit;
