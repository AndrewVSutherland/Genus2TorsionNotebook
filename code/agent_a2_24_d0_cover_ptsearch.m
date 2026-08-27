SetColumns(0); SetMemoryLimit(9*10^9);
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
Rt := Resultant(P4, P3, tt);
R2<Rr2,Mm2> := PolynomialRing(Q, 2);
DropT := function(poly) res := R2!0; mons := Monomials(poly); cof := Coefficients(poly);
    for u in [1..#mons] do ex := Exponents(mons[u]); res +:= (R2!cof[u])*Rr2^ex[1]*Mm2^ex[3]; end for; return res; end function;
Plane := DropT(Rt);
Grm := 0;
for g in Factorization(Plane) do if Evaluate(g[1],[1/3,2/9]) eq 0 and TotalDegree(g[1]) ge 3 then Grm := g[1]; end if; end for;
printf "Grm deg=%o  B on it=%o\n", TotalDegree(Grm), Evaluate(Grm,[1/3,2/9]) eq 0;

// point search: for small-height m0, solve Grm(r, m0)=0 for rational r
Pr<xr> := PolynomialRing(Q);
HR := function(HH) v := []; for den in [1..HH] do for num in [-HH..HH] do
    if GCD(num,den) eq 1 then Append(~v, Q!num/den); end if; end for; end for; return Sort(Setseq(Seqset(v))); end function;
Hm := 40;
pts := {};
for m0 in HR(Hm) do
    gpoly := Evaluate(Grm, [xr, m0]);   // univariate in r
    if gpoly eq 0 then continue; end if;
    for rt in Roots(gpoly) do Include(~pts, <rt[1], m0>); end for;
end for;
// also scan r0
for r0 in HR(Hm) do
    gpoly := Evaluate(Grm, [r0, xr]);
    if gpoly eq 0 then continue; end if;
    for rt in Roots(gpoly) do Include(~pts, <r0, rt[1]>); end for;
end for;
printf "rational points on Grm with height<=%o : %o\n", Hm, #pts;
for p in pts do printf "  (r,m)=(%o, %o)\n", p[1], p[2]; end for;
print "DONE"; quit;
