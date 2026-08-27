
//////////////////////////////////////////////////////////////////////
//  Derive the d=0 "Z/24 with split q3 through x=0" family Phi(r,t)=0.
//
//  d=0 slice: p = r(r+t)/2, Q = x^2, f = q*(x^4+q), f(0)=C^2 (C=q(0)).
//  Ansatz: 3-torsion contact with q3 = x*(x-beta) (one root at x=0):
//     f = h3^2 + kappa*q3^3,   q3^3 = x^6 -3b x^5 +3b^2 x^4 - b^3 x^3
//                                     (degrees 3..6 only!).
//  h3 = m x^3 + n x^2 + j x + l.  Matching x^0,x^1,x^2 (q3^3 has no such
//  terms) gives l,j,n RATIONALLY from f (l=C since f0=C^2):
//     l^2=f0 -> l=C ;  f1=2 j l ;  f2=j^2+2 n l.
//  x^6: kappa=f6-m^2 ; x^5: f5+3 kappa b = 2 m n -> b ;
//  x^4,x^3: two equations in m alone -> Resultant_m = Phi(r,t).
//
//  Phi(r,t)=0 is (a component of) the d=0 Z/24 locus.  We factor it,
//  find its genus/rational parametrization, and test 2-rank along it.
//
//  Usage: magma -b agent_a2_24_d0_derive.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB := 6; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals();

K<r,t> := RationalFunctionField(Q, 2);
Px<x> := PolynomialRing(K);

// A(8) f on the d=0 slice (p = r(r+t)/2)
p := r*(r+t)/2;
e := t^2 - 2*p*t/r; d := e + 2*p - r^2; lambda := r/t;
A := r^2 - lambda;
B := 2*r*p - 2*lambda*(p + r*t) + 2*r*lambda;
C := p^2 + 2*p*r^2 - r^4 - r^3*t - r*p^2/t
     - lambda*(r^2 + e) + 2*lambda*(r*p + r^2*t - 3*p*t + r*t^2);
printf "check d = %o (should be 0)\n", d;
qq := A*x^2 + B*x + C;
Qpoly := x^2 + d;                // = x^2 since d=0
f := qq*(Qpoly^2 + qq);
fc := [Coefficient(f, i) : i in [0..6]];   // f0..f6
printf "check f0 - C^2 = %o (should be 0)\n", fc[1] - C^2;

f0:=fc[1]; f1:=fc[2]; f2:=fc[3]; f3:=fc[4]; f4:=fc[5]; f5:=fc[6]; f6:=fc[7];

// h3 low coeffs (rational in r,t)
l := C;                              // l^2 = f0
j := f1/(2*l);
n := (f2 - j^2)/(2*l);

// eliminate m: kappa = f6-m^2, beta from x^5, then x^4 & x^3 eqns in m
Km<m> := PolynomialRing(K);
kap := f6 - m^2;
// x^5: f5 + 3*kap*beta = 2*m*n  => beta = (2 m n - f5)/(3 kap)
betaNum := 2*m*n - f5;
beta3k := 3*kap;                     // beta = betaNum/beta3k
// x^4:  f4 - 3*kap*beta^2 = n^2 + 2*m*j ; beta=betaNum/(3 kap), so
//   f4 - 3*kap*betaNum^2/(9 kap^2) = n^2+2mj => f4 - betaNum^2/(3 kap) = n^2+2mj
//   times 3 kap:  3 kap (f4 - n^2 - 2 m j) - betaNum^2 = 0   (kap divided out)
Eq4 := 3*kap*(f4 - n^2 - 2*m*j) - betaNum^2;
// x^3:  f3 + kap*beta^3 = 2ml+2nj ; kap*betaNum^3/(27 kap^3) = betaNum^3/(27 kap^2)
//   times 27 kap^2:  27 kap^2 (f3 - 2 m l - 2 n j) + betaNum^3 = 0   (kap divided out)
Eq3 := 27*kap^2*(f3 - 2*m*l - 2*n*j) + betaNum^3;

printf "deg_m Eq4=%o  deg_m Eq3=%o\n", Degree(Eq4), Degree(Eq3);
Res := Resultant(Eq4, Eq3);          // in K = Q(r,t)
PhiK := Res;
// clear to a bivariate polynomial
R2<rr,tt> := PolynomialRing(Q, 2);
num := Numerator(PhiK);
// map K-element numerator (a poly in r,t) into R2
function ToBiv(el)
    N := Numerator(el);
    // N is in the underlying poly ring of K; coerce monomial-wise
    mons := Monomials(N); cofs := Coefficients(N);
    res := R2!0;
    for i in [1..#mons] do
        ex := Exponents(mons[i]);
        res +:= (R2!cofs[i]) * rr^ex[1] * tt^ex[2];
    end for;
    return res;
end function;
Phi := ToBiv(PhiK);
printf "Phi(r,t): total degree %o\n", TotalDegree(Phi);
fac := Factorization(Phi);
printf "Phi factors (%o):\n", #fac;
for g in fac do
    printf "  [mult %o, deg %o]  %o\n", g[2], TotalDegree(g[1]), g[1];
end for;

// analyze the non-degenerate components (deg >= 2)
print "==== component genus / parametrization ====";
for g in fac do
    G := g[1];
    if TotalDegree(G) lt 2 then continue; end if;
    printf "-- component deg %o --\n", TotalDegree(G);
    // does curve B (1/3,-1) lie on it?
    printf "   B=(1/3,-1) on it? %o\n", Evaluate(G, [1/3, -1]) eq 0;
    PC := ProjectiveClosure(Curve(AffineSpace(R2), G));
    try
        gg := Genus(PC);
        printf "   genus = %o\n", gg;
        if gg eq 0 then
            bp, prm := IsRational(PC);   // rational parametrization
            printf "   IsRational = %o\n", bp;
        end if;
    catch ee
        printf "   genus/param error: %o\n", ee`Object;
    end try;
end for;
print "DONE";
quit;
