
//////////////////////////////////////////////////////////////////////
//  Symbolic reduction of the [4,16] second-halving cover of M_1(8,4)
//  to a surface in A^4, and local analysis at the <0,0> branch.
//
//  P_R = (-R, Y_R) is divisible by 2 iff
//      D(x) := f - c4*(x+R)*q^2  =  ell^2
//  for rational q = x^2+ax+b, ell quadratic.  Eliminating ell exactly:
//  D (deg <= 4) is a square of a quadratic iff
//      E1 = 8*d4^2*d1 - d3*(4*d4*d2 - d3^2) = 0,
//      E0 = 64*d4^3*d0 - (4*d4*d2 - d3^2)^2 = 0,   d4 != 0,
//  together with the squareclass condition d4 = square (so that
//  ell = sqrt(d4)*(x^2 + s x + t) is rational).
//
//  So the reduced cover is the surface
//      Sigma' : E1 = E0 = 0   in  A^4_{R,w,a,b},
//  and rational [4,16] points are rational points of Sigma' with
//  d4 in (Q*)^2 (plus the first-cover tangent condition on (R,w)).
//
//  The 7-adic tower at <0,0> (agent_m18_416_p7_smooth_scan.m) showed a
//  persistent never-smooth branch over (R,w) == (0,0) mod 7 with
//  square-structure signature.  Here we:
//    1. compute and factor E1, E0 (identify boundary factors);
//    2. expand at the origin (R,w,a,b)=(0,0,0,0): lowest-degree forms;
//    3. probe the singular locus of Sigma' through the origin.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();
R4<R,w,a,b> := PolynomialRing(Q, 4, "grevlex");
K4 := FieldOfFractions(R4);
PX<x> := PolynomialRing(K4);

t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
A := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
B := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);
f := x*A*B;
c4 := R + 2 + 4*t;
q := x^2 + a*x + b;
D := f - c4*(x+R)*q^2;
printf "deg D = %o (expect <= 4)\n", Degree(D);

d4 := Coefficient(D,4); d3 := Coefficient(D,3); d2 := Coefficient(D,2);
d1 := Coefficient(D,1); d0 := Coefficient(D,0);

E1 := 8*d4^2*d1 - d3*(4*d4*d2 - d3^2);
E0 := 64*d4^3*d0 - (4*d4*d2 - d3^2)^2;

function ClearPrim(g)
    gn := Numerator(g);
    den := LCM([Denominator(co) : co in Coefficients(gn)]);
    gn := den*gn;
    gg := GCD([Z!co : co in Coefficients(gn)]);
    return R4!(gn/gg);
end function;

E1n := ClearPrim(E1);
E0n := ClearPrim(E0);
d4n := ClearPrim(d4);
printf "E1: total degree %o, terms %o\n", TotalDegree(E1n), #Terms(E1n);
printf "E0: total degree %o, terms %o\n", TotalDegree(E0n), #Terms(E0n);
printf "d4 (cleared): %o\n", d4n;

print "Factorization of E1:";
fac1 := Factorization(E1n);
for ff in fac1 do
    printf "  <deg %o, terms %o, mult %o>", TotalDegree(ff[1]), #Terms(ff[1]), ff[2];
    if #Terms(ff[1]) le 8 then printf "  %o", ff[1]; end if;
    printf "\n";
end for;
print "Factorization of E0:";
fac0 := Factorization(E0n);
for ff in fac0 do
    printf "  <deg %o, terms %o, mult %o>", TotalDegree(ff[1]), #Terms(ff[1]), ff[2];
    if #Terms(ff[1]) le 8 then printf "  %o", ff[1]; end if;
    printf "\n";
end for;

// core factors: strip boundary pieces not vanishing at the origin and
// known artifact loci; keep every factor that vanishes at the origin.
function VanishesAtOrigin(g)
    return Evaluate(g, [0,0,0,0]) eq 0;
end function;
print "Factors of E1 vanishing at origin (the <0,0> branch lives here):";
E1van := [ff[1] : ff in fac1 | VanishesAtOrigin(ff[1])];
for g in E1van do
    printf "  deg %o terms %o\n", TotalDegree(g), #Terms(g);
end for;
print "Factors of E0 vanishing at origin:";
E0van := [ff[1] : ff in fac0 | VanishesAtOrigin(ff[1])];
for g in E0van do
    printf "  deg %o terms %o\n", TotalDegree(g), #Terms(g);
end for;

// lowest-degree homogeneous parts at the origin for the vanishing factors
function LowestForm(g)
    dmin := Min([TotalDegree(m) : m in Terms(g)]);
    lf := R4!0;
    for tm in Terms(g) do
        if TotalDegree(tm) eq dmin then lf +:= tm; end if;
    end for;
    return lf, dmin;
end function;
for g in E1van do
    lf, dmin := LowestForm(g);
    printf "E1 factor deg %o: lowest form degree %o = %o\n",
        TotalDegree(g), dmin, lf;
    printf "  lowest form factorization: %o\n", Factorization(lf);
end for;
for g in E0van do
    lf, dmin := LowestForm(g);
    printf "E0 factor deg %o: lowest form degree %o = %o\n",
        TotalDegree(g), dmin, lf;
    printf "  lowest form factorization: %o\n", Factorization(lf);
end for;

print "DONE";
quit;
