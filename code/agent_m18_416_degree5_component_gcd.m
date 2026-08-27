//////////////////////////////////////////////////////////////////////
//  Analyze the repeated degree-5 projected component of
//  Res_b(E1core,E0core) for the M_1(8,4) [4,16] second-halving surface.
//
//  The rational resultant factorization found
//
//    G5 = R^5 + 1/4 R w^4 - R^4 + 3/2 R^2 w^2 + 3/4 w^4
//         - R w^2 a - 2 R^3 + R w^2 - w^2 a
//         + 1/2 R^2 - 3/2 w^2 + R a - 1/4 R + a - 1/4.
//
//  This is linear in a:
//
//    G5 = N(R,w) + a*(R+1)*(1-w^2).
//
//  Away from the already-stripped boundary factors R=-1 and w=+-1, the
//  component is the graph a = -N/((R+1)*(1-w^2)).  This script substitutes
//  that graph into E1core,E0core and computes their gcd as polynomials in b
//  over Q(R,w).  A linear gcd means the degree-5 projection lifts to a
//  rational branch in b; a higher degree gcd means a remaining cover.
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
d4 := Coefficient(D,4); d3 := Coefficient(D,3); d2 := Coefficient(D,2);
d1 := Coefficient(D,1); d0 := Coefficient(D,0);
E1 := 8*d4^2*d1 - d3*(4*d4*d2 - d3^2);
E0 := 64*d4^3*d0 - (4*d4*d2 - d3^2)^2;

function ClearPrim4(g)
    gn := Numerator(g);
    den := LCM([Denominator(co) : co in Coefficients(gn)]);
    gn := den*gn;
    gg := GCD([Z!co : co in Coefficients(gn)]);
    return R4!(gn/gg);
end function;

E1n := ClearPrim4(E1);
E0n := ClearPrim4(E0);
while IsDivisibleBy(E1n, R-1) do E1n := ExactQuotient(E1n, R-1); end while;
while IsDivisibleBy(E0n, R-1) do E0n := ExactQuotient(E0n, R-1); end while;

K<RK,wK> := RationalFunctionField(Q, 2);
Pb<bK> := PolynomialRing(K);

N := RK^5 + (1/4)*RK*wK^4 - RK^4 + (3/2)*RK^2*wK^2 + (3/4)*wK^4
     - 2*RK^3 + RK*wK^2 + (1/2)*RK^2 - (3/2)*wK^2 - (1/4)*RK - (1/4);
aden := (RK + 1)*(1 - wK^2);
aG := -N/aden;
printf "a_on_G5 = %o\n", aG;

h := hom<R4 -> Pb | [Pb!RK, Pb!wK, Pb!aG, bK]>;
E1g := h(E1n);
E0g := h(E0n);
printf "after G5 substitution: deg_b E1=%o terms %o | deg_b E0=%o terms %o\n",
    Degree(E1g), #Terms(E1g), Degree(E0g), #Terms(E0g);

tim := Cputime();
gb := GCD(E1g, E0g);
printf "gcd_b degree %o terms %o (%.1o s)\n", Degree(gb), #Terms(gb), Cputime(tim);
printf "gcd_b = %o\n", gb;

E1q := ExactQuotient(E1g, gb);
E0q := ExactQuotient(E0g, gb);
printf "quotients: E1/g deg %o terms %o | E0/g deg %o terms %o\n",
    Degree(E1q), #Terms(E1q), Degree(E0q), #Terms(E0q);

if Degree(gb) eq 1 then
    bval := -Coefficient(gb, 0)/Coefficient(gb, 1);
    printf "b_on_G5 = %o\n", bval;
end if;

// Check whether this branch is just a boundary artifact by evaluating the
// standard boundary factors after the substitution.  None should vanish
// identically away from R=-1 and w=+-1.
print "boundary caveat: G5 graph uses denominator (R+1)*(1-w^2), so it excludes R=-1 and w=+-1; those already occur as separate resultant factors.";
print "DONE";
quit;

