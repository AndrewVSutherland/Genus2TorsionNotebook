
//////////////////////////////////////////////////////////////////////
//  Cross-product analysis of the second-stage conditions S_A, S_B for
//  P_R halving in M_1(8,4), on the double cover g^2 = G(R,m).
//
//  From agent_m18_416_sstage_symbolic.m:
//    V_Apm = 2*alpha_A +- 2*g*h_A,   V_Apm+ * V_Apm- = d_A   (exactly),
//    V_Bpm = 2*alpha_B +- 2*g*h_B,   product = d_B           (exactly),
//    h_A = P3^3/(m^3*R*(R+1)),  h_B = P3^2/(m^2*R*(R+1)),
//    P3 = R^3 - R + m^2/2 = (m^2 - K)/2.
//  S_A <=> V_A+ or V_A- is a square; same for B.  Both S_A and S_B
//  hold iff some pair (V_A*, V_B*) are both squares, iff their product
//  is a square AND one of them is.
//
//  This script computes the four cross-products
//    V_A+ * V_B+ = 4*(E+ + g*Y+),   V_A+ * V_B- = 4*(E- + g*Y-),
//    E+- = alpha_A*alpha_B +- G*h_A*h_B,  Y+- = alpha_B*h_A +- alpha_A*h_B,
//  factors E,Y over Q(R,m), and tests squareness of all four products
//  (and small-class twists) in the function field Q(R,m)(g).
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

K2<R,m> := RationalFunctionField(Q, 2);
PX<x> := PolynomialRing(K2);

K := -2*R*(R^2-1);
w := (m^2 + K)/(m^2 - K);
W := w^2;
t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
c4 := R + 2 + 4*t;
Apol := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
Bpol := c4*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);
XR := -c4*R;
At := PX![K2!co : co in Coefficients(c4^2*Evaluate(Apol, x/c4))];
Bt := PX![K2!co : co in Coefficients(c4*Evaluate(Bpol, x/c4))];
G := 2*(R^2-1)*(R*(2*R+1) - W*(R+2));

alphaA := XR + Coefficient(At,1)/2;
alphaB := XR + Coefficient(Bt,1)/2;
dA := Discriminant(At);
dB := Discriminant(Bt);
NA := Evaluate(At, XR);
NB := Evaluate(Bt, XR);
okA, hA := IsSquare(NA/G); assert okA;
okB, hB := IsSquare(NB/G); assert okB;

// exact product identities
assert (2*alphaA)^2 - 4*G*hA^2 eq dA;
assert (2*alphaB)^2 - 4*G*hB^2 eq dB;
print "exact identities V_A+ V_A- = d_A, V_B+ V_B- = d_B verified";

// cross data
Ep := alphaA*alphaB + G*hA*hB;
Em := alphaA*alphaB - G*hA*hB;
Yp := alphaB*hA + alphaA*hB;
Ym := alphaB*hA - alphaA*hB;

procedure ShowFactored(name, val)
    if val eq 0 then printf "%o = 0  (!!)\n", name; return; end if;
    num := Numerator(val); den := Denominator(val);
    printf "%o:\n  num: %o\n  den: %o\n", name, Factorization(num), Factorization(den);
end procedure;

print "\n== cross-product components ==";
ShowFactored("E+ (g-even part of V_A+V_B+/4)", Ep);
ShowFactored("E- (g-even part of V_A+V_B-/4)", Em);
ShowFactored("Y+ (g-odd part of V_A+V_B+/4)", Yp);
ShowFactored("Y- (g-odd part of V_A+V_B-/4)", Ym);

// ---- norm obstruction to function-level identities ----
// For X + Y*g in Q(R,m)(g), squareness requires N = X^2 - G*Y^2 to be a
// square in Q(R,m).  N(V_A+-) = d_A, N(V_B+-) = d_B, N(V_A V_B) = d_A d_B,
// and twisting by a base function c multiplies the norm by c^2 (same
// class).  So a function-level square identity exists iff the norm class
// is trivial.  Check the norm classes:
print "\n== norm classes (trivial <=> function-level identity possible) ==";
for pr in [<"d_A", dA>, <"d_B", dB>, <"d_A*d_B", dA*dB>] do
    val := pr[2];
    issq, _ := IsSquare(val);
    printf "  %o is a square in Q(R,m): %o\n", pr[1], issq;
end for;
print "conclusion: if all false, NO base twist makes any V or V*V a";
print "function-field square: the second-stage wall is arithmetic";
print "(point-dependent), not a constant-class identity.";

// ---- reduced 2-cover curves per R-fiber ----
// Eliminating g from {g^2 = G, y^2 = 2 alpha + 2 g h} using 4 h^2 G =
// 4 alpha^2 - d gives the plane curve
//     y^4 - 4*alpha(m)*y^2 + d(m) = 0
// over the m-line: the S-condition 2-cover of the fiber E_R: g^2 = G(m).
print "\n== reduced 2-cover curves (per R-fiber, over the m-line) ==";
print "S_A-cover:  y^4 - 4*alpha_A(m)*y^2 + d_A(m) = 0";
print "S_B-cover:  y^4 - 4*alpha_B(m)*y^2 + d_B(m) = 0";
printf "deg_m: alpha_A num %o / den %o, d_A num %o / den %o\n",
    Degree(Numerator(alphaA), 2), Degree(Denominator(alphaA), 2),
    Degree(Numerator(dA), 2), Degree(Denominator(dA), 2);
printf "deg_m: alpha_B num %o / den %o, d_B num %o / den %o\n",
    Degree(Numerator(alphaB), 2), Degree(Denominator(alphaB), 2),
    Degree(Numerator(dB), 2), Degree(Denominator(dB), 2);
print "Full halving <=> a rational point on the fiber product of the two";
print "covers over E_R -- per-R Selmer/local analysis is the next tool.";
print "DONE";
quit;
