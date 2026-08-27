// verify_moduli_identities.m -- the symbolic identities and small
// computations quoted without proof in Sections 2, 3 and 4 of the paper:
//
//   (1) Lemma 2.3(b) [Zarhin's halving formula]: on
//         y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2),
//       the Mumford divisor D_0 = (x^2 - s2 x + s4, (s1 s2 - s3)x - s1 s4)
//       (s_i the elementary symmetric functions of a,b,c,d) satisfies
//       2 D_0 = [(0,0) - infty].  Verified generically over Q(rho,sigma,tau)
//       with a = 1 (no loss: the family is invariant under the weighted
//       scaling (a,b,c,d; x,y) -> (la,lb,lc,ld; l^2 x, l^5 y)).
//   (2) Lemma 2.3(d): the four expressions displayed there are EXACTLY the
//       values of the Mumford polynomial q = x^2 - s2 x + s4 at the branch
//       points:  q(-a^2) = a(a+b)(a+c)(a+d), and symmetrically.
//   (3) Lemma 4.5: E_u ~ V_t for t = (u+1)/(u-1), via
//         (x0,y0) |-> ( (x0-u)/(x0+u), 2(t+1) y0 / (x0+u)^2 ).
//       (The map printed in an earlier manuscript draft,
//       (t-1)^4(t+1) y0/(x0+2(t^2-1)), is not correct; the identity below
//       proves the displayed map above.)
//   (4) Proposition 4.4: over Q(s), the displayed x(P'_1), x(P'_2) lie on
//       E_u for the displayed u(s); x(P'_2) is a square (so P'_2 is in the
//       image of the dual 2-isogeny); the displayed x(P'_1 + P'_2) is the
//       x-coordinate of P'_1 + P'_2 (for the right choice of signs of the
//       y-coordinates); (a : b : c : d) = (phi(x(P'_1)) : phi(x(P'_1+P'_2))
//       : 1 : t) with phi the x-map of (3); and -- the conclusion -- the
//       displayed a,b,c,d satisfy the four square conditions of
//       Lemma 2.3(c) identically in Q(s), so [2,2,4,4] c J_s(Q) for every
//       nondegenerate specialization.  One specialization is checked
//       exactly.
//   (5) Corollary 4.3: E_7 : y^2 = x(x+1)(x+49) has Mordell-Weil rank 1.
//   (6) The [2,2,2,8] section: Disc(q) = (ab+ac+ad+bc+bd+cd)^2 - 4abcd,
//       so the doubled-point condition is the quartic K3 surface S of
//       Eq.(eq:2228-k3); the displayed rational curve
//         a = 4t^2(t+1)/(t^2+t+1)^2, b = t/(t+1), c = -1, d = -t
//       satisfies S identically in Q(t); two specializations are checked
//       to give genus-2 curves with [2,2,2,8] embedded in the torsion.
//       Density mechanism of the section's theorem: on the displayed curve
//       d/c = t, and the generic fiber S_t = S cap {d = tc} is a plane
//       quartic of geometric genus 1 over Q(t); the displayed elliptic
//       model A_t : y^2 = x(x-t-1)(x-(t+1)/t) contains (1,1) identically,
//       of infinite order (specialization t=2, hence generically); at
//       sample fibers the Jacobian of S_t is isomorphic to A_t, so each
//       such fiber carries infinitely many rational points and the section
//       sweeps a dense set.
//   (7) Theorem 3.3 (infinitude of S^o(Q)): on the hyperplane
//       H: u-v-a+b = 0 (chart u = a-b+v), the quartic a^4+b^4+c^4-u^4-v^4
//       is irreducible as a polynomial but splits, MODULO the quadric
//       a^2+b^2+c^2-u^2-v^2, as (smooth quadric)*(rank-2 quadric); the
//       rank-2 factor is the plane pair (a+v)(b-v), each plane cutting the
//       quadric surface in a double line inside the degenerate locus
//       {c = 0}, while the smooth factor cuts a smooth genus-1 (2,2)-curve
//       through P_1 whose Jacobian is E : y^2 = x^3 - 21x - 20 itself;
//       E has conductor 288 and Mordell-Weil rank exactly 1, and
//       Q = (-3, 4) generates E(Q) modulo torsion; the displayed map
//       psi : E -> S satisfies both surface equations, the hyperplane
//       equation and the (2,2)-curve equation identically in the function
//       field of E; psi(2Q) = P_1 exactly as points of P^4; and
//       psi(4Q) != P_1 (so psi is nonconstant), with psi(4Q) in S^o.
// Run from this directory:  magma -b verify_moduli_identities.m
SetColumns(0);
SetMemoryLimit(8*10^9);
Q := Rationals();
t0c := Cputime();

// helper: does [2,2,4,4]-type G0 embed in the abelian group with invariant
// factors I (ascending divisibility)?
function Embeds(G0, I)
    if #G0 gt #I then return false; end if;
    J := I[#I-#G0+1..#I];
    return &and[ IsDivisibleBy(J[k], G0[k]) : k in [1..#G0] ];
end function;

// ---- (1) Lemma 2.3(b): the halving formula, generically ----
K3<rho,sigma,tau> := RationalFunctionField(Q, 3);
Pz<z> := PolynomialRing(K3);
s1 := 1 + rho + sigma + tau;
s2 := rho + sigma + tau + rho*sigma + rho*tau + sigma*tau;
s3 := rho*sigma + rho*tau + sigma*tau + rho*sigma*tau;
s4 := rho*sigma*tau;
Jgen := Jacobian(HyperellipticCurve(
    z*(z+1)*(z+rho^2)*(z+sigma^2)*(z+tau^2)));
D0 := Jgen![z^2 - s2*z + s4, (s1*s2 - s3)*z - s1*s4];
assert 2*D0 eq Jgen![z, 0];
printf "(1) Lemma 2.3(b): 2 D_0 = [(0,0) - infty] generically\n";

// ---- (2) Lemma 2.3(d): the values q(-a_i^2) ----
K4<a,b,c,d> := RationalFunctionField(Q, 4);
e2 := a*b+a*c+a*d+b*c+b*d+c*d; e4 := a*b*c*d;
q := func< x | x^2 - e2*x + e4 >;
assert q(-a^2) eq a*(a+b)*(a+c)*(a+d);
assert q(-b^2) eq b*(b+a)*(b+c)*(b+d);
assert q(-c^2) eq c*(c+a)*(c+b)*(c+d);
assert q(-d^2) eq d*(d+a)*(d+b)*(d+c);
printf "(2) Lemma 2.3(d): the four displayed expressions are q(-a_i^2) exactly\n";

// ---- (3) Lemma 4.5: the model isomorphism E_u ~ V_t ----
Ku<u> := FunctionField(Q);
t := (u+1)/(u-1);
R2<x0,y0> := PolynomialRing(Ku, 2);
FF := FieldOfFractions(R2);
fE := x0*(x0+1)*(x0+u^2);
phi := (x0 - u)/(x0 + u);
lhs := (phi^2 - 1)*(phi^2 - t^2);       // the V_t equation at the image
// y_V = 2(t+1) y0/(x0+u)^2: check lhs * (x0+u)^4 = (2(t+1))^2 y0^2-side
assert FF!lhs eq FF!fE * (2*(t+1)/(x0+u)^2)^2;
printf "(3) Lemma 4.5: (x0,y0) -> ((x0-u)/(x0+u), 2(t+1)y0/(x0+u)^2) maps E_u to V_t\n";

// ---- (4) Proposition 4.4 ----
Ks<s> := FunctionField(Q);
us := (1/3) * (s^2-1/2)^-1 * (s^2+1/2)^-1 * (s^2-7/2) * (s^2+3/2);
x1 := (-1/3) * (s^2+1/2)^-2 * (s^2-1/2)^-1 * (s^2+3/2) * (s^2-7/2)^2;
x2 := (1/9) * s^2 * (s^2-1/2)^-2 * (s^2-7/2)^2;
x12 := (-1/3) * (s^2-1/2)^-1 * (s^2+3/2)
       * (s^3-s^2+(5/2)*s-3/2)^-2 * (s^3-3*s^2+(5/2)*s+3/2)^2;
fEs := func< x | x*(x+1)*(x+us^2) >;
ok1, y1 := IsSquare(fEs(x1)); assert ok1;       // P'_1 in E_u(Q(s))
ok2, y2 := IsSquare(fEs(x2)); assert ok2;       // P'_2 in E_u(Q(s))
assert IsSquare(x2);                            // P'_2 in image of dual isogeny
E := EllipticCurve([0, 1+us^2, 0, us^2, 0]);    // y^2 = x(x+1)(x+u^2)
S12 := E![x1, y1] + E![x2, -y2];                // sign choice for P'_2
assert S12[1] eq x12;
av := (s - 2)*(s^2 - s + 1/2)*(s^2 - 3/2)*(s^4 + s^2 + 9/4);
bv := (s^2 - 2*s + 1/2)*(s^2 - 2*s + 9/2)*(s^2 - 3/2)*(s^2 + 1);
cv := -(s - 2)*(2*s^2 - 2*s + 1)*(s^4 + s^2 + 9/4);
dv := (s - 2)*(2*s^2 - 2*s + 1)*(2*s^2 - 3)*(s^2 + 1);
ts := (us+1)/(us-1);
phis := func< x | (x - us)/(x + us) >;
assert av/cv eq phis(x1) and bv/cv eq phis(x12) and dv/cv eq ts;
A := av^2; B := bv^2; Cc := cv^2; Dd := dv^2;
assert IsSquare((Cc-A)*(Cc-B)) and IsSquare((Cc-A)*(Dd-A))
   and IsSquare((Cc-B)*(Dd-B)) and IsSquare((Dd-A)*(Dd-B));
printf "(4) Prop 4.4: all displayed identities hold over Q(s); Lemma 2.3(c) conditions are identically squares\n";
// one exact specialization
Px<x> := PolynomialRing(Q);
sv := Q!3;
abcd := [ Evaluate(f0, sv) : f0 in [av, bv, cv, dv] ];
fs := x * &*[ x + w^2 : w in abcd ];
assert Discriminant(fs) ne 0;
ds := LCM([ Denominator(co) : co in Coefficients(fs) ]);
Is := Invariants(AbelianGroup(TorsionSubgroup(Jacobian(
    HyperellipticCurve(ds^2*fs)))));
assert Embeds([2,2,4,4], Is);
printf "    s = 3 fiber: torsion %o contains [2,2,4,4]\n", Is;

// ---- (5) Corollary 4.3: rank E_7 = 1 ----
E7 := EllipticCurve([0, 50, 0, 49, 0]);          // y^2 = x(x+1)(x+49)
r, exact := Rank(E7);
assert exact and r eq 1;
printf "(5) Cor 4.3: E_7 has Mordell-Weil rank exactly 1\n";

// ---- (6) Theorem 4.7: the [2,2,2,8] K3 surface ----
assert e2^2 - 4*e4 eq (a*b+a*c+a*d+b*c+b*d+c*d)^2 - 4*a*b*c*d;  // Disc(q)
Kt<tt> := FunctionField(Q);
at := 4*tt^2*(tt+1)/(tt^2+tt+1)^2;
bt := tt/(tt+1); ct := Kt!-1; dt := -tt;
assert (at*bt+at*ct+at*dt+bt*ct+bt*dt+ct*dt)^2 eq 4*at*bt*ct*dt;
assert dt/ct eq tt;                       // the fibration parameter on S
for tv in [ Q!2, Q!5 ] do
    q4 := [ Evaluate(w, tv) : w in [at, bt, ct, dt] ];
    ft := x * &*[ x + w^2 : w in q4 ];
    assert Discriminant(ft) ne 0;
    dt2 := LCM([ Denominator(co) : co in Coefficients(ft) ]);
    It := Invariants(AbelianGroup(TorsionSubgroup(Jacobian(
        HyperellipticCurve(dt2^2*ft)))));
    assert Embeds([2,2,2,8], It);
    printf "(6) t = %o fiber on S: torsion %o contains [2,2,2,8]\n", tv, It;
end for;
// the elliptic fibration and the density mechanism of the section's theorem
P2t<A2,B2,C2> := ProjectiveSpace(Kt, 2);
Ft := (A2*B2 + A2*C2 + A2*(tt*C2) + B2*C2 + B2*(tt*C2) + C2*(tt*C2))^2
      - 4*A2*B2*C2*(tt*C2);              // S cap {d = t c}, d eliminated
SFib := Curve(P2t, Numerator(Ft));
assert GeometricGenus(SFib) eq 1;         // generic fiber is genus 1
printf "(6) generic fiber S_t (plane quartic, d = tc) has geometric genus 1 over Q(t)\n";
// A_t : y^2 = x(x-t-1)(x-(t+1)/t) contains (1,1) identically
assert (Kt!1) eq 1*(1 - tt - 1)*(1 - (tt+1)/tt);
// (1,1) has infinite order at t = 2 (hence generically)
A2E := EllipticCurve([0, -(2+1) - (2+1)/2, 0, (2+1)*(2+1)/2, 0]);
assert Order(A2E![1,1]) eq 0;
// sample fibers: Jac(S_t) is isomorphic to A_t, so rank >= 1 fiberwise
for tv in [ Q!2, Q!5, Q!-3/2 ] do
    av := Evaluate(at, tv); bv := Evaluate(bt, tv); cv := Q!-1;
    P2Q<AA,BB,CC> := ProjectiveSpace(Q, 2);
    FQ2 := (AA*BB + AA*CC + AA*(tv*CC) + BB*CC + BB*(tv*CC) + CC*(tv*CC))^2
           - 4*AA*BB*CC*(tv*CC);
    Cfib := Curve(P2Q, FQ2);
    pfib := Cfib![av, bv, cv];
    Efib := EllipticCurve(Cfib, pfib);
    Atv := EllipticCurve([0, -(tv+1) - (tv+1)/tv, 0, (tv+1)^2/tv, 0]);
    assert IsIsomorphic(Efib, Atv);
    assert Order(Atv![1,1]) eq 0;
    printf "(6) t = %o fiber: Jac(S_t) = A_t, section of infinite order\n", tv;
end for;

// ---- (7) Theorem 3.3: infinitude of S^o(Q) via E : y^2 = x^3 - 21x - 20 ----
// (a) the hyperplane section H cap S, on H : u - v - a + b = 0 (u = a-b+v)
R4H<Ah,Bh,Ch,Vh> := PolynomialRing(Q, 4);
Uh := Ah - Bh + Vh;
Q0H := Ah^2 + Bh^2 + Ch^2 - Uh^2 - Vh^2;          // the restricted quadric
F4H := Ah^4 + Bh^4 + Ch^4 - Uh^4 - Vh^4;          // the restricted quartic
assert IsIrreducible(F4H);       // NOT a factorization of polynomials ...
qSm := 4*Ah^2 + 4*Bh^2 + Ch^2 - 4*Ah*Vh + 4*Bh*Vh - 4*Vh^2;
qSing := (Ah + Vh)*(Bh - Vh);
// ... but modulo the quadric the quartic splits into the two quadrics:
assert F4H - qSm*qSing eq (Ch^2 - 3*Ah*Bh + 3*Ah*Vh - 3*Bh*Vh + 3*Vh^2)*Q0H;
GramRank := func< qq | Rank(Matrix(Q, 4, 4,
    [ Q!Derivative(Derivative(qq, i), j) : j in [1..4], i in [1..4] ])) >;
assert GramRank(qSing) eq 2 and GramRank(qSm) eq 4;   // singular resp. smooth
// each plane of qSing cuts the quadric surface in a double line inside the
// degenerate locus {c = 0}:
assert Evaluate(Q0H, [-Vh, Bh, Ch, Vh]) eq Ch^2;      // the plane a + v = 0
assert Evaluate(Q0H, [Ah, Vh, Ch, Vh]) eq Ch^2;       // the plane b - v = 0
// the smooth factor cuts a smooth genus-1 (2,2)-curve through P_1:
P3H<a3,b3,c3,v3> := ProjectiveSpace(Q, 3);
CsmH := Curve(Scheme(P3H, [ Evaluate(g0, [a3,b3,c3,v3]) : g0 in [Q0H, qSm] ]));
assert IsNonsingular(CsmH) and GeometricGenus(CsmH) eq 1;
p1h := CsmH![120, 143, 266, 241];         // P_1, with u = 120-143+241 = 218
printf "(7) Thm 3.3: on H, quartic = (smooth quadric)*(rank-2 quadric) mod the quadric;\n";
printf "    plane pair (a+v)(b-v) cuts double lines in {c=0}; smooth part cuts a genus-1 (2,2)-curve through P_1\n";
// (b) the elliptic curve E = 288b1 and its Mordell-Weil group
E288 := EllipticCurve([0, 0, 0, -21, -20]);           // y^2 = x^3 - 21x - 20
assert Conductor(E288) eq 288;
r288, ex288 := Rank(E288);
assert ex288 and r288 eq 1;
MW288, mw288 := MordellWeilGroup(E288);
assert Invariants(MW288) eq [2, 2, 0];                // E(Q) = Z/2 x Z/2 x Z
Qg := E288![-3, 4];
gfree := [ mw288(MW288.i) : i in [1..Ngens(MW288)] | Order(MW288.i) eq 0 ][1];
assert Order(Qg - gfree) gt 0 or Order(Qg + gfree) gt 0;  // Q = +-gfree + torsion
// the Jacobian of the genus-1 curve is E itself (whence the parametrization):
assert IsIsomorphic(EllipticCurve(CsmH, p1h), E288);
printf "    E has conductor 288 and rank exactly 1; Q = (-3,4) generates E(Q)/torsion;\n";
printf "    the genus-1 (2,2)-curve has Jacobian isomorphic to E\n";
// (c) the map psi : E -> S, as identities in the function field of E.
//     NOTE: an earlier manuscript draft displayed psi as
//       [3t(t-1) : 1-3t : 2w(3t^2+1) : 6t^2-3t+1 : 3t^2-3t+2],
//     which sends 2Q to [-120:-143:266:241:218] -- the image of P_1 under
//     the G-action (negate a,b; swap u,v), not P_1 itself.  The corrected
//     psi below (first two coordinates negated, last two swapped) hits P_1
//     on the nose.
Psi := function(x0v, y0v)
    wv := (3*x0v + y0v - 6)/(y0v + 18);
    tv2 := -(4*x0v^2 - 7*x0v + 12*y0v + 16)/((x0v - 8)*(3*x0v + 2*y0v + 12));
    return [ 3*tv2*(1 - tv2), 3*tv2 - 1, 2*wv*(3*tv2^2 + 1),
             3*tv2^2 - 3*tv2 + 2, 6*tv2^2 - 3*tv2 + 1 ];
end function;
FF288<xE, yE> := FunctionField(E288);
cg := Psi(xE, yE);
assert cg[1]^2 + cg[2]^2 + cg[3]^2 eq cg[4]^2 + cg[5]^2;   // psi lands on S,
assert cg[1]^4 + cg[2]^4 + cg[3]^4 eq cg[4]^4 + cg[5]^4;
assert cg[4] - cg[5] - cg[1] + cg[2] eq 0;                  // on H,
assert 4*cg[1]^2 + 4*cg[2]^2 + cg[3]^2
       - 4*cg[1]*cg[5] + 4*cg[2]*cg[5] - 4*cg[5]^2 eq 0;    // and on the (2,2)-curve
printf "    psi maps E onto the genus-1 curve on S (identities in Q(E))\n";
// (d) psi(2Q) = P_1 exactly, and psi is nonconstant
P4S := ProjectiveSpace(Q, 4);
P1pt := P4S![120, 143, 266, 218, 241];
tq := 2*Qg;                                       // 2Q = (105/16, -715/64)
im2 := Psi(tq[1], tq[2]);
assert P4S!im2 eq P1pt;
fq := 4*Qg;
im4 := Psi(fq[1], fq[2]);
assert P4S!im4 ne P1pt;                           // psi nonconstant
sq2 := [ z^2 : z in im2 ]; sq4 := [ z^2 : z in im4 ];
assert 0 notin sq2 and #Seqset(sq2) eq 5;         // psi(2Q) in S^o
assert 0 notin sq4 and #Seqset(sq4) eq 5;         // psi(4Q) in S^o too
printf "    psi(2Q) = P_1 = [120:143:266:218:241] EXACTLY; psi(4Q) != P_1 and lies in S^o\n";

printf "ALL MODULI IDENTITIES VERIFIED (%.1o s)\n", Cputime() - t0c;
quit;
