//////////////////////////////////////////////////////////////////////
// Geometry of the first-radicand subcover of the orbit-12 halving
// cover over the complete rational Clebsch--Klein chart.
//////////////////////////////////////////////////////////////////////

Q := Rationals();
P<t,m> := PolynomialRing(Q,2);

R1 := 1+t*(t+2)*m;
R2 := t*m*(m-t-2);
R3 := -1+m+t*(t+1)*m^2;
R4 := 1+t-m-t*m^2;
R5 := -(1+t)*(1+t*m^2);
RR := [R1,R2,R3,R4,R5];

assert &+RR eq 0;
assert &+[r^3 : r in RR] eq 0;

// For the fixed marked pair {R1^2,R2^2}, the common exact
// Stoll--Zarhin squareclass is
//
//   G0 = -(R1^2-R3^2)(R1^2-R4^2)(R1^2-R5^2).
//
G0 := -&*[R1^2-RR[j]^2 : j in [3..5]];

Afac := t^2*m^2-t^2*m+t*m^2-2*t*m+m-2;
Bfac := t^2*m+t*m^2+2*t*m-t+m;
Cfac := t^2*m^2+t^2*m+t*m^2+2*t*m+t+2;
G0factored := m*(m-1)*t*(t+1)*(t-m+2)*(t*m+1)
    *(t*m+m-1)*(t*m+t+1)*Afac*Bfac*Cfac;
assert G0 eq -G0factored;

print "PAIRWISE_LINEAR_FACTORS";
for j in [3..5] do
    print j, Factorization(R1-RR[j]), Factorization(R1+RR[j]);
end for;
print "G0_FACTORIZATION", Factorization(G0);
print "DEGREES_T_M", Degree(G0,t), Degree(G0,m), TotalDegree(G0);

print "DISCRIMINANT_IN_M", Factorization(Discriminant(G0,m));
print "DISCRIMINANT_IN_T", Factorization(Discriminant(G0,t));

// The two coordinate projections of the complete rational chart both
// give squarefree degree-12 hyperelliptic fibers, hence genus 5.  This
// certifies that the apparent elliptic structure below is an S_3
// quotient, not a hidden birational rewriting of either chart fiber.
Ktau<tau> := FunctionField(Q);
PM<mu> := PolynomialRing(Ktau);
Cm := HyperellipticCurve(PM!Evaluate(G0,[tau,mu]));
assert Genus(Cm) eq 5;
Krho<rho> := FunctionField(Q);
PT<theta> := PolynomialRing(Krho);
Ct := HyperellipticCurve(PT!Evaluate(G0,[theta,rho]));
assert Genus(Ct) eq 5;
print "GENERIC_CHART_FIBER_GENERA", Genus(Cm), Genus(Ct);

// Pass to the quotient by permutations of R3,R4,R5.  Normalize R1=1,
// set q=R2/R1 and s=(R3*R4+R3*R5+R4*R5)/R1^2.  The CK equations give
// e1=-(1+q), e3=(1+q)(q-s), and the first radicand becomes a quartic
// in q.  All denominators below are even powers, so these are exact
// squareclass identities on R1 ne 0.
K := FieldOfFractions(P);
q := K!(R2/R1);
s := K!((R3*R4+R3*R5+R4*R5)/R1^2);
gquot := q*(s-q)*(2-q^2+(q+2)*s);
assert K!(G0/R1^6) eq gquot;

// Cubic splitting polynomial of the normalized complementary roots.
U<u> := PolynomialRing(K);
splitpoly := u^3+(1+q)*u^2+s*u-(1+q)*(q-s);
for r in [R3/R1,R4/R1,R5/R1] do
    assert Evaluate(splitpoly,K!r) eq 0;
end for;

// Quartic-to-cubic map.  If v^2=gquot and q*s*(s+1) ne 0, put
// X=2*s*(s+1)/q and Y=2*s*(s+1)*v/q^2.  The quotient is the elliptic
// curve displayed below.
F<ss> := FunctionField(Q);
EF<X> := PolynomialRing(F);
aquart := 2*ss*(ss+1);
bquart := ss^2-2*ss-2;
cquart := -2*ss;
dquart := F!1;
weier := X^3+bquart*X^2+aquart*cquart*X+aquart^2*dquart;
assert weier eq
    (X-2*(ss+1))*(X^2+ss^2*X-2*ss^2*(ss+1));
E := EllipticCurve([F|0,bquart,0,aquart*cquart,aquart^2*dquart]);
discE := Discriminant(E);
c4E := cInvariants(E)[1];
assert discE eq 256*ss^2*(ss+1)^4*(ss^2+8*ss+8);
assert c4E eq 16*(ss^4+8*ss^3+12*ss^2+8*ss+4);
print "S3_QUOTIENT_QUARTIC",
      "v^2=q(s-q)(2-q^2+(q+2)s)";
print "S3_QUOTIENT_WEIERSTRASS", E;
print "WEIERSTRASS_FACTORIZATION", Factorization(weier);
print "WEIERSTRASS_DISCRIMINANT",
      Factorization(Numerator(Discriminant(E)));
print "WEIERSTRASS_J_INVARIANT", jInvariant(E);

// This is a rational elliptic surface.  Its geometric singular fibers
// are I_2 at s=0, I_4 at s=-1 and s=infinity, and two I_1 fibers at the
// roots of s^2+8s+8.  Shioda--Tate therefore gives geometric MW rank
// 10-2-(1+3+3)=1.  The visible order-four section is rational.
P4 := E![F|0,2*ss*(ss+1),1];
T2 := E![F|2*(ss+1),0,1];
assert 2*P4 eq T2 and 2*T2 eq E!0;

// A non-torsion geometric section is already visible on the quartic:
// q=-1, v=i(s+1).  Complex conjugation negates it.  Since the geometric
// rank is one, the Q(s)-rank is zero.  We certify non-torsion after the
// good specialization s=2 (Order returns 0 for an infinite-order point).
ZI<ii0> := PolynomialRing(Q);
KI<ii> := NumberField(ii0^2+1);
E2I := EllipticCurve([KI|0,-2,0,-48,144]);
QI := E2I![KI|-12,36*ii,1];
assert Order(QI) eq 0;

// Generic torsion injects into the good specialization s=1, whose
// rational torsion is Z/4.  Together with P4 this proves
// E(Q(s)) = Z/4; all four sections are boundary sections q=0,s,infinity.
Eone := EllipticCurve([Q|0,-3,0,-8,16]);
Tone, tone_map := TorsionSubgroup(Eone);
assert Invariants(Tone) eq [4];
print "ELLIPTIC_SURFACE_FIBERS", "I2@0 I4@-1 I4@infinity I1+I1";
print "ELLIPTIC_SURFACE_GEOMETRIC_MW_RANK", 1;
print "ELLIPTIC_SURFACE_RATIONAL_MW", "Z/4 (rank 0)";
print "ELLIPTIC_SURFACE_NONBOUNDARY_RATIONAL_SECTIONS", 0;

// Fixed-q labelled model.  Write r1=1,r2=qq,r3=xx.  The remaining
// two roots have the displayed sum and product.  Their discriminant
// must be a square, as must G0.  Clearing denominators by multiplying
// numerator and denominator preserves squareclasses.
Qqx<qq,xx> := PolynomialRing(Q,2);
Sfix := -(1+qq+xx);
Pfix := (qq+xx)*(1+qq)*(1+xx)/(1+qq+xx);
Dfix := Numerator(Sfix^2-4*Pfix)*Denominator(Sfix^2-4*Pfix);
Gfix0 := -(1-xx^2)*((1+Pfix)^2-Sfix^2);
Gfix := Numerator(Gfix0)*Denominator(Gfix0);

// Delete the common even factor (qq+xx+1)^2 from the G0 branch
// polynomial.  The resultant locates every collision between the two
// remaining branch sets.
Gfixsf := &*[z[1] : z in Factorization(Gfix) | IsOdd(z[2])];
branchres := Resultant(Dfix,Gfixsf,xx);
assert Factorization(branchres) eq [
    <qq,11>, <qq-1,2>, <qq+1,4>, <qq+2,1>,
    <qq^2-2*qq-4,1>, <qq^2+3*qq+Q!8/3,2>
];
print "FIXED_Q_BRANCH_RESULTANT", Factorization(branchres);

// The only non-boundary rational constant-q fiber singled out by this
// resultant is q=-2.  Write x for the first complementary root.
// Eliminating the other two roots gives
//
//  s=(x^2-x^3+2)/(x-1),
//  G0 ~ -x(x-1)(x+1)(x-2),
//
// where ~ means equality up to a nonzero rational square.  Hence the
// first-radicand condition already lies on the genus-one quartic Cminus2.
PX<x> := PolynomialRing(Q);
sminus2 := (x^2-x^3+2)/(x-1);
fminus2 := -x*(x-1)*(x+1)*(x-2);
splitminus2 := (x-1)*(x^3+x^2-x-9);
assert 4*(sminus2+2) eq 4*fminus2/(x-1)^2;

// The second quartic is the square condition for the discriminant of
// the remaining quadratic (so that the other two complementary roots
// are rational).  The fiber product has genus 1+1+2=4, but it is not
// needed: Cminus2 itself has only four rational points, all ramification
// points with G0=0.
Cminus2 := HyperellipticCurve(fminus2);
assert Genus(Cminus2) eq 1;
pminus2 := Cminus2![0,0,1];
Eminus2, phiminus2 := EllipticCurve(Cminus2,pminus2);
Eminus2min, minmapminus2 := MinimalModel(Eminus2);
assert aInvariants(Eminus2min) eq [ Q!0,Q!1,Q!0,Q!-4,Q!-4 ];
rlo, rhi := RankBounds(Eminus2min);
assert rlo eq 0 and rhi eq 0;
Tminus2, tmapminus2 := TorsionSubgroup(Eminus2min);
assert Invariants(Tminus2) eq [ 2,2 ];

// Enumerate the entire Mordell--Weil group (rank zero, four torsion
// points) and pull it back to the quartic.  The resulting x-coordinates
// are exactly {-1,0,1,2}; every one has fminus2=0 and is CK boundary.
xminus2 := {};
for z in Tminus2 do
    ep := tmapminus2(z);
    cp := Inverse(phiminus2)(Inverse(minmapminus2)(ep));
    Include(~xminus2,Q!(cp[1]/cp[3]));
end for;
assert xminus2 eq { Q!-1,Q!0,Q!1,Q!2 };
assert &and[Evaluate(fminus2,a) eq 0 : a in xminus2];

// The simultaneous square conditions before imposing the other three
// orbit-12 radicands form a V_4 cover of the x-line.  Its three quotient
// genera are 1,1,2, hence its smooth projective genus is 4.
fthirdminus2 := -x*(x+1)*(x-2)*(x^3+x^2-x-9);
assert Genus(HyperellipticCurve(splitminus2)) eq 1;
assert Genus(HyperellipticCurve(fthirdminus2)) eq 2;

print "Q_MINUS_2_FIRST_RADICAND_QUARTIC", fminus2;
print "Q_MINUS_2_SPLIT_COMPLETION_QUARTIC", splitminus2;
print "Q_MINUS_2_FIBER_PRODUCT_GENUS", 4;
print "Q_MINUS_2_ELLIPTIC_MODEL", Eminus2min;
print "Q_MINUS_2_RANK_BOUNDS", rlo, rhi;
print "Q_MINUS_2_TORSION", Invariants(Tminus2);
print "Q_MINUS_2_ALL_RATIONAL_X", xminus2;
print "Q_MINUS_2_NO_SMOOTH_POINT", true;

quit;
