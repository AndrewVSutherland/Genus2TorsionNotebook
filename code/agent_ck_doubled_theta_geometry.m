//////////////////////////////////////////////////////////////////////
// Exact geometry of the three doubled-theta loci on the
// Clebsch--Klein [2,2,2,10] surface.
//
// The odd model is
//
//   Z^2 = g(X) = Product_i (1-r_i^2 X),
//   T = (0,1)-infinity,  Order(T)=5,
//
// with Sum r_i = Sum r_i^3 = 0.  Up to sign and S_5, the three
// nontrivial order-10 cases are
//
//   A: T+S12,  B: 2T+S01,  C: 2T+S12.
//
// This script:
//   * constructs their Mumford quadratics by interpolation;
//   * factors the discriminants on the rational CK chart;
//   * computes the normalization genera of the unique odd cores;
//   * certifies that none of the three cores has a smooth rational
//     CK point, using a genus-2 Chabauty quotient for A, a rank-zero
//     elliptic quotient for B, and a real obstruction for C.
//
// Typical run:
//   magma code/agent_ck_doubled_theta_geometry.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(4*10^9);

Q := Rationals();

//////////////////////////////////////////////////////////////////////
// 1. Actual discriminants on a rational CK chart.
//////////////////////////////////////////////////////////////////////

A2<t,m> := AffineSpace(Q,2);
R := CoordinateRing(A2);
PX<X> := PolynomialRing(R);

rr := [
    1+t*(t+2)*m,
    t*m*(m-t-2),
    -1+m+t*(t+1)*m^2,
    1+t-m-t*m^2,
    -(1+t)*(1+t*m^2)
];
assert &+rr eq 0 and &+[r^3 : r in rr] eq 0;

aa := [r^2 : r in rr];
A1 := &+aa;
g := &*[PX | 1-a*X : a in aa];

function QuadraticDiscriminant(num,den)
    ok,q := IsDivisibleBy(num,den);
    assert ok and Degree(q) eq 2;
    return q, Discriminant(q);
end function;

// A: the quadratic through (0,1), W_1, W_2.
ellA := (1-aa[1]*X)*(1-aa[2]*X);
uA,dA := QuadraticDiscriminant(
    g-ellA^2, X*(1-aa[1]*X)*(1-aa[2]*X));

// B: the quadratic tangent at (0,1) and through W_1.
cB := A1*aa[1]/2-aa[1]^2;
ellB := 1-A1*X/2+cB*X^2;
uB,dB := QuadraticDiscriminant(
    g-ellB^2, X^2*(1-aa[1]*X));

// C: the cubic tangent at (0,1) and through W_1,W_2.
cC := aa[1]+aa[2]-A1/2;
ellC := (1-aa[1]*X)*(1-aa[2]*X)*(1+cC*X);
uC,dC := QuadraticDiscriminant(
    ellC^2-g, X^2*(1-aa[1]*X)*(1-aa[2]*X));

print "CK_DOUBLED_THETA_GEOMETRY";
for rec in [<"A_T_plus_S12",dA>,<"B_2T_plus_S01",dB>,
            <"C_2T_plus_S12",dC>] do
    label,d := Explode(rec);
    fac := Factorization(R!d);
    odd := [f[1] : f in fac | IsOdd(f[2])];
    assert #odd eq 1;
    core := odd[1];
    curve := Curve(A2,core);
    print "CHART_CORE",label,
          "factor_degrees_exponents",[<TotalDegree(f[1]),f[2]> : f in fac],
          "bidegree",<Degree(core,t),Degree(core,m)>,
          "total_degree",TotalDegree(core),
          "irreducible",IsIrreducible(core),
          "genus",Genus(curve),
          "points_height_1000",PointSearch(ProjectiveClosure(curve),1000);
end for;

//////////////////////////////////////////////////////////////////////
// 2. Symmetric formulas.
//////////////////////////////////////////////////////////////////////

// For A and C, put
//   s=r1+r2, p=r1*r2,
//   q=r3*r4+r3*r5+r4*r5, u=r3*r4*r5.
// The CK equations give u=s*(p-q).  Up to a nonzero scalar, the
// Mumford quadratics are the following.
S3<s,p,q> := PolynomialRing(Q,3);
PS<Y> := PolynomialRing(S3);
d := q-p;
UA := d*(2+(p+q-2*s^2)*Y-s^2*d*Y^2);
UC := p*d*(2+2*q*Y+p*d*Y^2);
discA := Discriminant(UA);
discC := Discriminant(UC);
assert discA eq d^2*((p+q-2*s^2)^2+8*s^2*d);
assert discC eq 4*p^2*d^2*(q^2-2*p*d);
assert q^2-2*p*d eq (q-p)^2+p^2;
print "SYMMETRIC_A",UA,"disc",Factorization(discA);
print "SYMMETRIC_C",UC,"disc",Factorization(discC),
      "core_sum_of_squares",(q-p)^2+p^2;

// For B, single out r1=r.  For the other four roots put
// e1=-r, e2=q, e3=-r*q, e4=w.  Again up to scalar:
S2<qb,w> := PolynomialRing(Q,2);
PB<Z> := PolynomialRing(S2);
UB := w*(w*Z^2+2*qb*Z+2);
discB := Discriminant(UB);
assert discB eq 4*w^2*(qb^2-2*w);
print "SYMMETRIC_B",UB,"disc",Factorization(discB);

//////////////////////////////////////////////////////////////////////
// 3. Case A: genus-5 discriminant cover and genus-2 quotient.
//////////////////////////////////////////////////////////////////////

// On the smooth open, s and q-p are nonzero: s=0 or q=p forces
// u=r3*r4*r5=0.  Set s=1.  The core is parameterized by
//
//   p=(v+1)^2, q=1+2v-v^2, u=2v^2.
//
// Splitting r1,r2 rationally gives
//
//   v=-(k^2-k+1)/(k^2+1),
//   r1=1/(k^2+1), r2=k^2/(k^2+1).
//
// Full splitting of the complementary cubic forces its discriminant
// D12(k) to be a square.
Pk<k> := PolynomialRing(Q);
D12 := 8*k^12-80*k^11+256*k^10-640*k^9+1097*k^8-1568*k^7+
       1702*k^6-1568*k^5+1097*k^4-640*k^3+256*k^2-80*k+8;
Kk := FieldOfFractions(Pk);
kK := Kk!k;
PKz<zK> := PolynomialRing(Kk);
vK := -(kK^2-kK+1)/(kK^2+1);
qK := 1+2*vK-vK^2;
uK := 2*vK^2;
cubicK := zK^3+zK^2+qK*zK-uK;
assert Discriminant(cubicK) eq (Kk!D12)/(kK^2+1)^6;
assert IsIrreducible(D12);
C5 := HyperellipticCurve(D12);
assert Genus(C5) eq 5;

// D12 is reciprocal.  With h=k+1/k and y'=y/k^3 its quotient is:
Ph<h> := PolynomialRing(Q);
F2 := 8*h^6-80*h^5+208*h^4-240*h^3+145*h^2-48*h+4;
assert Kk!D12 eq kK^6*Evaluate(F2,kK+1/kK);
C2 := HyperellipticCurve(F2);
J2 := Jacobian(C2);
r2lo,r2hi := RankBounds(J2);
assert r2lo eq 0 and r2hi eq 1;
TG2,mpG2 := TorsionSubgroup(J2);
DA := J2![h^2,2-12*h];
assert Order(DA) eq 0;
ptsA := Chabauty(DA : ptC := C2![0,2,1]);
expectedA := {C2![0,2,1],C2![0,-2,1]};
assert ptsA eq expectedA;
print "CASE_A_GENUS5", "irreducible",true,"genus",Genus(C5);
print "CASE_A_GENUS2_QUOTIENT",F2,"rank_bounds",<r2lo,r2hi>,
      "torsion",Invariants(TG2),"infinite_divisor",DA,
      "chabauty_points",ptsA;
print "CASE_A_CONCLUSION",
      "only h=0, but h=k+1/k=0 requires k^2+1=0; no rational lift";

//////////////////////////////////////////////////////////////////////
// 4. Case B: rank-zero elliptic discriminant quotient.
//////////////////////////////////////////////////////////////////////

// On the open, r1 and w are nonzero.  Normalize r1=1.  The core
// q^2=2w leaves the complementary splitting polynomial
//
//   z^4+z^3+Q*z^2+Q*z+Q^2/2.
//
// Complete rational splitting forces its discriminant to be a square.
PQ<QB> := PolynomialRing(Q);
quarticB := Z^4+Z^3+qb*Z^2+qb*Z+qb^2/2;
discQuarticB := Discriminant(quarticB);
assert discQuarticB eq qb^3*(32*qb^3+56*qb^2-35*qb-16)/4;

FB := QB*(32*QB^3+56*QB^2-35*QB-16);
CB := HyperellipticCurve(FB);
P0B := CB![0,0,1];
EB,phiB := EllipticCurve(CB,P0B);
rBlo,rBhi := RankBounds(EB);
assert rBlo eq 0 and rBhi eq 0;
TB,mpTB := TorsionSubgroup(EB);
assert Invariants(TB) eq [3];
ptsB := {@ P0B @};
for tb in TB do
    et := mpTB(tb);
    if et ne EB!0 then
        Include(~ptsB,et @@ phiB);
    end if;
end for;
expectedB := {@ CB![0,0,1], CB![-4,-316,17], CB![-4,316,17] @};
assert ptsB eq expectedB;
assert { p[1]/p[3] : p in ptsB } eq { Q!0, Q!(-4/17) };

Pz<z> := PolynomialRing(Q);
specialB := 289*z^4+289*z^3-68*z^2-68*z+8;
assert IsIrreducible(specialB);
print "CASE_B_ELLIPTIC",MinimalModel(EB),"rank_bounds",<rBlo,rBhi>,
      "torsion",Invariants(TB),"rational_points",ptsB;
print "CASE_B_SPECIALIZATIONS",
      "Q=0 is boundary; Q=-4/17 gives irreducible quartic",specialB;
print "CASE_B_CONCLUSION","no completely split open specialization";

//////////////////////////////////////////////////////////////////////
// 5. Case C: real obstruction.
//////////////////////////////////////////////////////////////////////

// The omitted discriminant factors p=0 and q=p force respectively
// r1*r2=0 and r3*r4*r5=s*(p-q)=0.  Thus they are CK boundary.
// On the open, the remaining equation is positive definite over R.
print "CASE_C_CORE",(q-p)^2+p^2;
print "CASE_C_CONCLUSION",
      "(q-p)^2+p^2=0 has only p=q=0 over Q/R, hence boundary";

print "FINAL", "no smooth rational CK point on any of A,B,C";
quit;
