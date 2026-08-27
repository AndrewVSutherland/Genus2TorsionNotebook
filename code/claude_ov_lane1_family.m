// claude_ov_lane1_family.m -- Lane 1, overnight 2026-07-25.
//
// THE u = -1/2 SLICE OF THE CONTACT-7 THREE-ROOT SURFACE IS AN ELLIPTIC CURVE
// OF RANK 1, hence an INFINITE family of [2,2,14] genus-2 Jacobians.
//
// Contact-7 chart: h = 1 - (7/2)x + a x^2 + b x^3,  f = (h^2+(x-1)^7)/x^2
// (monic quintic), marked class [P-inf] of order 7 at P = (1,h(1)).
// A rational root r of f <=> r = 1-v^2 with v a rational root of
//   Q(v) = v^5 + c4 v^4 + (2c4-1)v^3 + (c4+c0-1/2)v^2 + 2c0 v + c0,
//   c4 = b+2, c0 = 5/2-(a+b).
// Three rational roots (s,t,u) <=> R(s,t,u) = 0, the S3-symmetric tridegree
// (3,3,3) three-root surface.
//
// LANE-1 FINDING.  Put u = c and symmetrise in (s,t) via p = s+t, q = st:
// R becomes a PLANE CUBIC in (p,q) through (p,q)=(0,0) with linear part
// c^2(2c+1) p, so the origin is a SINGULAR point of the fibre cubic exactly
// for c in {0,-1/2}; c = 0 is degenerate (root r = 1 = the marked point), so
// c = -1/2 is the unique nondegenerate nodal fibre.  Its cubic
//     -p^3 + 3 p q^2 + 2 p q + 2 q^3 = 0
// is nodal at the origin, hence rational: p = m q gives
//     q = 2m/((m+1)^2(m-2)),   p = 2m^2/((m+1)^2(m-2)),
// and s,t are rational iff p^2-4q is a square iff
//     w^2 = -m^4 + 6m^2 + 4m = -m(m+2)(m^2-2m-2),
// with then
//     s = (m^2+w)/((m+1)^2(m-2)),   t = (m^2-w)/((m+1)^2(m-2)),  u = -1/2.
// The quartic is the elliptic curve  E : y^2 = x^3 + 6x^2 - 16  (conductor
// 288) via x = 4/m, y = 4w/m^2;  E(Q) = Z/2 + Z  (rank 1).
//
// Usage:  magma -b KMAX:=10 TORSMAX:=6 code/claude_ov_lane1_family.m

SetColumns(0);
if not assigned MemGB  then MemGB  := 8;  elif Type(MemGB)  eq MonStgElt then MemGB  := StringToInteger(MemGB);  end if;
SetMemoryLimit(MemGB*10^9);
if not assigned KMAX   then KMAX   := 10; elif Type(KMAX)   eq MonStgElt then KMAX   := StringToInteger(KMAX);   end if;
if not assigned TORSMAX then TORSMAX := 6; elif Type(TORSMAX) eq MonStgElt then TORSMAX := StringToInteger(TORSMAX); end if;

Q := Rationals(); Z := Integers();
P<x> := PolynomialRing(Q);

/* ---------------- 0. rebuild R(s,t,u) symbolically ---------------- */
S3<sv,tv,uv> := PolynomialRing(Q,3);
FS := FieldOfFractions(S3);
Afun := func<w0 | (2*w0^5+4*w0^4+6*w0^3+8*w0^2+10*w0+5)/(2*(w0+1)^2)>;
bST := (Afun(FS!tv)-Afun(FS!sv))/(FS!(sv^2-tv^2));
aST := Afun(FS!sv) - bST*(1-FS!(sv^2));
Enum := Afun(FS!uv) - aST - bST*(1-FS!(uv^2));
NN := S3!Numerator(Enum);
ok, Rfull := IsDivisibleBy(NN, S3!((uv-sv)*(uv-tv)));
assert ok;
// strip integer content
dl := LCM([Denominator(c) : c in Coefficients(Rfull)]);
Rfull := S3!(dl*Rfull);
gc := GCD([Z!c : c in Coefficients(Rfull)]);
Rfull := S3!(Rfull div (S3!gc));
printf "R rebuilt: tridegree = (%o,%o,%o)\n", Degree(Rfull,1), Degree(Rfull,2), Degree(Rfull,3);
printf "R S3-symmetric: cyc=%o transp=%o\n",
   Evaluate(Rfull,[tv,uv,sv]) eq Rfull, Evaluate(Rfull,[tv,sv,uv]) eq Rfull;
e1 := sv+tv+uv; e2 := sv*tv+sv*uv+tv*uv; e3 := sv*tv*uv;
Rrec := 2*(e2-2)*e3^2 + (4*e2^2+4*e1*e2-2*e1-4*e2-1)*e3 + e2*(2*(e1+e2)^2+e1);
printf "R equals recorded e-form: %o\n", Rrec eq Rfull;
printf "R = %o\n", Rfull;

elevens := [
 [Q| -10,-10/7,-1/2], [Q| -5,-15/8,-15/22], [Q| -3,-3/4,-3/5], [Q| -15/8,-15/19,-1/2],
 [Q| -5/18,-10/49,4/17], [Q| -4/9,-4/25,4/17], [Q| -511/61,-511/625,-1/2],
 [Q| -165/41,-33/16,-165/289], [Q| -164/297,-1/2,164/361], [Q| -17/50,-34/189,34/121],
 [Q| -1/2,-13/49,13/50] ];
printf "eleven recorded points all on R: %o\n",
   &and[ Evaluate(Rfull, pt) eq 0 : pt in elevens ];

/* ---------------- 1. the parametrisation satisfies R identically --------- */
Fm<m> := FunctionField(Q);
Pw<Wv> := PolynomialRing(Fm);
Kw<wq> := ext<Fm | Wv^2 - (-m^4+6*m^2+4*m)>;
den := (m+1)^2*(m-2);
sPar := (Kw!(m^2)+wq)/Kw!den; tPar := (Kw!(m^2)-wq)/Kw!den; uPar := Kw!(-1/2);
printf "PARAMETRISATION identically on R: %o\n", Evaluate(Rfull,[sPar,tPar,uPar]) eq 0;
printf "  s+t = %o\n  s*t = %o\n", sPar+tPar, sPar*tPar;

/* ---------------- 2. walk E(Q) ---------------- */
Gfun := func<u0 | -(u0^5 - u0^3 - u0^2/2)/(u0+1)^2>;
E := EllipticCurve([0,6,0,0,-16]);
printf "E = %o  conductor %o\n", E, Conductor(E);
Gp := E![4,-12]; Tp := E![-2,0];
printf "generator %o (infinite order: %o), 2-torsion %o\n", Gp, Order(Gp) eq 0, Tp;

results := [* *];
korder := [0] cat &cat[[j,-j] : j in [1..KMAX]];
for k in korder do
 for eps in [0,1] do
  Pt := k*Gp + eps*Tp;
  if Pt eq E!0 then continue; end if;
  xc := Pt[1]/Pt[3]; yc := Pt[2]/Pt[3];
  if xc eq 0 then continue; end if;
  mm := 4/xc; ww := yc*mm^2/4;
  assert ww^2 eq -mm^4+6*mm^2+4*mm;
  d := (mm+1)^2*(mm-2);
  if d eq 0 then printf "k=%o eps=%o : m=%o degenerate denominator\n",k,eps,mm; continue; end if;
  s0 := (mm^2+ww)/d; t0 := (mm^2-ww)/d; u0 := Q!(-1/2);
  assert Evaluate(Rfull,[s0,t0,u0]) eq 0;
  vs := [s0,t0,u0];
  if #{v : v in vs} ne 3 then printf "k=%o eps=%o m=%o : repeated v\n",k,eps,mm; continue; end if;
  if #{v^2 : v in vs} ne 3 then printf "k=%o eps=%o m=%o : repeated v^2\n",k,eps,mm; continue; end if;
  if (-1 in vs) or (0 in vs) then printf "k=%o eps=%o m=%o : v in {-1,0}\n",k,eps,mm; continue; end if;
  c4 := (Gfun(s0)-Gfun(t0))/(s0^2-t0^2); c0 := Gfun(s0) - c4*s0^2;
  b := c4 - 2; a := 9/2 - c0 - c4;
  h := 1 - (7/2)*x + a*x^2 + b*x^3;
  num := h^2 + (x-1)^7;
  if num mod x^2 ne 0 then printf "k=%o eps=%o : x^2 nondivisor\n",k,eps; continue; end if;
  f := num div x^2;
  if Degree(f) ne 5 or Discriminant(f) eq 0 then printf "k=%o eps=%o m=%o : f degenerate (deg %o)\n",k,eps,mm,Degree(f); continue; end if;
  ok3 := &and[ Evaluate(f, 1-v^2) eq 0 : v in vs ];
  ft := Sort([Degree(tt[1]) : tt in Factorization(f)]);
  mden := LCM([Denominator(c) : c in Coefficients(f)]);
  F := P![ Z!(Coefficient(f,i)*mden^(6-i)) : i in [0..5] ];
  C := HyperellipticCurve(F); J := Jacobian(C);
  tr := #Invariants(TwoTorsionSubgroup(J));
  Ppt := C ! [mden*1, mden^3*Evaluate(h,1)];
  D := J ! (Ppt - PointsAtInfinity(C)[1]);
  ordD := Order(D);
  hgt := Maximum([Abs(Numerator(cc)) + Abs(Denominator(cc)) : cc in Coefficients(f)]);
  printf "k=%o eps=%o m=%o\n  s=%o\n  t=%o\n  u=-1/2  (a,b)=(%o,%o)\n  roots-ok=%o ftype=%o 2rank=%o markedord=%o log10ht=%o\n",
    k, eps, mm, s0, t0, a, b, ok3, ft, tr, ordD, Ilog(10, Maximum(2,hgt));
  tt := [Z| ];
  if Abs(k) le TORSMAX then
    tt := [Z| ii : ii in Invariants(TorsionSubgroup(J))];
    printf "  TORSION=%o\n", tt;
    printf "  F = %o\n", F;
  end if;
  Append(~results, <k,eps,mm,s0,t0,ft,tr,ordD,tt>);
 end for;
end for;

printf "\n=== SUMMARY ===\n";
for r in results do
  printf "k=%o eps=%o m=%o ftype=%o 2rank=%o markedord=%o tors=%o\n", r[1],r[2],r[3],r[6],r[7],r[8],r[9];
end for;
printf "LANE1_FAMILY_DONE\n";
quit;
