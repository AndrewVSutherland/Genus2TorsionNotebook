// claude_ari_surface_verify.m
// Verify Ari's (A,B,C)-surface for M(2,2,2,12):
//  PART A: symbolic identities (email formulas vs document S1..S4; sextic factorizations)
//  PART B: recover both known (2,2,2,12) curves from (A:B:C) points; Igusa match
//  PART C: the two special-fiber genus-1 curves E2 (t=2 line pair) and Eh (t=1/2, z-split)
//          -> ranks; positive rank of either = infinitely many surface points
SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals();

print "==== PART A: symbolic identities ====";
P3<A,B,C> := PolynomialRing(Q, 3);
FF := FieldOfFractions(P3);

// email formulas
sE := (A^3 - A*B^2 + 2*A*C^2)/(2*A^2*B - 2*B^3 + 2*B*C^2);
mE := (-2*A^2*B^2 + 2*B^4 + 2*A^2*C^2 - 6*B^2*C^2 + 4*C^4)/(2*A^3*B - 2*A*B^3 + 2*A*B*C^2);
nE := FF!A/FF!B;

// document parametrization chain: (u,v,w) -> (n,m,s), then u=A/C, v=(A^2-B^2+2C^2)/(AC), w=B/C
uu := FF!A/FF!C;
vv := (A^2 - B^2 + 2*C^2)/(A*C);
ww := FF!B/FF!C;
nD := uu/ww;
mD := vv*(ww^2 + uu^2 - uu*vv)/(ww*(ww^2 - uu^2 - uu*vv));
sD := -vv*uu^2/(ww*(ww^2 + uu^2 - uu*vv));
printf "chain n eq email n: %o\n", nD eq nE;
printf "chain m eq email m: %o\n", mD eq mE;
printf "chain s eq email s: %o\n", sD eq sE;

// email sextics = factored forms
EmailF := A^4*B^2 - 2*A^2*B^4 + B^6 - A^4*C^2 + 3*A^2*B^2*C^2 - 2*B^4*C^2 - A^2*C^4 + B^2*C^4;
EmailG := A^6 - 2*A^4*B^2 + A^2*B^4 + 2*A^4*C^2 - 3*A^2*B^2*C^2 + B^4*C^2 + 2*A^2*C^4 - 2*B^2*C^4 + C^6;
Fsx := (B^2-A^2)*(B^2-C^2)*(B^2-A^2-C^2);
Gsx := (B^2-A^2-C^2)*(B^2*(A^2+C^2) - (A^4+A^2*C^2+C^4));
printf "EmailF eq (B2-A2)(B2-C2)(B2-A2-C2): %o\n", EmailF eq Fsx;
printf "EmailG eq (B2-A2-C2)(B2(A2+C2)-(A4+A2C2+C4)): %o\n", EmailG eq Gsx;

// document S-quantities
S1f := func< s,m,n | m*(n-s)*(n*m + n*s - m*s - 2*s^2) >;
S2f := func< s,m,n | 2*(m + 2*s)*(n - 2*s)*(n*m + n*s - m*s - 2*s^2) >;
S3f := func< s,m,n | -s*n*(n*m + n*s - m*s - 2*s^2) >;
S4f := func< s,m,n | -s*m*(m + 2*s)*n*(n - 2*s)*(n - s)*(n*m - 2*n*s - 2*m*s) >;

function SquareClass(r)  // r in FF nonzero: return unit u in Q with r = u * (square), and whether square
  nu := Numerator(r); de := Denominator(r);
  pol := nu*de;   // r = pol / de^2
  fac := Factorization(pol);
  pr := #fac eq 0 select Parent(pol)!1 else &*[ t[1]^t[2] : t in fac ];
  unit := pol div pr;   // Factorization drops the unit: reconstruct (lab convention S1)
  assert unit*pr eq pol;
  oddpart := #fac eq 0 select Parent(pol)!1 else &*[ t[1]^(t[2] mod 2) : t in fac ];
  isq := (oddpart eq 1) and IsSquare(Q!unit);
  return isq, Q!unit, oddpart;
end function;

s1v := S1f(sE,mE,nE); s2v := S2f(sE,mE,nE); s3v := S3f(sE,mE,nE); s4v := S4f(sE,mE,nE);
isq, un, od := SquareClass(s1v/FF!Fsx);
printf "S1/Fsextic is square: %o (unit %o, oddpart %o)\n", isq, un, od;
isq, un, od := SquareClass(s4v/FF!Gsx);
printf "S4/Gsextic is square: %o (unit %o, oddpart %o)\n", isq, un, od;
isq, un, od := SquareClass(s2v);
printf "S2 is square identically: %o (unit %o, oddpart %o)\n", isq, un, od;
isq, un, od := SquareClass(s3v);
printf "S3 is square identically: %o (unit %o, oddpart %o)\n", isq, un, od;

print "==== PART B: the two known curves ====";
DiscPolys := [P3| C, B, B-C, B+C, B-A, A, A-C, A+C, A+B, B^2-A^2-2*C^2, B^2-A^2-C^2,
  B^2-A^2-A*C-C^2, B^2-A^2+A*C-C^2, B^4-A^2*B^2+2*A^2*C^2-2*B^2*C^2+C^4,
  A^4-A^2*B^2+A^2*C^2-B^2*C^2+C^4, A^4-2*A^2*B^2+B^4+2*A^2*C^2-B^2*C^2];

function EmailSMN(Av, Bv, Cv)
  s := (Av^3 - Av*Bv^2 + 2*Av*Cv^2)/(2*Av^2*Bv - 2*Bv^3 + 2*Bv*Cv^2);
  m := (-2*Av^2*Bv^2 + 2*Bv^4 + 2*Av^2*Cv^2 - 6*Bv^2*Cv^2 + 4*Cv^4)/(2*Av^3*Bv - 2*Av*Bv^3 + 2*Av*Bv*Cv^2);
  n := Av/Bv;
  return s, m, n;
end function;

Px<x> := PolynomialRing(Q);
function QuinticModel(s, m, n)
  a := (m + 2*s)*(n - s)*(n*m - 2*n*s + 4*s^2)*(n*m + n*s - m*s - 2*s^2)*(n*m + 4*n*s - 4*s^2);
  b := -8*(m + 2*s)*(n - 2*s)*s^2*(n - s)^2*(n*m + n*s - m*s - 2*s^2);
  c := -s*(n - 2*s)*(n*m + n*s - m*s - 2*s^2)*(n*m + 4*n*s - 4*s^2)*(n*m + 4*n*s - 2*m*s - 4*s^2);
  d := -s*(n - s)*(m + 2*s)^2*(n - 2*s)^2*(n*m + 4*n*s - 4*s^2);
  return x*(x+a)*(x+b)*(x+c)*(x+d), [a,b,c,d];
end function;

// recorded data
smn2 := [Q| 2208, -8303, -7200];        // curve #2
smn1 := [Q| 336396, -689185, -166464];  // curve #1
f2rec := 36*x^6+36750*x^5-462983772*x^4-301623595823*x^3+1518598238654317*x^2+397058962729817115*x-1282993930035013443975;
h2rec := x^2+x;
f1rec := 756900*x^6 + 737595570*x^5 + 150572203590*x^4 - 15854483576121*x^3 - 530648977741620*x^2 + 32014154874551031*x + 830742747091037849;
h1rec := x^2+1;
C2rec := HyperellipticCurve(f2rec, h2rec);
C1rec := HyperellipticCurve(f1rec, h1rec);
G2rec2 := G2Invariants(C2rec);
G2rec1 := G2Invariants(C1rec);

ABCcands := [ [120,218,143], [408,437,143] ];
recsmn := [ smn2, smn1 ];
recG2  := [ G2rec2, G2rec1 ];
for i in [1..2] do
  T := ABCcands[i];
  Av := Q!T[1]; Bv := Q!T[2]; Cv := Q!T[3];
  printf "--- candidate (A:B:C) = (%o:%o:%o) for curve #%o ---\n", T[1], T[2], T[3], 3-i;
  sv, mv, nv := EmailSMN(Av, Bv, Cv);
  printf "email (s,m,n) = (%o, %o, %o)\n", sv, mv, nv;
  rs := recsmn[i][1]; rm := recsmn[i][2]; rn := recsmn[i][3];
  projmatch := (sv*rm - mv*rs eq 0) and (sv*rn - nv*rs eq 0) and (mv*rn - nv*rm eq 0);
  printf "projective match with recorded (s:m:n): %o\n", projmatch;
  Fv := Evaluate(Fsx, [Av,Bv,Cv]); Gv := Evaluate(Gsx, [Av,Bv,Cv]);
  okF, yv := IsSquare(Fv); okG, zv := IsSquare(Gv);
  printf "F value square: %o (y=%o); G value square: %o (z=%o)\n", okF, yv, okG, zv;
  disc_ok := forall{ dp : dp in DiscPolys | Evaluate(dp, [Av,Bv,Cv]) ne 0 };
  printf "all 16 discriminant factors nonzero: %o\n", disc_ok;
  tval := Evaluate(B^2-A^2-C^2, [Av,Bv,Cv]) / (Av*Cv);
  printf "conic-pencil parameter t = (B^2-A^2-C^2)/(AC) = %o\n", tval;
  fq, abcd := QuinticModel(rs, rm, rn);
  Cq := HyperellipticCurve(fq);
  match := G2Invariants(Cq) eq recG2[i];
  printf "quintic model from recorded (s,m,n) matches recorded curve (G2Invariants): %o\n", match;
  fq2, _ := QuinticModel(sv, mv, nv);
  match2 := G2Invariants(HyperellipticCurve(fq2)) eq recG2[i];
  printf "quintic model from email-affine (s,m,n) matches too: %o\n", match2;
end for;

print "==== PART C: special-fiber genus-1 curves ====";
// E2: t=2 member (line pair B = A+C).  y1^2 = 4a^2+10ac+4c^2, z1^2 = 4a^2+2ac+4c^2
P3s<aa,cc,yy,zz> := ProjectiveSpace(Q, 3);
CC2 := Curve(P3s, [yy^2 - (4*aa^2+10*aa*cc+4*cc^2), zz^2 - (4*aa^2+2*aa*cc+4*cc^2)]);
printf "E2 quadric intersection nonsingular: %o\n", IsNonsingular(CC2);
pt2 := CC2![1,0,2,2];
E2, phi2 := EllipticCurve(CC2, pt2);
E2m := MinimalModel(E2);
printf "E2 = %o\n", aInvariants(E2m);
printf "E2 conductor = %o\n", Conductor(E2m);
r1, r2 := RankBounds(E2m);
printf "E2 RankBounds: [%o, %o]\n", r1, r2;
T2 := TorsionSubgroup(E2m);
printf "E2 torsion: %o\n", Invariants(T2);

// Eh: t=1/2 member: conic 2B^2 = 2A^2+AC+2C^2, quartic y4^2 = 2(A+2C)(2A+C)
P2h<ha,hc,hb> := ProjectiveSpace(Q, 2);
ConH := Conic(P2h, 2*hb^2 - 2*ha^2 - ha*hc - 2*hc^2);
pth := ConH![1,0,1];
parh := Parametrization(ConH, pth);
dp := DefiningPolynomials(parh);   // [A(u,v), C(u,v), B(u,v)]
Puv := Parent(dp[1]);
printf "Eh conic parametrization: A=%o, C=%o, B=%o\n", dp[1], dp[2], dp[3];
q4uv := 2*(dp[1] + 2*dp[2])*(2*dp[1] + dp[2]);
uP := Puv.1; vP := Puv.2;
q4 := Evaluate(q4uv, [x, 1]);
printf "Eh quartic q(u) = %o  (degree %o)\n", q4, Degree(q4);
// leading/infinity behavior: also check value at (1:0)
qinf := Evaluate(q4uv, [1, 0]);
printf "q at infinity (u:v)=(1:0): %o (square: %o)\n", qinf, IsSquare(qinf);
// Jacobian via classical binary quartic invariants: q = a u^4 + b u^3 + c u^2 + d u + e
cf := Coefficients(q4);  // ascending
while #cf lt 5 do Append(~cf, 0); end while;
e0 := cf[1]; d0 := cf[2]; c0 := cf[3]; b0 := cf[4]; a0 := cf[5];
Iinv := 12*a0*e0 - 3*b0*d0 + c0^2;
Jinv := 72*a0*c0*e0 - 27*a0*d0^2 - 27*b0^2*e0 + 9*b0*c0*d0 - 2*c0^3;
Eh := MinimalModel(EllipticCurve([0,0,0,-27*Iinv,-27*Jinv]));
printf "Eh Jacobian = %o\n", aInvariants(Eh);
printf "Eh conductor = %o\n", Conductor(Eh);
rh1, rh2 := RankBounds(Eh);
printf "Eh RankBounds: [%o, %o]\n", rh1, rh2;
Th := TorsionSubgroup(Eh);
printf "Eh torsion: %o\n", Invariants(Th);

// small direct point searches on both fibers to harvest explicit points
print "-- direct point harvest on E2 (line B=A+C): a/c height <= 300 --";
cnt2 := 0;
for hgt in [1..300] do
  for aN in [-hgt..hgt] do
    for s in [ [aN, hgt], [hgt, aN] ] do
      av := s[1]; cv := s[2];
      if Abs(av) ne hgt and Abs(cv) ne hgt then continue; end if;
      if GCD(av, cv) ne 1 then continue; end if;
      F1v := 2*(2*av+cv)*(av+2*cv);
      if F1v lt 0 or not IsSquare(F1v) then continue; end if;
      F2v := 2*(2*av^2+av*cv+2*cv^2);
      if F2v lt 0 or not IsSquare(F2v) then continue; end if;
      Avv := av; Bvv := av+cv; Cvv := cv;
      printf "E2 point (a:c)=(%o:%o) -> (A:B:C)=(%o:%o:%o) disc_ok=%o\n", av, cv, Avv, Bvv, Cvv,
        forall{ dpp : dpp in DiscPolys | Evaluate(dpp, [Q!Avv,Q!Bvv,Q!Cvv]) ne 0 };
      cnt2 +:= 1;
    end for;
  end for;
end for;
printf "E2 direct points found (height<=300): %o\n", cnt2;

print "-- direct point harvest on Eh (conic 2B^2=2A^2+AC+2C^2): u height <= 300 --";
cnth := 0;
for hgt in [1..300] do
  for aN in [-hgt..hgt] do
    for s in [ [aN, hgt], [hgt, aN] ] do
      uv := s[1]; vv2 := s[2];
      if Abs(uv) ne hgt and Abs(vv2) ne hgt then continue; end if;
      if GCD(uv, vv2) ne 1 then continue; end if;
      Avv := Q!Evaluate(dp[1], [uv, vv2]); Cvv := Q!Evaluate(dp[2], [uv, vv2]); Bvv := Q!Evaluate(dp[3], [uv, vv2]);
      if Avv eq 0 and Bvv eq 0 and Cvv eq 0 then continue; end if;
      den := LCM([Denominator(Avv), Denominator(Bvv), Denominator(Cvv)]);
      Ai := Integers()!(Avv*den); Bi := Integers()!(Bvv*den); Ci := Integers()!(Cvv*den);
      qv := 2*(Ai + 2*Ci)*(2*Ai + Ci);
      if qv lt 0 or not IsSquare(qv) then continue; end if;
      g := GCD([Ai, Bi, Ci]);
      Avv := Ai div g; Bvv := Bi div g; Cvv := Ci div g;
      printf "Eh point (u:v)=(%o:%o) -> (A:B:C)=(%o:%o:%o) disc_ok=%o\n", uv, vv2, Avv, Bvv, Cvv,
        forall{ dpp : dpp in DiscPolys | Evaluate(dpp, [Q!Avv,Q!Bvv,Q!Cvv]) ne 0 };
      cnth +:= 1;
    end for;
  end for;
end for;
printf "Eh direct points found (height<=300): %o\n", cnth;

print "ALL_DONE";
quit;
