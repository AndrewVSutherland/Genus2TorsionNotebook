// claude_ari_orbit_structure.m — the symmetric orbit model of Ari's surface.
// Empirical discovery: an orbit (= one (2,2,2,12) curve) is a triple {X1,X2,X3} and pair {B1,B2} with
//   (I)  B1^2+B2^2 = X1^2+X2^2+X3^2
//   (II) B1^2*B2^2 = e2(X1^2, X2^2, X3^2)
// i.e. {B1^2, B2^2} = roots of T^2 - e1(x)T + e2(x), x = squares of the X's.
// Representations: (A,C) = any 2-subset of {X}, B = either B.
SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals(); Px<x> := PolynomialRing(Q);

orb1X := [143, 408, 1015]; orb1B := [437, 1013];   // curve #1
orb2X := [120, 143, 266];  orb2B := [218, 241];    // curve #2

for pair in [ <orb1X, orb1B, "curve#1">, <orb2X, orb2B, "curve#2"> ] do
  X := pair[1]; Bp := pair[2];
  xs := [ t^2 : t in X ];
  e1 := &+xs; e2 := xs[1]*xs[2] + xs[1]*xs[3] + xs[2]*xs[3];
  lhs1 := Bp[1]^2 + Bp[2]^2; lhs2 := (Bp[1]*Bp[2])^2;
  printf "%o: e1 identity %o, e2 identity %o\n", pair[3], lhs1 eq e1, lhs2 eq e2;
end for;

print "==== solve NEWCURVE_3 orbit from its one representation (5364,19661,4165) ====";
Av := Q!5364; Cv := Q!4165; B1 := Q!19661;
k := B1^2 - Av^2 - Cv^2;
u := (Av^2*Cv^2 + B1^2*k)/k;     // X3^2
v := u - k;                       // B2^2
oku, X3 := IsSquare(u);
okv, B2 := IsSquare(v);
printf "X3^2 = %o square: %o\n", u, oku;
printf "B2^2 = %o square: %o\n", v, okv;
if oku and okv then
  printf "X-triple = {5364, 4165, %o}, B-pair = {19661, %o}\n", X3, B2;
  // verify all 6 representations lie on the surface and give the same curve
  F := func< A0,B0,C0 | (B0^2-A0^2)*(B0^2-C0^2)*(B0^2-A0^2-C0^2) >;
  G := func< A0,B0,C0 | (B0^2-A0^2-C0^2)*(B0^2*(A0^2+C0^2) - (A0^4+A0^2*C0^2+C0^4)) >;
  Xs := [Q!5364, Q!4165, X3]; Bs := [B1, B2];
  for i in [1..3] do for j in [i+1..3] do for bb in Bs do
    okF, _ := IsSquare(F(Xs[i], bb, Xs[j]));
    okG, _ := IsSquare(G(Xs[i], bb, Xs[j]));
    printf "rep (A,B,C)=(%o,%o,%o): F sq %o, G sq %o\n", Xs[i], bb, Xs[j], okF, okG;
  end for; end for; end for;
end if;

print "==== converse test: symmetric-model scan (X1<X2<X3<=260), do solutions land on the surface? ====";
// conditions: disc = e1^2-4e2 square, both roots squares. Then verify F,G squares at a representation.
nsol := 0;
for X3v in [3..260] do
  x3 := X3v^2;
  for X2v in [2..X3v-1] do
    x2 := X2v^2;
    for X1v in [1..X2v-1] do
      x1 := X1v^2;
      if GCD([X1v, X2v, X3v]) ne 1 then continue; end if;
      e1 := x1+x2+x3; e2 := x1*x2+x1*x3+x2*x3;
      dd := e1^2 - 4*e2;
      if dd le 0 then continue; end if;
      okd, sq := IsSquare(dd);
      if not okd then continue; end if;
      r1 := (e1+sq) div 2; r2 := (e1-sq) div 2;
      if (e1+sq) mod 2 ne 0 then continue; end if;
      ok1, b1 := IsSquare(r1);
      if not ok1 then continue; end if;
      ok2, b2 := IsSquare(r2);
      if not ok2 then continue; end if;
      nsol +:= 1;
      // check surface membership of representation (X1, b1, X2) etc.
      Fv := (Q!r1-x1)*(Q!r1-x2)*(Q!r1-x1-x2);
      Gv := (Q!r1-x1-x2)*(Q!r1*(x1+x2) - (Q!x1^2+x1*x2+Q!x2^2));
      okF, _ := IsSquare(Fv); okG, _ := IsSquare(Gv);
      printf "SYMSOL {%o,%o,%o} B={%o,%o}: F sq %o G sq %o\n", X1v, X2v, X3v, b1, b2, okF, okG;
    end for;
  end for;
end for;
printf "symmetric solutions found: %o\n", nsol;
print "ALL_DONE";
quit;
