// Lane 2: the c=2 slice, in the SYMMETRIC coordinates a = w1+w2, b = w1*w2.
// (z1 = 2; the other four z's are w1..w4 with e1=4, e4=1, 3*e3 = 2*e2.)
// Splitting the four w's into the pair {w1,w2} (sum a, product b) and the residual pair
// (sum 4-a, product 1/b), the relation 3e3 = 2e2 becomes the PLANE CUBIC
//        Cab :  -3ab^2 + 10b^2 + 2a^2 b - 8ab + 3a - 2 = 0,
// which is INVARIANT under the involution iota : (a,b) -> (4-a, 1/b) (swap the two pairs).
//   d1 := a^2 - 4b   square  <=>  w1,w2 rational  -> three rational z -> [2,2,14]
//   d2 := (4-a)^2 - 4/b square <=> w3,w4 rational -> (the iota-image of the same condition)
//   BOTH square  <=>  all five z rational  ->  [2,2,2,14], ORDER 112.
//
// FIXED 2026-07-25: Cab has GENUS 0 (the original script called EllipticCurve on it and
// crashed).  Genus 0 means Cab is a RATIONAL curve -- parametrize it by T and enumerate.
// The genus-1 object is the double cover u^2 = d1 (that is claude_ov_c7z_mwscan.m).
SetColumns(0);
if not assigned BT then BT := 400; elif Type(BT) eq MonStgElt then BT := StringToInteger(BT); end if;
SetMemoryLimit(4*10^9);
Q := Rationals();
R2<a,b> := PolynomialRing(Q,2);
A2 := AffineSpace(R2);
F := -3*a*b^2 + 10*b^2 + 2*a^2*b - 8*a*b + 3*a - 2;
C := Curve(A2,F);
PC := ProjectiveClosure(C);
printf "Cab genus = %o (RATIONAL curve -- no elliptic model exists)\n", Genus(PC);

// the recorded three-root hits that live on the slice z1=2, in (a,b)
known := [ [-1/9,-7/3], [-8/7,19/4], [-61/450,625/114], [297/133,361/525], [49/36,50/63] ];
kab := [];
for p in known do
  aa := p[1]+p[2]; bb := p[1]*p[2];
  printf "  hit (w1,w2)=(%o,%o) -> (a,b)=(%o,%o)  on Cab: %o   d1=%o (sq %o)  d2=%o (sq %o)\n",
     p[1],p[2],aa,bb, Evaluate(F,[aa,bb]) eq 0, aa^2-4*bb, IsSquare(aa^2-4*bb),
     (4-aa)^2-4/bb, IsSquare((4-aa)^2-4/bb);
  Append(~kab,[aa,bb]);
end for;

// ---- parametrization of Cab ----
pt0 := PC ! [-22/9, 7/27, 1];
prm := Parametrization(PC, Place(pt0));
eqs := DefiningEquations(prm);
PP<X> := PolynomialRing(Q);
KT := FieldOfFractions(PP); TT := KT!X;
at := Evaluate(eqs[1],[TT,KT!1])/Evaluate(eqs[3],[TT,KT!1]);
bt := Evaluate(eqs[2],[TT,KT!1])/Evaluate(eqs[3],[TT,KT!1]);
printf "a(T) = %o\nb(T) = %o\n", at, bt;
assert Evaluate(F,[at,bt]) eq 0;

// locate the known hits on the parameter line
for p in kab do
  num := PP ! Numerator(at - p[1]);
  rts := [ r[1] : r in Roots(num) ];
  printf "  (a,b)=(%o,%o)  <->  T = %o\n", p[1], p[2], [ t : t in rts | Evaluate(bt,t) eq p[2] ];
end for;

// ---- enumerate Cab(Q) by T of bounded height and test the two square conditions ----
// square classes (verified in claude_ov_c7z_xdecide.m):
//   d1 ~ G1 = T(T-3)(T^2-24T+36),   d2 ~ G2 = (T-3)(7T-12)(37T^2-96T+36)
G1 := X*(X-3)*(X^2-24*X+36);
G2 := (X-3)*(7*X-12)*(37*X^2-96*X+36);
for tv in [Q| 1/2, 5, -7, 11/3, 100 ] do
  assert IsSquare(Evaluate(at,tv)^2 - 4*Evaluate(bt,tv))
      eq (Evaluate(G1,tv) ge 0 and IsSquare(Evaluate(G1,tv)));
end for;
printf "square classes of d1,d2 confirmed against G1,G2\n";
n3 := 0; n5 := 0; hits3 := []; hits5 := [];
for q in [1..BT] do
  q4 := q^4;
  for p in [-BT..BT] do
    if GCD(Abs(p),q) ne 1 then continue; end if;
    tv := p/q;
    v1 := Evaluate(G1,tv)*q4;  v2 := Evaluate(G2,tv)*q4;
    s1 := v1 ge 0 and IsSquare(Integers()!v1);
    s2 := v2 ge 0 and IsSquare(Integers()!v2);
    if s1 and s2 then n5 +:= 1; Append(~hits5, tv);
    elif s1 then n3 +:= 1; if n3 le 12 then Append(~hits3, tv); end if;
    end if;
  end for;
end for;
printf "CAB_ENUM_DONE BT=%o  three-root T-values=%o (first: %o)  split-all T-values=%o : %o\n",
   BT, n3, hits3, n5, hits5;
printf "(all split-all T are DEGENERATE: see claude_ov_c7z_chabauty.m)\n";
quit;
