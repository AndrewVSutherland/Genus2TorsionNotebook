// ===========================================================================
// Lane 2 (claude_ov_c7z) -- CLOSE THE ONE REMAINING GAP IN THE c=2 PROOF.
//
// The proof runs  order-112 config on the slice z1=2  ==>  rational point (a,b)
// of the genus-0 plane cubic  Cab_2 : -3ab^2+10b^2+2a^2b-8ab+3a-2 = 0  with
// d1 = a^2-4b and d2 = (4-a)^2-4/b both squares, and then parametrizes Cab_2 by
// T in P^1 and decides the (Z/2)^2-cover X : u^2=G1(T), v^2=G2(T).
//
// GAP: a rational parametrization P^1 -> Cab_2 hits every rational PLACE, but a
// rational SINGULAR point of the plane model whose branches are conjugate is a
// rational point of Cab_2 that is NOT in the image of P^1(Q).  A genus-0 plane
// cubic has exactly one singular point, so at most one point can be missed.
// This script finds every singular point of Cab_2 over Q and tests it directly.
//
// It also re-checks: (i) the poles of a(T), b(T) really are non-configurations,
// (ii) the T=3 fibre of X (where X -> H is undefined) is the only place the
//      quotient map loses information, and its two places are irrational.
// ===========================================================================
SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals();
R2<aa,bb> := PolynomialRing(Q,2);
A2 := AffineSpace(R2);
Fab := -3*aa*bb^2 + 10*bb^2 + 2*aa^2*bb - 8*aa*bb + 3*aa - 2;
Cab := Curve(A2, Fab);
PC := ProjectiveClosure(Cab);
printf "Cab_2 : %o = 0 ; degree %o ; genus %o\n", Fab, Degree(PC), Genus(PC);

printf "\n---- singular points of the projective model ----\n";
Sg := SingularSubscheme(PC);
printf "dim(Sing) = %o ; degree = %o\n", Dimension(Sg), Degree(Sg);
sp := RationalPoints(Sg);
printf "rational singular points: %o\n", sp;
for P in sp do
  printf "  P = %o\n", P;
  pls := Places(PC ! P);
  printf "    places above P: %o ; residue degrees %o\n", #pls, [Degree(p) : p in pls];
  if P[3] eq 0 then printf "    -> at infinity of the (a,b) plane: a or b infinite, NOT a configuration\n"; continue; end if;
  av := P[1]/P[3]; bv := P[2]/P[3];
  printf "    (a,b) = (%o,%o)\n", av, bv;
  if bv eq 0 then printf "    -> b = 0 means a root w_i = 0, contradicts e4 = 1.  NOT a configuration\n"; continue; end if;
  d1 := av^2 - 4*bv; d2 := (4-av)^2 - 4/bv;
  printf "    d1 = %o (square %o) ; d2 = %o (square %o)\n", d1, IsSquare(d1), d2, IsSquare(d2);
  if IsSquare(d1) and IsSquare(d2) then
    printf "    *** would be an ORDER-112 point -- check degeneracy!\n";
    PW<w> := PolynomialRing(Q);
    quart := (w^2-av*w+bv)*(w^2-(4-av)*w+1/bv);
    printf "    quartic = %o ; roots %o\n", quart, Roots(quart);
  else
    printf "    -> NOT an order-112 point\n";
  end if;
end for;

printf "\n---- the parametrization and which rational points it reaches ----\n";
PP<X> := PolynomialRing(Q);
KT := FieldOfFractions(PP); TT := KT!X;
pt0 := PC ! [-22/9, 7/27, 1];
prm := Parametrization(PC, Place(pt0));
eqs := DefiningEquations(prm);
at := Evaluate(eqs[1],[TT,KT!1])/Evaluate(eqs[3],[TT,KT!1]);
bt := Evaluate(eqs[2],[TT,KT!1])/Evaluate(eqs[3],[TT,KT!1]);
printf "a(T) = %o\nb(T) = %o\n", at, bt;
assert Evaluate(Fab,[at,bt]) eq 0;
printf "poles of a(T): %o\n", [ r[1] : r in Roots(PP!Denominator(at)) ];
printf "poles of b(T): %o\n", [ r[1] : r in Roots(PP!Denominator(bt)) ];
printf "zeros of b(T): %o\n", [ r[1] : r in Roots(PP!Numerator(bt)) ];

printf "\n---- the T=3 fibre of X (where X -> H is undefined) ----\n";
G1 := X*(X-3)*(X^2-24*X+36);
G2 := (X-3)*(7*X-12)*(37*X^2-96*X+36);
GH := X*(7*X-12)*(X^2-24*X+36)*(37*X^2-96*X+36);
printf "G1*G2 = (X-3)^2 * GH ? %o\n", G1*G2 eq (X-3)^2*GH;
u1 := Evaluate(ExactQuotient(G1, X-3), 3);
u2 := Evaluate(ExactQuotient(G2, X-3), 3);
printf "G1/(X-3) at X=3 = %o (square %o) ; G2/(X-3) at X=3 = %o (square %o)\n",
   u1, IsSquare(u1), u2, IsSquare(u2);
printf "  -> the two places of X over T=3 have residue field Q(sqrt(%o)) : rational = %o\n",
   u1, IsSquare(u1);
printf "GH(3) = %o (square %o) -> H has no rational point over X=3\n", Evaluate(GH,3), IsSquare(Evaluate(GH,3));
printf "  (so the affine singular point (T,u,v)=(3,0,0) is the ONLY rational point of the\n";
printf "   affine model X that does not map to a rational point of H, and it was tested.)\n";

printf "\n---- independent recount: all T of height <= 3000 with G1,G2 both squares ----\n";
Z := Integers();
c1 := [ Z!Coefficient(G1,i) : i in [0..4] ];
c2 := [ Z!Coefficient(G2,i) : i in [0..4] ];
found := [];
BH := 3000;
for q in [1..BH] do
  q2v := q*q; q3v := q2v*q; q4v := q3v*q;
  for p in [-BH..BH] do
    if GCD(Abs(p),q) ne 1 then continue; end if;
    p2v := p*p; p3v := p2v*p; p4v := p3v*p;
    n1 := c1[1]*q4v + c1[2]*p*q3v + c1[3]*p2v*q2v + c1[4]*p3v*q + c1[5]*p4v;
    if n1 lt 0 or not IsSquare(n1) then continue; end if;
    n2 := c2[1]*q4v + c2[2]*p*q3v + c2[3]*p2v*q2v + c2[4]*p3v*q + c2[5]*p4v;
    if n2 lt 0 or not IsSquare(n2) then continue; end if;
    Append(~found, p/q);
  end for;
end for;
printf "T with BOTH G1,G2 squares, |p|,q <= %o : %o\n", BH, found;
printf "GAP_DONE\n";
quit;
