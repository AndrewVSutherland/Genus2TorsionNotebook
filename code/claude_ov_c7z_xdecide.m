// ===========================================================================
// Lane 2 (claude_ov_c7z) -- DECIDE the order-112 locus on the contact-7 c=2 slice.
//
// Chain (re-derived from scratch in Parts 0-2 below; nothing is taken on trust):
//   contact-7 chart   h = 1 - (7/2)x + A x^2 + B x^3,  f = (h^2 + (x-1)^7)/x^2
//   rational root of f  <->  r = 1 - v^2,  v a root of
//        Q(v) = (v+1)^2 (c4 v^2 + c0) + v^5 - v^3 - v^2/2 ,  c4 = B+2, c0 = 5/2-(A+B).
//   z_i := 1/(1+v_i)  ==>  P(z) = prod(z-z_i) = z^5 - 6z^4 + e2 z^3 - e3 z^2 + e4 z - 2
//        with the SINGLE remaining universal relation  e3 = 2 e4 - 2.
//   k rational z_i  <->  quintic factor type  <->  2-rank;  k=5  <->  [2,2,2,14], order 112.
//   SLICE z_1 = 2:  the other four are roots of the pencil
//        P_m(w) = w^4 - 4 w^3 + m w^2 - (2/3) m w + 1.
//   Split into two rational quadratics (w^2-aw+b)(w^2-(4-a)w+1/b):
//        Cab :  -3ab^2 + 10b^2 + 2a^2 b - 8ab + 3a - 2 = 0        (genus 0)
//   ORDER-112 locus of the slice  =  X : u^2 = a^2-4b , v^2 = (4-a)^2-4/b  over Cab.
// ===========================================================================
SetColumns(0);
if not assigned MemGB then MemGB := 12; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned BH then BH := 1200; elif Type(BH) eq MonStgElt then BH := StringToInteger(BH); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals();
PP<X> := PolynomialRing(Q);
KT := FieldOfFractions(PP);
TT := KT ! X;

printf "======== PART 0: universal relations, re-derived symbolically ========\n";
Rc<c4,c0> := PolynomialRing(Q,2);
Pv<v> := PolynomialRing(Rc);
Qv := (v+1)^2*(c4*v^2+c0) + v^5 - v^3 - v^2/2;
printf "Q(v) = %o\n", Qv;
printf "Q(-1) = %o   Q'(-1) = %o\n", Evaluate(Qv,-1), Evaluate(Derivative(Qv),-1);
Fc := FieldOfFractions(Rc);
Pz<z> := PolynomialRing(Fc);
Qsub := &+[ (Fc!Coefficient(Qv,i)) * (1/z - 1)^i : i in [0..5] ];
Pzz := Pz ! ( -2 * z^5 * Qsub );
printf "P(z) = %o\n", Pzz;
cf := [ Coefficient(Pzz,i) : i in [0..5] ];
e1 := -cf[5]; e2 := cf[4]; e3 := -cf[3]; e4 := cf[2]; e5 := -cf[1];
printf "monic=%o | e1=%o e5=%o | RELATION e3-2e4+2 = %o\n", cf[6] eq 1, e1, e5, e3-2*e4+2;
assert cf[6] eq 1 and e1 eq 6 and e5 eq 2 and (e3-2*e4+2) eq 0;
printf "PART0_OK: universal relations are  e1=6, e5=2, e3=2e4-2  (2 free params e2,e4)\n";

printf "\n======== PART 1: the c=2 slice and the curve Cab ========\n";
R2<aa,bb> := PolynomialRing(Q,2);
Fab := 3*(aa + (4-aa)*bb^2) - 2*(bb^2 + 1 + 4*aa*bb - aa^2*bb);
printf "Cab : %o = 0\n", Fab;
assert Fab eq -3*aa*bb^2 + 10*bb^2 + 2*aa^2*bb - 8*aa*bb + 3*aa - 2;
A2 := AffineSpace(R2);
Cab := Curve(A2, Fab);
PCab := ProjectiveClosure(Cab);
printf "Cab irreducible = %o, genus = %o\n", IsIrreducible(Fab), Genus(PCab);
// hand cross-check: solving Fab=0 for aa gives disc_a = (b-1)^2 (9b^2-14b+9)
disc_a := (-3*bb^2 - 8*bb + 3)^2 - 4*(2*bb)*(10*bb^2 - 2);
printf "disc_a(Cab as quadratic in a) = %o ; equals (b-1)^2(9b^2-14b+9)? %o\n",
    disc_a, disc_a eq (bb-1)^2*(9*bb^2-14*bb+9);

printf "\n======== PART 2: parametrize Cab and build g1, g2 ========\n";
pt0 := PCab ! [-22/9, 7/27, 1];
prm := Parametrization(PCab, Place(pt0));
eqs := DefiningEquations(prm);
at := Evaluate(eqs[1],[TT,KT!1])/Evaluate(eqs[3],[TT,KT!1]);
bt := Evaluate(eqs[2],[TT,KT!1])/Evaluate(eqs[3],[TT,KT!1]);
printf "a(T) = %o\nb(T) = %o\n", at, bt;
assert Evaluate(Fab,[at,bt]) eq 0;
d1 := at^2 - 4*bt;
d2 := (4-at)^2 - 4/bt;
printf "d1(T) = %o\nd2(T) = %o\n", d1, d2;

// square-class extraction WITHOUT losing the unit constant (Factorization unit trap)
function SqClass(fn)
  Nn := PP ! Numerator(fn); Dd := PP ! Denominator(fn);
  G := Nn*Dd;
  fac := Factorisation(G);
  unit := G div &*[ t[1]^t[2] : t in fac ];
  assert unit * &*[ t[1]^t[2] : t in fac ] eq G;
  sf := unit;
  for t in fac do if IsOdd(t[2]) then sf *:= t[1]; end if; end for;
  cofac := ExactQuotient(G, sf);
  ok := IsSquare(cofac);
  assert ok;
  return sf;
end function;

G1 := SqClass(d1);
G2 := SqClass(d2);
printf "G1 = %o\n", G1;
printf "G2 = %o\n", G2;
printf "G1 factorisation: %o   lc=%o (square %o)\n", Factorisation(G1), LeadingCoefficient(G1), IsSquare(LeadingCoefficient(G1));
printf "G2 factorisation: %o   lc=%o (square %o)\n", Factorisation(G2), LeadingCoefficient(G2), IsSquare(LeadingCoefficient(G2));
printf "gcd(G1,G2) = %o\n", GCD(G1,G2);
nchk := 0;
for tv in [Q| 1/2, 5, -7, 11/3, 100, -4/9, 17, 23/5 ] do
  q1 := Evaluate(d1,tv)/Evaluate(G1,tv);
  q2 := Evaluate(d2,tv)/Evaluate(G2,tv);
  assert IsSquare(q1) and IsSquare(q2);
  nchk +:= 1;
end for;
printf "square-class of d1,d2 verified at %o sample T values\n", nchk;

// clean integral representatives (square class unchanged)
function CleanForm(G)
  d := LCM([ Denominator(c) : c in Coefficients(G) ]);
  Gi := PP ! (d*G);
  cnt := GCD([ Numerator(c) : c in Coefficients(Gi) ]);
  Gi := PP ! (Gi/cnt);
  // G = Gi * cnt/d ; adjust class by the squarefree part of cnt*d
  s := SquarefreeFactorization(cnt*d);
  return s*Gi;
end function;
G1c := CleanForm(G1); G2c := CleanForm(G2);
printf "clean G1 = %o\nclean G2 = %o\n", G1c, G2c;
assert IsSquare(KT!(G1c/G1)) and IsSquare(KT!(G2c/G2));

printf "\n======== PART 3: genus of X : u^2=G1(T), v^2=G2(T) ========\n";
KA := FunctionField(Q); TA := KA.1;
PA<YA> := PolynomialRing(KA);
F1 := FunctionField(YA^2 - Evaluate(G1c,TA));
printf "  genus(u^2=G1) = %o   (deg over Q(T) = %o)\n", Genus(F1), Degree(F1);
PB<YB> := PolynomialRing(F1);
F2 := FunctionField(YB^2 - (F1!Evaluate(G2c,TA)));
gX := Genus(F2);
printf "  GENUS(X) via function-field tower = %o  (deg over Q(T) = %o)\n", gX, Degree(F2)*Degree(F1);

printf "\n======== PART 4: the three (Z/2)^2 quotients ========\n";
for pr in [ <G1c,"E1 (u^2=G1)">, <G2c,"E2 (v^2=G2)"> ] do
  Hc := HyperellipticCurve(pr[1]);
  Em := MinimalModel(EllipticCurve(Hc));
  Gm, mm, p1, p2 := MordellWeilGroup(Em);
  printf "%o  genus %o -> %o  cond %o  MW %o proven=(%o,%o)\n",
     pr[2], Genus(Hc), Em, Conductor(Em), Invariants(Gm), p1, p2;
end for;
GH := CleanForm(SqClass(KT!(G1*G2)));
printf "H : y^2 = %o   (deg %o)\n", GH, Degree(GH);
printf "H factorisation: %o  lc=%o\n", Factorisation(GH), LeadingCoefficient(GH);
CH := HyperellipticCurve(GH);
printf "H genus = %o ; GH squarefree = %o\n", Genus(CH), IsSquarefree(GH);
printf "H discriminant = %o\n", Factorisation(Numerator(Q!Discriminant(CH)));
JH := Jacobian(CH);
tS := TorsionSubgroup(JH);
printf "H : TorsionSubgroup(Jac) = %o\n", Invariants(tS);
tt := Cputime();
lo, hi := RankBounds(JH);
printf "H : RankBounds(Jac) = %o .. %o   [%o s]\n", lo, hi, Cputime(tt);

printf "\n======== PART 5: rational points ========\n";
ptsH := Points(CH : Bound := 10000);
printf "H : Points(Bound:=10000) = %o\n", ptsH;

// direct search on X : both G1c(p/q)*q^4 and G2c(p/q)*q^4 squares (integer arithmetic)
Z := Integers();
d1i := LCM([ Denominator(c) : c in Coefficients(G1c) ]);
d2i := LCM([ Denominator(c) : c in Coefficients(G2c) ]);
c1 := [ Z!(d1i*Coefficient(G1c,i)) : i in [0..4] ];
c2 := [ Z!(d2i*Coefficient(G2c,i)) : i in [0..4] ];
printf "integral forms: d1i=%o c1=%o ; d2i=%o c2=%o\n", d1i, c1, d2i, c2;
found := [];
for q in [1..BH] do
  q2v := q*q; q3v := q2v*q; q4v := q3v*q;
  for p in [-BH..BH] do
    if GCD(Abs(p),q) ne 1 then continue; end if;
    p2v := p*p; p3v := p2v*p; p4v := p3v*p;
    n1 := d1i*(c1[1]*q4v + c1[2]*p*q3v + c1[3]*p2v*q2v + c1[4]*p3v*q + c1[5]*p4v);
    if n1 lt 0 or not IsSquare(n1) then continue; end if;
    n2 := d2i*(c2[1]*q4v + c2[2]*p*q3v + c2[3]*p2v*q2v + c2[4]*p3v*q + c2[5]*p4v);
    if n2 lt 0 or not IsSquare(n2) then continue; end if;
    Append(~found, p/q);
  end for;
end for;
printf "X : direct search |p|<=%o, 1<=q<=%o  ->  T with BOTH squares: %o\n", BH, BH, found;
printf "XDECIDE_DONE\n";
quit;
