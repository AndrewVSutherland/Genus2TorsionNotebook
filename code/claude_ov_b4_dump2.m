// claude_ov_b4_dump2.m -- Lane 4 (route B4): dump the REMAINING raise-locus
// data that claude_ov_b4_dump.m did not export.
//
// Over K = Q(E), E: yy^2 = -u^3+3u^2-2u+1 (= 92.a1, rank 1, trivial torsion),
// the Richelot codomain of the (6,2)-stream member at (u,yy) is
//        f' = h1 * q4'   (generic factor type [2,4], 2-rank 1).
// A 2-rank RAISE needs one of
//   (a) q4' = product of two rational quadratics        <=> some resolvent
//       root z with z - p' a square; generically the unique rational root
//       gives A1 (already dumped).  An EXTRA rational resolvent root exists
//       only where disc(Res1) = disc(q4') is a square  -> dumped here as W.
//   (b) h1 splits AND q4' is reducible                  <=> D1 square (dumped)
// so {A1 square} u {W square} u {D1 square} is a NECESSARY superset for a
// 2-rank raise, and each is a double cover of E.
//
// Also dumped: q1 (the linear coefficient of the depressed q4'); the
// biquadratic degeneration q1 = 0 is the one case where the A1 criterion
// changes shape, so we need its (finite) zero locus.
//
// Run: nohup magma -b claude_ov_b4_dump2.m > results/claude_ov_b4_dump2.log 2>&1 &

SetColumns(0);
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
Fu<u> := FunctionField(Q);
PY<Y> := PolynomialRing(Fu);
K<yy> := ext<Fu | Y^2 - (-u^3+3*u^2-2*u+1)>;

uu := K!u;
a2f := (uu-2)^4;
b1f := 2*(uu^5-5*uu^4+9*uu^3-8*uu^2+4*uu-2);
vv := (-b1f + 4*(uu-1)*yy) / (2*a2f);

RT<Tt> := PolynomialRing(K);
RTX<xx> := PolynomialRing(RT);
Ft := xx^6 + 2*xx^5 + (2*Tt+3)*xx^4 + 2*xx^3 + (Tt^2+1)*xx^2 + 2*Tt*(1-Tt)*xx + Tt^2;
_, rem := Quotrem(Ft, xx^2 + (RT!uu)*xx + (RT!vv));
gg := GCD(Coefficient(rem,1), Coefficient(rem,0));
tK := -Coefficient(gg,0)/Coefficient(gg,1);

PK<X> := PolynomialRing(K);
F := X^6 + 2*X^5 + (2*tK+3)*X^4 + 2*X^3 + (tK^2+1)*X^2 + 2*tK*(1-tK)*X + tK^2;
G1 := X^2 + uu*X + vv;
Q4 := F div G1;
printf "SETUP ok, Q4 monic=%o\n", LeadingCoefficient(Q4) eq 1;

PZ<Z> := PolynomialRing(K);
function ResData(q4)
    a3 := Coefficient(q4,3); a2c := Coefficient(q4,2);
    a1 := Coefficient(q4,1); a0 := Coefficient(q4,0);
    s := a3/4; p := a2c - 3*a3^2/8;
    q := a1 - a3*a2c/2 + a3^3/8;
    r := a0 - a1*a3/4 + a2c*a3^2/16 - 3*a3^4/256;
    return s, p, q, r, Z^3 - p*Z^2 - 4*r*Z + (4*p*r - q^2);
end function;

s0, p0, q0, r0, Res0 := ResData(Q4);
A := Roots(Res0)[1][1] - p0;

PS<S> := PolynomialRing(K);
L<aL> := ext<K | S^2 - A>;
PL<XL> := PolynomialRing(L);
sL := L!s0; pL := L!p0; qL := L!q0; AL := L!A;
cc := ((pL + AL) + qL/aL)/2;
bb := ((pL + AL) - qL/aL)/2;
G1L := PL![L!Coefficient(G1,i) : i in [0..2]];
G2L := (XL + sL)^2 + aL*(XL + sL) + bb;
G3L := (XL + sL)^2 - aL*(XL + sL) + cc;
M := Matrix(L, 3, 3, [Coefficient(G1L,2), Coefficient(G1L,1), Coefficient(G1L,0),
                      Coefficient(G2L,2), Coefficient(G2L,1), Coefficient(G2L,0),
                      Coefficient(G3L,2), Coefficient(G3L,1), Coefficient(G3L,0)]);
Delta := Determinant(M);
H1 := Derivative(G2L)*G3L - G2L*Derivative(G3L);
H2 := Derivative(G3L)*G1L - G3L*Derivative(G1L);
H3 := Derivative(G1L)*G2L - G1L*Derivative(G2L);
fpL := -H1*H2*H3/Delta;
fp := PK![K!Coefficient(fpL,i) : i in [0..6]];
printf "CODOMAIN built\n";

facfp := Factorization(fp);
printf "CODOMAIN FACTOR DEGREES OVER K: %o\n", [Degree(t[1]) : t in facfp];
quarts := [t[1] : t in facfp | Degree(t[1]) eq 4];
q4p := quarts[1]/LeadingCoefficient(quarts[1]);
s1, p1, q1, r1, Res1 := ResData(q4p);
A1 := Roots(Res1)[1][1] - p1;
W := Discriminant(q4p);
printf "IS W=disc(q4') A SQUARE IN K? %o\n", IsSquare(W);
printf "IS q1 IDENTICALLY ZERO? %o\n", q1 eq 0;
LW := ext<K | S^2 - W>;
printf "GENUS of {disc(q4') is a square} cover = %o\n", Genus(LW);
printf "constant field of that cover: %o\n", ConstantField(LW);

procedure Dump(name, el)
    e := Eltseq(el);         // el = e[1] + e[2]*yy
    printf "%o_C0_NUM %o\n", name, Numerator(e[1]);
    printf "%o_C0_DEN %o\n", name, Denominator(e[1]);
    printf "%o_C1_NUM %o\n", name, Numerator(e[2]);
    printf "%o_C1_DEN %o\n", name, Denominator(e[2]);
end procedure;

Dump("W", W);
Dump("Q1", q1);
// square-class-reduced version of A1 for the record
printf "A1_EQ_CHECK %o\n", A1 eq Roots(Res1)[1][1] - p1;
printf "DUMP2_DONE\n";
quit;
