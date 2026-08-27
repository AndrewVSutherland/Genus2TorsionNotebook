// claude_ov_b4_dump.m -- Lane 4 (route B4): export the raise-locus functions.
//
// On the (6,2) component K = Q(E), E: Y^2 = -u^3+3u^2-2u+1 (cond 92, rank 1),
// every member's Richelot image has factor type [2,4] generically, and the
// image's 2-rank rises to 2 exactly on the double cover  s^2 = A1  (image
// quartic splits into two rational quadratics), which has genus 3.
// A second, weaker locus is  s^2 = d1  (image quadratic splits, type [1,1,4],
// 2-rank stays 1) -- also genus 3.
//
// Since E(Q) = Z.G with G = (u,Y) = (0,1), the rational points of BOTH covers
// lie over {nG}, so testing "A1(nG) is a rational square" for |n| <= N is an
// EXHAUSTIVE search of the cover in that range.  This script dumps A1 and d1
// in the form c0(u) + c1(u)*Y so a cheap mod-p Legendre sieve can run over n.
//
// Run: magma -b claude_ov_b4_dump.m

SetColumns(0);
if not assigned MemGB then MemGB := 24; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
Fu<u> := FunctionField(Q);
PY<Y> := PolynomialRing(Fu);
K<yy> := ext<Fu | Y^2 - (-u^3+3*u^2-2*u+1)>;

uu := K!u;
a2f := (uu-2)^4;
b1f := 2*(uu^5-5*uu^4+9*uu^3-8*uu^2+4*uu-2);
c0f := uu^2*(uu^2-uu+1)^2;
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

facfp := Factorization(fp);
quarts := [t[1] : t in facfp | Degree(t[1]) eq 4];
quads  := [t[1] : t in facfp | Degree(t[1]) eq 2];
q4p := quarts[1]/LeadingCoefficient(quarts[1]);
s1, p1, q1, r1, Res1 := ResData(q4p);
A1 := Roots(Res1)[1][1] - p1;
h1 := quads[1];
d1 := Coefficient(h1,1)^2 - 4*Coefficient(h1,2)*Coefficient(h1,0);

// genus of the "image quartic has a rational root" cover (degree-4 over E)
if IsIrreducible(q4p) then
    L4 := ext<K | q4p>;
    printf "GENUS of {image quartic has a rational root} cover = %o\n", Genus(L4);
end if;
// genus of the "seed quartic has a rational root" cover, for the record
if IsIrreducible(Q4) then
    L5 := ext<K | Q4>;
    printf "GENUS of {SEED quartic has a rational root} cover = %o\n", Genus(L5);
end if;

procedure Dump(name, el)
    e := Eltseq(el);         // el = e[1] + e[2]*yy
    printf "%o_C0_NUM %o\n", name, Numerator(e[1]);
    printf "%o_C0_DEN %o\n", name, Denominator(e[1]);
    printf "%o_C1_NUM %o\n", name, Numerator(e[2]);
    printf "%o_C1_DEN %o\n", name, Denominator(e[2]);
end procedure;

Dump("A", A);
Dump("A1", A1);
Dump("D1", d1);
printf "DUMP_DONE\n";
quit;
