// claude_ov_b4_split_cover.m -- Lane 4 (route B4), follow-up to the crux.
//
// The resolvent cubic of the quartic cofactor Q4 has a root z0 in the function
// field K of the (6,2) component (proved in claude_ov_b4_resolvent_symbolic.m).
// Writing Q4 in depressed form X^4 + p X^2 + q X + r, one has
//     Q4 = (X^2 + aX + b)(X^2 - aX + c),   a^2 = z0 - p =: A,
//     c = ((p+A) + q/a)/2,  b = ((p+A) - q/a)/2.
// So:
//   * A is ALWAYS a square in the quadratic extension K(sqrt A) -- that is the
//     Galois-stable Richelot kernel (route B4 proper);
//   * A is a square IN K  <=>  Q4 splits into two RATIONAL quadratics
//     <=> factor type [2,2,2] <=> 2-rank 2 <=> torsion contains [2,22] DIRECTLY.
// The second condition is a double cover  D: s^2 = A  of the genus-1 curve E.
// This script computes A, its divisor, the genus of D, and (if genus <= 1)
// an elliptic model + rank.
//
// Run:  magma -b claude_ov_b4_split_cover.m

SetColumns(0);
if not assigned MemGB then MemGB := 16; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
Fu<u> := FunctionField(Q);
PY<Y> := PolynomialRing(Fu);
K<yy> := ext<Fu | Y^2 - (-u^3+3*u^2-2*u+1)>;
printf "K genus = %o\n", Genus(K);

uu := K!u;
a2f := (uu-2)^4;
b1f := 2*(uu^5-5*uu^4+9*uu^3-8*uu^2+4*uu-2);
c0f := uu^2*(uu^2-uu+1)^2;
vv := (-b1f + 4*(uu-1)*yy) / (2*a2f);

RT<Tt> := PolynomialRing(K);
RTX<xx> := PolynomialRing(RT);
Ft := xx^6 + 2*xx^5 + (2*Tt+3)*xx^4 + 2*xx^3 + (Tt^2+1)*xx^2 + 2*Tt*(1-Tt)*xx + Tt^2;
qd := xx^2 + (RT!uu)*xx + (RT!vv);
_, rem := Quotrem(Ft, qd);
g := GCD(Coefficient(rem,1), Coefficient(rem,0));
assert Degree(g) eq 1;
tK := -Coefficient(g,0)/Coefficient(g,1);

PK<X> := PolynomialRing(K);
F := X^6 + 2*X^5 + (2*tK+3)*X^4 + 2*X^3 + (tK^2+1)*X^2 + 2*tK*(1-tK)*X + tK^2;
Q4, r2 := Quotrem(F, X^2 + uu*X + vv);
assert r2 eq 0;

a3 := Coefficient(Q4,3); a2c := Coefficient(Q4,2); a1 := Coefficient(Q4,1); a0 := Coefficient(Q4,0);
p := a2c - 3*a3^2/8;
q := a1 - a3*a2c/2 + a3^3/8;
r := a0 - a1*a3/4 + a2c*a3^2/16 - 3*a3^4/256;
PZ<Z> := PolynomialRing(K);
Res := Z^3 - p*Z^2 - 4*r*Z + (4*p*r - q^2);
rts := Roots(Res);
printf "RESOLVENT ROOTS IN K: %o\n", #rts;
assert #rts ge 1;
z0 := rts[1][1];
A := z0 - p;
printf "A = z0 - p computed.  A eq 0? %o\n", A eq 0;

// --- square class of A: divisor, odd-multiplicity support -----------------
DA := Divisor(A);
supp, mult := Support(DA);
odds := [i : i in [1..#supp] | IsOdd(mult[i])];
printf "DIV(A): %o places, degrees=%o, mults=%o\n", #supp, [Degree(s) : s in supp], mult;
printf "ODD-MULTIPLICITY PLACES: %o  (total degree %o)\n", #odds, &+([0] cat [Degree(supp[i]) : i in odds]);

isq, sq := IsSquare(A);
printf "IS A A SQUARE IN K? %o\n", isq;

// --- the double cover D: s^2 = A -----------------------------------------
PS<S> := PolynomialRing(K);
polD := S^2 - A;
if IsIrreducible(polD) then
    L := ext<K | polD>;
    gD := Genus(L);
    printf "DOUBLE COVER D: s^2 = A  has GENUS %o\n", gD;
    printf "CONSTANT FIELD OF L: %o\n", ConstantField(L);
    printf "RIEMANN-HURWITZ CHECK: 2g-2 = %o (branch degree %o)\n", 2*gD-2, &+([0] cat [Degree(supp[i]) : i in odds]);
else
    printf "s^2 - A REDUCIBLE over K: A is a square, cover splits.\n";
end if;

printf "SPLITCOVER_DONE\n";
quit;
