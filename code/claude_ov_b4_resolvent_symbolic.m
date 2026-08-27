// claude_ov_b4_resolvent_symbolic.m  -- Lane 4 (route B4), THE CRUX COMPUTATION.
//
// Over the (6,2) component of the quadratic-factor incidence of Flynn's
// order-11 family, decide whether the quartic Weierstrass cofactor Q4 has
// Galois group inside D4, i.e. whether its resolvent cubic has a root in the
// function field of the component.  If it does IDENTICALLY, every member of
// the stream carries a Galois-stable Richelot (2,2)-kernel over Q and the
// seed supply for route B4 is INFINITE (the component is the rank-1
// conductor-92 elliptic curve).  If not, the condition cuts out a 3:1 cover
// whose genus we then compute.
//
// Component:  a2(u) v^2 + b1(u) v + c0(u) = 0, smooth model
//             E: Y^2 = -u^3+3u^2-2u+1  (= w^3+3w^2+2w+1 at w=-u; cond 92, rank 1)
//             v = (-b1(u) + 4(u-1)Y) / (2 a2(u)).
//
// Run:  magma -b claude_ov_b4_resolvent_symbolic.m

SetColumns(0);
if not assigned MemGB then MemGB := 8; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
Fu<u> := FunctionField(Q);
PY<Y> := PolynomialRing(Fu);
K<yy> := ext<Fu | Y^2 - (-u^3+3*u^2-2*u+1)>;
printf "FUNCTION FIELD K over Q(u): genus=%o\n", Genus(K);

uu := K!u;
a2 := (uu-2)^4;
b1 := 2*(uu^5-5*uu^4+9*uu^3-8*uu^2+4*uu-2);
c0 := uu^2*(uu^2-uu+1)^2;
vv := (-b1 + 4*(uu-1)*yy) / (2*a2);
printf "COMPONENT CHECK a2*v^2+b1*v+c0 = %o\n", a2*vv^2+b1*vv+c0;

// --- recover t as an element of K ---------------------------------------
PT<T> := PolynomialRing(K);
PX<x> := PolynomialRing(PolynomialRing(K));   // placeholder, rebuilt below

// F_t(x) with t symbolic: work in K[T][x]
RT<Tt> := PolynomialRing(K);
RTX<xx> := PolynomialRing(RT);
Ft := xx^6 + 2*xx^5 + (2*Tt+3)*xx^4 + 2*xx^3 + (Tt^2+1)*xx^2 + 2*Tt*(1-Tt)*xx + Tt^2;
qd := xx^2 + (RT!uu)*xx + (RT!vv);
_, rem := Quotrem(Ft, qd);
e1 := Coefficient(rem, 1);
e0 := Coefficient(rem, 0);
printf "deg_t e1 = %o, deg_t e0 = %o\n", Degree(e1), Degree(e0);
g := GCD(e1, e0);
printf "GCD degree in t = %o\n", Degree(g);
assert Degree(g) eq 1;
tK := -Coefficient(g,0)/Coefficient(g,1);
printf "t recovered in K (degree data): ok\n";

// --- build Q4 over K -----------------------------------------------------
PK<X> := PolynomialRing(K);
F := X^6 + 2*X^5 + (2*tK+3)*X^4 + 2*X^3 + (tK^2+1)*X^2 + 2*tK*(1-tK)*X + tK^2;
qd2 := X^2 + uu*X + vv;
Q4, r2 := Quotrem(F, qd2);
printf "DIVIDES EXACTLY: %o\n", r2 eq 0;
printf "Q4 degree = %o, monic = %o\n", Degree(Q4), LeadingCoefficient(Q4) eq 1;

// factor type of Q4 over K (i.e. generic factor type of the stream)
facQ4 := Factorization(Q4);
printf "Q4 FACTOR DEGREES OVER K: %o\n", [Degree(t[1]) : t in facQ4];

// --- resolvent cubic -----------------------------------------------------
a3 := Coefficient(Q4,3); a2c := Coefficient(Q4,2); a1 := Coefficient(Q4,1); a0 := Coefficient(Q4,0);
p := a2c - 3*a3^2/8;
q := a1 - a3*a2c/2 + a3^3/8;
r := a0 - a1*a3/4 + a2c*a3^2/16 - 3*a3^4/256;
PZ<Z> := PolynomialRing(K);
Res := Z^3 - p*Z^2 - 4*r*Z + (4*p*r - q^2);
facR := Factorization(Res);
printf "RESOLVENT FACTOR DEGREES OVER K: %o\n", [t[1] : t in facR] ne [] select [Degree(t[1]) : t in facR] else [];

rts := Roots(Res);
printf "NUMBER OF ROOTS OF RESOLVENT IN K: %o\n", #rts;
if #rts gt 0 then
    z0 := rts[1][1];
    printf "IDENTITY: the resolvent cubic HAS a root in K.\n";
    printf "CHECK Res(z0) = %o\n", Evaluate(Res, z0);
    printf "SEED SUPPLY INFINITE: every (6,2)-member has Galois group of Q4 inside D4.\n";
else
    printf "NO IDENTITY: the D4 condition is a proper 3:1 cover of E.\n";
end if;

printf "SYMBOLIC_DONE\n";
quit;
