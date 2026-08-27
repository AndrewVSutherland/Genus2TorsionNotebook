// claude_ov_b4_codomain_symbolic.m -- Lane 4 (route B4), the decisive
// symbolic analysis of the Richelot IMAGE family.
//
// Setting: K = function field of the (6,2) component of the quadratic-factor
// incidence of Flynn's order-11 family (= the rank-1 conductor-92 elliptic
// curve E: Y^2 = -u^3+3u^2-2u+1).  Over K the sextic is F = G1 * Q4 with
// G1 = x^2+ux+v rational and Q4 an irreducible quartic whose resolvent cubic
// has a root z0 in K (proved in claude_ov_b4_resolvent_symbolic.m).  Hence
// Q4 = G2*G3 over L = K(sqrt A), A = z0 - p, with G2,G3 conjugate: the
// partition {G1,G2,G3} is GALOIS-STABLE and the Richelot isogeny is defined
// over Q for EVERY member.
//
// This script computes the Richelot codomain
//        f' = -H1 * H2 * H3 / Delta,   H_i = G_j' G_k - G_j G_k'
// (NB: Delta*H1H2H3 is the quadratic TWIST BY A of the correct model, and A is
//  a non-square exactly because the kernel is Galois-stable but not rational;
//  the sign is the H_i cyclic-order convention.  Both were pinned against
//  Magma's RichelotIsogenousSurfaces and against isogeny-invariance of
//  #J(F_p) -- see code/claude_ov_b4_validate.m.  Factor DEGREES are unaffected
//  by the twist, so the genus conclusions below do not depend on it.)
// as an element of K[x], factors it over K (=> the GENERIC factor type, hence
// the generic 2-rank, of the whole image family), and -- if that type is
// [2,4] -- computes the cover of E cut out by "the image quartic splits into
// two rational quadratics", i.e. the locus where the 2-rank RISES to 2 and the
// torsion becomes [2,22].  Its genus decides whether the image family contains
// infinitely many or only finitely many generic [2,22] curves.
//
// Run:  magma -b claude_ov_b4_codomain_symbolic.m

SetColumns(0);
if not assigned MemGB then MemGB := 24; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
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
_, rem := Quotrem(Ft, xx^2 + (RT!uu)*xx + (RT!vv));
gg := GCD(Coefficient(rem,1), Coefficient(rem,0));
assert Degree(gg) eq 1;
tK := -Coefficient(gg,0)/Coefficient(gg,1);

PK<X> := PolynomialRing(K);
F := X^6 + 2*X^5 + (2*tK+3)*X^4 + 2*X^3 + (tK^2+1)*X^2 + 2*tK*(1-tK)*X + tK^2;
G1 := X^2 + uu*X + vv;
Q4, r2 := Quotrem(F, G1);
assert r2 eq 0;
printf "STEP1 ok: F = G1 * Q4 over K\n";

// ---- helper: resolvent cubic data of a monic quartic over a field ---------
function ResData(q4, PZ)
    a3 := Coefficient(q4,3); a2c := Coefficient(q4,2);
    a1 := Coefficient(q4,1); a0 := Coefficient(q4,0);
    s := a3/4;
    p := a2c - 3*a3^2/8;
    q := a1 - a3*a2c/2 + a3^3/8;
    r := a0 - a1*a3/4 + a2c*a3^2/16 - 3*a3^4/256;
    Z := PZ.1;
    return s, p, q, r, Z^3 - p*Z^2 - 4*r*Z + (4*p*r - q^2);
end function;

PZ<Z> := PolynomialRing(K);
s0, p0, q0, r0, Res0 := ResData(Q4, PZ);
rts0 := Roots(Res0);
printf "STEP2 resolvent roots of Q4 in K: %o\n", #rts0;
assert #rts0 ge 1;
A := rts0[1][1] - p0;
printf "STEP2 A ne 0: %o\n", A ne 0;

// ---- the quadratic extension L = K(sqrt A) -------------------------------
PS<S> := PolynomialRing(K);
assert IsIrreducible(S^2 - A);
L<aL> := ext<K | S^2 - A>;
printf "STEP3 L = K(sqrt A) built; genus(L) = %o\n", Genus(L);

PL<XL> := PolynomialRing(L);
sL := L!s0; pL := L!p0; qL := L!q0; AL := L!A;
cc := ((pL + AL) + qL/aL)/2;
bb := ((pL + AL) - qL/aL)/2;
G1L := PL![L!Coefficient(G1,i) : i in [0..2]];
G2L := (XL + sL)^2 + aL*(XL + sL) + bb;
G3L := (XL + sL)^2 - aL*(XL + sL) + cc;
FL  := PL![L!Coefficient(F,i) : i in [0..6]];
printf "STEP4 G1*G2*G3 = F over L: %o\n", G1L*G2L*G3L eq FL;

// ---- Richelot ------------------------------------------------------------
M := Matrix(L, 3, 3, [Coefficient(G1L,2), Coefficient(G1L,1), Coefficient(G1L,0),
                      Coefficient(G2L,2), Coefficient(G2L,1), Coefficient(G2L,0),
                      Coefficient(G3L,2), Coefficient(G3L,1), Coefficient(G3L,0)]);
Delta := Determinant(M);
printf "STEP5 Delta ne 0: %o\n", Delta ne 0;
H1 := Derivative(G2L)*G3L - G2L*Derivative(G3L);
H2 := Derivative(G3L)*G1L - G3L*Derivative(G1L);
H3 := Derivative(G1L)*G2L - G1L*Derivative(G2L);
fpL := -H1*H2*H3/Delta;
printf "STEP5 codomain degree = %o\n", Degree(fpL);

// descend to K
coefsK := [];
descends := true;
for i in [0..Degree(fpL)] do
    c := Coefficient(fpL, i);
    ok, cK := IsCoercible(K, c);
    if not ok then descends := false; break; end if;
    Append(~coefsK, cK);
end for;
printf "STEP6 codomain descends to K: %o\n", descends;
assert descends;
fp := PK!coefsK;

// ---- generic factor type of the image family ------------------------------
facfp := Factorization(fp);
printf "STEP7 CODOMAIN FACTOR DEGREES OVER K: %o\n", Sort([Degree(t[1]) : t in facfp]);
printf "STEP7 multiplicities: %o\n", [t[2] : t in facfp];

degs := Sort([Degree(t[1]) : t in facfp]);
// 2-rank from the factor type (even-degree subsets modulo {empty, full})
neven := 0;
k := #degs;
for m in [0..2^k-1] do
    sm := &+([0] cat [degs[i] : i in [1..k] | ((m div 2^(i-1)) mod 2) eq 1]);
    if IsEven(sm) then neven +:= 1; end if;
end for;
printf "STEP7 GENERIC 2-RANK OF THE IMAGE FAMILY = %o\n", Ilog(2, neven) - 1;

// ---- the raise locus: does the image quartic split? -----------------------
quarts := [t[1] : t in facfp | Degree(t[1]) eq 4];
if #quarts eq 0 then
    printf "STEP8 no degree-4 factor: the image family already has 2-rank >= 2 generically!\n";
else
    q4p := quarts[1];
    q4p := q4p / LeadingCoefficient(q4p);
    s1, p1, q1, r1, Res1 := ResData(q4p, PZ);
    rts1 := Roots(Res1);
    printf "STEP8 image quartic resolvent roots in K: %o\n", #rts1;
    if #rts1 eq 0 then
        printf "STEP8 image quartic has NO K-rational resolvent root: the image family\n";
        printf "STEP8 is NOT further Richelot-walkable by the q-block, and the split\n";
        printf "STEP8 condition is a 3:1 cover.\n";
    else
        A1 := rts1[1][1] - p1;
        printf "STEP8 A1 ne 0: %o\n", A1 ne 0;
        isq, _ := IsSquare(A1);
        printf "STEP8 IS A1 A SQUARE IN K (=> image ALWAYS 2-rank 2)? %o\n", isq;
        if not isq then
            DA1 := Divisor(A1);
            supp, mult := Support(DA1);
            odds := [i : i in [1..#supp] | IsOdd(mult[i])];
            bdeg := &+([0] cat [Degree(supp[i]) : i in odds]);
            printf "STEP8 branch divisor degree of s^2 = A1 : %o\n", bdeg;
            PS1<S1> := PolynomialRing(K);
            if IsIrreducible(S1^2 - A1) then
                L1 := ext<K | S1^2 - A1>;
                printf "STEP8 RAISE-LOCUS COVER  s^2 = A1  HAS GENUS %o\n", Genus(L1);
                printf "STEP8 constant field: %o\n", ConstantField(L1);
            else
                printf "STEP8 s^2 - A1 reducible\n";
            end if;
        end if;
    end if;
end if;

// ---- the other raise route: does the image QUADRATIC block split? ---------
quads := [t[1] : t in facfp | Degree(t[1]) eq 2];
if #quads gt 0 then
    h1 := quads[1];
    d1 := Coefficient(h1,1)^2 - 4*Coefficient(h1,2)*Coefficient(h1,0);
    isq1 := IsSquare(d1);
    printf "STEP9 IS disc(image quadratic) A SQUARE IN K? %o\n", isq1;
    if not isq1 then
        PS2<S2> := PolynomialRing(K);
        if IsIrreducible(S2^2 - K!d1) then
            L2 := ext<K | S2^2 - K!d1>;
            printf "STEP9 COVER {image quadratic splits} HAS GENUS %o\n", Genus(L2);
        end if;
    end if;
end if;

printf "CODOMAIN_SYMBOLIC_DONE\n";
quit;
