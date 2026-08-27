// claude_ov_b4_c43param.m -- Lane 4 (route B4): parametrize the NEW (4,3)
// component of the quadratic-factor incidence of Flynn's order-11 family and
// carry out over it the SAME symbolic analysis that closed the (6,2) stream.
//
//   C43 : -4v^3 + (u^2-8u+8)v^2 + (2u^3-6u^2+10u-4)v + (u^2-u+1)^2 = 0
// has GEOMETRIC GENUS 0 with rational points, hence is birational to P^1:
// an INFINITE, height-dense seed stream (~H^2 members of height <= H, versus
// the (6,2) stream's ~log H).
//
// Magma's Parametrization(C, Place(p), P1) was run once (base point
// (41/9 : 49/9 : 1), denominator (s-44)^2(s-68)^2) and the resulting degree-4
// map was reduced by s = (68w/3 - 44)/(w/3 - 1) to the CLEAN normal form
// hard-coded below:
//        u(w) = (2w^3 - 2w^2 + 2w - 1)/w^2
//        v(w) = ((w^2 - w + 1)/w)^2
//        t(w) = -((w-1)(w^2-w+1)/w^2)^2          <-- t is MINUS a square
// (all three verified here from scratch: C43(u,v) = 0 and x^2+ux+v | F_{t}).
//
// Run: code/claude_magma_slot.sh -b MemGB:=12 code/claude_ov_b4_c43param.m \
//        > results/claude_ov_b4_c43param.log 2>&1 &

SetColumns(0);
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
Fs<w> := FunctionField(Q);

uu := (2*w^3 - 2*w^2 + 2*w - 1)/w^2;
vv := ((w^2 - w + 1)/w)^2;
tS := -((w-1)*(w^2-w+1)/w^2)^2;

chk := -4*vv^3 + (uu^2-8*uu+8)*vv^2 + (2*uu^3-6*uu^2+10*uu-4)*vv + (uu^2-uu+1)^2;
printf "CHECK C43(u(w),v(w)) = %o\n", chk;
assert chk eq 0;

PK<X> := PolynomialRing(Fs);
F := X^6 + 2*X^5 + (2*tS+3)*X^4 + 2*X^3 + (tS^2+1)*X^2 + 2*tS*(1-tS)*X + tS^2;
G1 := X^2 + uu*X + vv;
Q4, r2 := Quotrem(F, G1);
printf "DIVIDES EXACTLY: %o\n", r2 eq 0;
assert r2 eq 0;
facF := Factorization(F);
printf "GENERIC FACTOR TYPE of F over Q(w): %o\n", [Degree(f[1]) : f in facF];
for f in facF do printf "  FACTOR (%o): %o\n", Degree(f[1]), f[1]; end for;
printf "Q4 factor degrees over Q(w): %o\n", [Degree(f[1]) : f in Factorization(Q4)];
printf "F squarefree: %o\n", Discriminant(F) ne 0;

PZ<Z> := PolynomialRing(Fs);
function ResData(q4)
    a3 := Coefficient(q4,3); a2c := Coefficient(q4,2);
    a1 := Coefficient(q4,1); a0 := Coefficient(q4,0);
    sc := a3/4; p := a2c - 3*a3^2/8;
    q := a1 - a3*a2c/2 + a3^3/8;
    r := a0 - a1*a3/4 + a2c*a3^2/16 - 3*a3^4/256;
    return sc, p, q, r, Z^3 - p*Z^2 - 4*r*Z + (4*p*r - q^2);
end function;

PS<Sv> := PolynomialRing(Fs);

procedure Dump(name, el)
    printf "%o_NUM %o\n", name, Numerator(el);
    printf "%o_DEN %o\n", name, Denominator(el);
end procedure;

Dump("C43_U", uu);  Dump("C43_V", vv);  Dump("C43_T", tS);

// --------------------------------------------------------------------------
// Richelot analysis, once for each rational quadratic block of F.
// --------------------------------------------------------------------------
blocks := [];
quadsF := [f[1] : f in facF | Degree(f[1]) eq 2 and f[2] eq 1];
linsF  := [f[1] : f in facF | Degree(f[1]) eq 1 and f[2] eq 1];
for b in quadsF do Append(~blocks, b); end for;
for i in [1..#linsF-1] do for j in [i+1..#linsF] do
    Append(~blocks, linsF[i]*linsF[j]);
end for; end for;
printf "NUMBER OF RATIONAL QUADRATIC BLOCKS over Q(w): %o\n", #blocks;

for bi in [1..#blocks] do
    B := blocks[bi]/LeadingCoefficient(blocks[bi]);
    R4, rr := Quotrem(F, B);
    assert rr eq 0;
    R4 := R4/LeadingCoefficient(R4);
    printf "\n=== BLOCK %o : %o ===\n", bi, B;
    s0, p0, q0, r0, Res0 := ResData(R4);
    rts0 := Roots(Res0);
    printf "  RESOLVENT ROOTS of the quartic cofactor in Q(w): %o\n", #rts0;
    if #rts0 eq 0 then
        printf "  NO IDENTITY: Galois(quartic) is not generically inside D4.\n";
        printf "  Res0 irreducible over Q(w): %o\n", IsIrreducible(Res0);
        LD := ext<Fs | Res0>;
        printf "  D4 COVER: degree %o, genus %o, constant field %o\n", Degree(LD), Genus(LD), ConstantField(LD);
        dd := Discriminant(Res0);
        b2, _ := IsSquare(dd);
        printf "  disc(Res0) square in Q(w)? %o\n", b2;
        for i in [0..3] do
            printf "  RES0_B%o_C%o_NUM %o\n", bi, i, Numerator(Coefficient(Res0,i));
            printf "  RES0_B%o_C%o_DEN %o\n", bi, i, Denominator(Coefficient(Res0,i));
        end for;
        continue;
    end if;
    A := rts0[1][1] - p0;
    bA, _ := IsSquare(A);
    printf "  A ne 0: %o ; A a square in Q(w) (kernel rational)? %o\n", A ne 0, bA;
    if A eq 0 then printf "  A = 0 (biquadratic degeneration) -- skipped\n"; continue; end if;

    L<aL> := ext<Fs | Sv^2 - A>;
    PL<XL> := PolynomialRing(L);
    sL := L!s0; pL := L!p0; qL := L!q0; AL := L!A;
    cc := ((pL + AL) + qL/aL)/2;
    bb := ((pL + AL) - qL/aL)/2;
    G1L := PL![L!Coefficient(B,i) : i in [0..2]];
    G2L := (XL + sL)^2 + aL*(XL + sL) + bb;
    G3L := (XL + sL)^2 - aL*(XL + sL) + cc;
    printf "  KERNEL CHECK G1G2G3 = F: %o\n", G1L*G2L*G3L eq PL![L!Coefficient(F,i) : i in [0..6]];
    M := Matrix(L, 3, 3, [Coefficient(G1L,2), Coefficient(G1L,1), Coefficient(G1L,0),
                          Coefficient(G2L,2), Coefficient(G2L,1), Coefficient(G2L,0),
                          Coefficient(G3L,2), Coefficient(G3L,1), Coefficient(G3L,0)]);
    Delta := Determinant(M);
    printf "  Delta ne 0: %o\n", Delta ne 0;
    if Delta eq 0 then continue; end if;
    H1 := Derivative(G2L)*G3L - G2L*Derivative(G3L);
    H2 := Derivative(G3L)*G1L - G3L*Derivative(G1L);
    H3 := Derivative(G1L)*G2L - G1L*Derivative(G2L);
    fpL := -H1*H2*H3/Delta;
    printf "  CODOMAIN descends to Q(w): %o\n", &and[Coefficient(fpL,i) in Fs : i in [0..6]];
    fp := PK![Fs!Coefficient(fpL,i) : i in [0..6]];
    facfp := Factorization(fp);
    printf "  CODOMAIN FACTOR DEGREES OVER Q(w): %o\n", [Degree(f[1]) : f in facfp];
    for f in facfp do printf "    CODFACTOR (%o): %o\n", Degree(f[1]), f[1]; end for;

    quarts := [f[1] : f in facfp | Degree(f[1]) eq 4];
    quads  := [f[1] : f in facfp | Degree(f[1]) eq 2];
    lins   := [f[1] : f in facfp | Degree(f[1]) eq 1];
    printf "  #quartics=%o #quadratics=%o #linears=%o\n", #quarts, #quads, #lins;
    if #quarts eq 0 then
        printf "  *** CODOMAIN SPLITS FURTHER over Q(w) -- generic 2-rank may be >= 2 ***\n";
        continue;
    end if;
    q4p := quarts[1]/LeadingCoefficient(quarts[1]);
    s1, p1, q1, r1, Res1 := ResData(q4p);
    printf "  IS q1 IDENTICALLY ZERO? %o\n", q1 eq 0;
    rts1 := Roots(Res1);
    printf "  image quartic resolvent roots in Q(w): %o\n", #rts1;
    if #rts1 gt 0 then
        A1 := rts1[1][1] - p1;
        b1s, _ := IsSquare(A1);
        printf "  IS A1 A SQUARE IN Q(w) (=> image ALWAYS 2-rank 2)? %o\n", b1s;
        if not b1s then
            LA := ext<Fs | Sv^2 - A1>;
            printf "  RAISE-LOCUS COVER {A1 square} : genus %o, constant field %o\n", Genus(LA), ConstantField(LA);
        end if;
        Dump("C43_B" cat IntegerToString(bi) cat "_A1", A1);
    end if;
    if #quads gt 0 then
        h1 := quads[1];
        d1 := Coefficient(h1,1)^2 - 4*Coefficient(h1,2)*Coefficient(h1,0);
        bd, _ := IsSquare(d1);
        printf "  IS disc(image quadratic) A SQUARE IN Q(w)? %o\n", bd;
        if not bd then
            LDq := ext<Fs | Sv^2 - d1>;
            printf "  COVER {image quadratic splits} : genus %o, constant field %o\n", Genus(LDq), ConstantField(LDq);
        end if;
        Dump("C43_B" cat IntegerToString(bi) cat "_D1", d1);
    end if;
    W := Discriminant(q4p);
    bw, _ := IsSquare(W);
    printf "  IS disc(q4') A SQUARE IN Q(w)? %o\n", bw;
    if not bw then
        LW := ext<Fs | Sv^2 - W>;
        printf "  COVER {disc(q4') square} : genus %o, constant field %o\n", Genus(LW), ConstantField(LW);
    end if;
    Dump("C43_B" cat IntegerToString(bi) cat "_W", W);
end for;

printf "C43PARAM_DONE\n";
quit;
