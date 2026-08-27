// 41-hunt stage 2: identify the newform carrying the mod-41 Eisenstein
// system at level 12270.  AL signs forced by the kernel's U-eigenvalues:
// w2 = w3 = w5 = -1 (U_q = 1), w409 = +1 (U_409 = 409 = -1 mod 41).
// Slice the newspace by those AL signs, decompose only the slice, and
// find the piece(s) whose T_7 charpoly mod 41 has root 8 (= 1+7).
//
// magma -b code/gl2_eis41_identify.m > results/gl2_eis41_identify.log

SetColumns(0);
SetSeed(1);
SetMemoryLimit(120*10^9);
Level := 12270;
P := 41;
tt0 := Cputime();
M := ModularSymbols(Level, 2, 0);
S := CuspidalSubspace(M);
printf "DIM_S %o t=%o\n", Dimension(S), Cputime(tt0);
NS := NewSubspace(S);
printf "DIM_NEW %o t=%o\n", Dimension(NS), Cputime(tt0);
Q := Rationals();
dN := Dimension(NS);
W := VectorSpace(Q, dN);
for tq in [<2,-1>, <3,-1>, <5,-1>, <409,1>] do
    tl := Cputime();
    ALq := ChangeRing(Matrix(AtkinLehnerOperator(NS, tq[1])), Q);
    W := W meet Kernel(ALq - tq[2]*IdentityMatrix(Q, dN));
    printf "AL%o done t=%o dim %o\n", tq[1], Cputime(tl), Dimension(W);
end for;
printf "AL_SLICE_DIM %o t=%o\n", Dimension(W), Cputime(tt0);
// W is a vector space; find newform pieces inside: decompose NS but skip
// pieces not in the slice: cheaper: charpoly of T7 restricted to W
BW := BasisMatrix(W);
T7 := ChangeRing(Matrix(HeckeOperator(NS, 7)), Q);
T7W := Solution(BW, BW*T7);
cp7 := CharacteristicPolynomial(T7W);
printf "T7_SLICE_CHARPOLY_FACTORS:\n";
for fa in Factorization(cp7) do
    printf "  deg %o mult %o : %o\n", Degree(fa[1]), fa[2], fa[1];
    // check root 8 mod 41
    R41 := PolynomialRing(GF(P));
    f41 := R41![GF(P)!(Integers()!(Numerator(c)) * Modinv(Integers()!Denominator(c), P)) : c in Coefficients(fa[1])];
    if Evaluate(f41, GF(P)!8) eq 0 then
        printf "    ^^ HAS ROOT a_7 = 8 mod 41 (Eisenstein candidate, Hecke degree %o)\n", Degree(fa[1]);
    end if;
end for;
printf "EIS41_IDENTIFY_DONE t=%o\n", Cputime(tt0);
quit;
