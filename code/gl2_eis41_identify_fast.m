// 41-hunt stage 2 (fast): identify the carrier on the FULL cuspidal space.
// AL slice (w2,w3,w5,w409)=(-,-,-,+) of S (includes old parts), restrict
// T7 (54s on S), factor charpoly of the restriction, report factors with
// root 8 mod 41 and their degrees.  Newness cross-checked afterwards.
SetColumns(0);
SetSeed(1);
SetMemoryLimit(120*10^9);
Level := 12270; P := 41;
Q := Rationals();
tt0 := Cputime();
M := ModularSymbols(Level, 2, 0);
S := CuspidalSubspace(M);
n := Dimension(S);
printf "DIM_S %o t=%o\n", n, Cputime(tt0);
W := VectorSpace(Q, n);
for tq in [<2,-1>, <3,-1>, <5,-1>, <409,1>] do
    tl := Cputime();
    ALq := ChangeRing(Matrix(AtkinLehnerOperator(S, tq[1])), Q);
    W := W meet Kernel(ALq - tq[2]*IdentityMatrix(Q, n));
    printf "AL%o done t=%o dim %o\n", tq[1], Cputime(tl), Dimension(W);
end for;
printf "AL_SLICE_DIM %o t=%o\n", Dimension(W), Cputime(tt0);
tl := Cputime();
T7 := ChangeRing(Matrix(HeckeOperator(S, 7)), Q);
printf "T7_FULL done t=%o\n", Cputime(tl);
BW := BasisMatrix(W);
T7W := Solution(BW, BW*T7);
cp7 := CharacteristicPolynomial(T7W);
printf "T7_SLICE_CHARPOLY computed t=%o\n", Cputime(tt0);
R41 := PolynomialRing(GF(P));
for fa in Factorization(cp7) do
    dn := LCM([Denominator(c) : c in Coefficients(fa[1])]);
    f41ok := dn mod P ne 0;
    hasroot := false;
    if f41ok then
        f41 := R41![GF(P)!(Integers()!(dn*c)) : c in Coefficients(fa[1])] * (GF(P)!dn)^-1;
        hasroot := Evaluate(f41, GF(P)!8) eq 0;
    end if;
    printf "FACTOR deg %o mult %o root8mod41 %o : %o\n", Degree(fa[1]), fa[2], hasroot, fa[1];
end for;
printf "EIS41_IDENTIFY_DONE t=%o\n", Cputime(tt0);
quit;
