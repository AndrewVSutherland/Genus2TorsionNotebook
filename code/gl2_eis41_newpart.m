// 41-hunt stage 3: the NEW part of the AL slice, exactly.
// NS is T7- and AL-stable; restrict the (cheap) full-space T7 to
// NS meet slice and factor the 42x42 charpoly.
SetColumns(0);
SetSeed(1);
SetMemoryLimit(120*10^9);
if not assigned Level then Level := 12270; elif Type(Level) eq MonStgElt then Level := StringToInteger(Level); end if;
P := 41;
Q := Rationals();
tt0 := Cputime();
M := ModularSymbols(Level, 2, 0);
S := CuspidalSubspace(M);
n := Dimension(S);
W := VectorSpace(Q, n);
qbig := [f[1] : f in Factorization(Level) | f[1] gt 5][1];
if not assigned WBig then WBig := 1; elif Type(WBig) eq MonStgElt then WBig := StringToInteger(WBig); end if;
altq := [ <f[1], f[1] eq qbig select WBig else -1> : f in Factorization(Level) ];
for tq in altq do
    ALq := ChangeRing(Matrix(AtkinLehnerOperator(S, tq[1])), Q);
    W := W meet Kernel(ALq - tq[2]*IdentityMatrix(Q, n));
end for;
printf "SLICE %o t=%o\n", Dimension(W), Cputime(tt0);
NS := NewSubspace(S);
printf "DIM_NEW %o t=%o\n", Dimension(NS), Cputime(tt0);
// NS as subspace of the S coordinate space
BNS := ChangeRing(BasisMatrix(VectorSpace(NS)), Q);
BS := ChangeRing(BasisMatrix(VectorSpace(S)), Q);
CNS := Solution(BS, BNS);   // rows: NS basis in S-coordinates
VNS := VectorSpaceWithBasis(CNS);
WN := W meet VNS;
printf "NEW_SLICE_DIM %o t=%o\n", Dimension(WN), Cputime(tt0);
T7 := ChangeRing(Matrix(HeckeOperator(S, 7)), Q);
BWN := BasisMatrix(WN);
T7W := Solution(BWN, BWN*T7);
cp7 := CharacteristicPolynomial(T7W);
R41 := PolynomialRing(GF(P));
for fa in Factorization(cp7) do
    dn := LCM([Denominator(c) : c in Coefficients(fa[1])]);
    hasroot := false;
    if dn mod P ne 0 then
        f41 := R41![GF(P)!(Integers()!(dn*c)) : c in Coefficients(fa[1])] * (GF(P)!dn)^-1;
        hasroot := Evaluate(f41, GF(P)!8) eq 0;
    end if;
    printf "NEWFACTOR deg %o mult %o root8mod41 %o : %o\n", Degree(fa[1]), fa[2], hasroot, fa[1];
end for;
printf "EIS41_NEWPART_DONE t=%o\n", Cputime(tt0);
quit;
