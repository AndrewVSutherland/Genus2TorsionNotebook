// Which of {±SIL, ±SIL^t} satisfies M * PM = conj(PM)?  (the true
// complex-conjugation matrix on the Periods lattice frame)
SetColumns(0);
SetSeed(1);
SetMemoryLimit(8*10^9);
Prec := 40;
SetDefaultRealField(RealField(Prec + 10));
CC := ComplexField(Prec + 10);
Lv := 2190;
trtargets := [<7,2>, <11,-6>, <13,-2>, <17,6>, <19,10>];
M := ModularSymbols(Lv, 2, 0);
S := CuspidalSubspace(M);
NS := NewSubspace(S);
D := NewformDecomposition(NS);
target := 0;
for i in [1..#D] do
    if Dimension(D[i]) ne 4 then continue; end if;
    ok := true;
    for tt in trtargets do
        if Trace(HeckeOperator(D[i], tt[1])) ne 2*tt[2] then ok := false; break; end if;
    end for;
    if ok then target := i; break; end if;
end for;
Af := D[target];
Q := Rationals();
BAf := ChangeRing(BasisMatrix(VectorSpace(Af)), Q);
Lb := ChangeRing(BasisMatrix(Lattice(Af)), Q);
Cm := Solution(BAf, Lb);
SIm := ChangeRing(StarInvolution(Af), Q);
SIL := Cm * SIm * Cm^-1;
nterms := Ceiling(Lv * Prec * Log(10)/(2*Pi(RealField(20)))) + 200;
Pv := Periods(Af, nterms);
PMb := Matrix(CC, 4, 2, &cat[ [CC!Pv[i][1], CC!Pv[i][2]] : i in [1..4] ]);
PMc := Matrix(CC, 4, 2, [Conjugate(x) : x in Eltseq(PMb)]);
nrm := Max([Abs(x) : x in Eltseq(PMb)]);
T7 := ChangeRing(HeckeOperator(Af, 7), Q);
T7L := Cm * T7 * Cm^-1;
for tt in [<"SIL",SIL>, <"SILt",Transpose(SIL)>, <"SIm_basis",SIm>, <"SImt_basis",Transpose(SIm)>,
           <"T7_basis",T7>, <"T7t_basis",Transpose(T7)>, <"T7_lattice",T7L>, <"T7t_lattice",Transpose(T7L)>] do
    Mx := ChangeRing(ChangeRing(tt[2], Q), CC);
    X := Mx*PMb;                       // 4x2
    top := Submatrix(PMc, 1,1, 2,2);
    Rm := top^-1 * Submatrix(X, 1,1, 2,2);
    res := Max([Abs(x) : x in Eltseq(Submatrix(PMc,3,1,2,2)*Rm - Submatrix(X,3,1,2,2))]);
    printf "VARIANT %o residual %o R=%o\n", tt[1], RealField(6)!(res/nrm),
        Matrix(ComplexField(8),2,2,Eltseq(Rm));
end for;
printf "STAR_PERIOD_DONE\n";
quit;
