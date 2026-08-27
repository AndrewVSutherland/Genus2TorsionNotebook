// Check: is StarInvolution(Af) integral in the Basis frame, and what is it
// in the lattice (= Periods) frame?  Decides the correct SId for v8.
SetColumns(0);
SetSeed(1);
SetMemoryLimit(8*10^9);
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
error if target eq 0, "piece not found";
Af := D[target];
Q := Rationals();
SIm := ChangeRing(StarInvolution(Af), Q);
BAf := ChangeRing(BasisMatrix(VectorSpace(Af)), Q);
Lb := ChangeRing(BasisMatrix(Lattice(Af)), Q);
Cm := Solution(BAf, Lb);
printf "SIM_BASIS_INTEGRAL %o\n", &and[Denominator(v) eq 1 : v in Eltseq(SIm)];
SIL := Cm * SIm * Cm^-1;
printf "SIL_LATTICE_INTEGRAL %o\n", &and[Denominator(v) eq 1 : v in Eltseq(SIL)];
printf "SIM %o\n", Eltseq(SIm);
printf "SIL %o\n", Eltseq(SIL);
printf "SAME %o\n", SIm eq SIL;
printf "STAR_FRAME_DONE\n";
quit;
