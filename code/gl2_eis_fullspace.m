// Eisenstein-kernel prover WITHOUT NewformDecomposition: intersect
// ker(T_l - (1+l)) mod P over the FULL cuspidal space (integral lattice),
// plus U_q targets at q | N.  Validation: Level 2190, P 37 must give dim >= 2.
// Production: Level 12270 = 2*3*5*409, P 41 (409 = 10*41 - 1).
//
// magma -b Level:=2190 P:=37 code/gl2_eis_fullspace.m
// magma -b Level:=12270 P:=41 code/gl2_eis_fullspace.m

SetColumns(0);
SetSeed(1);
if not assigned Level then Level := 2190; elif Type(Level) eq MonStgElt then Level := StringToInteger(Level); end if;
if not assigned P then P := 37; elif Type(P) eq MonStgElt then P := StringToInteger(P); end if;
if not assigned MemGB then MemGB := 80; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

printf "LEVEL %o P %o\n", Level, P;
tt0 := Cputime();
M := ModularSymbols(Level, 2, 0);
S := CuspidalSubspace(M);
printf "DIM_CUSPIDAL %o (t=%o)\n", Dimension(S), Cputime(tt0);
Q := Rationals();
BAf := ChangeRing(BasisMatrix(VectorSpace(S)), Q);
Lb := ChangeRing(BasisMatrix(Lattice(S)), Q);
Cm := Solution(BAf, Lb);
F := GF(P);
n := Dimension(S);
W := VectorSpace(F, n);
// good primes first (cheap, big cuts), then bad U_q
goodells := [ l : l in [7, 11, 13, 17, 19] | Level mod l ne 0 ];
badq := [ f[1] : f in Factorization(Level) ];
// bad targets: U_l = 1 for the small primes, U_q = q for the big one (w_q = +1 pattern)
for l in goodells do
    tl := Cputime();
    Tl := HeckeOperator(S, l);
    TlL := Cm * ChangeRing(Tl, Q) * Cm^-1;
    dn := LCM([Denominator(x) : x in Eltseq(TlL)]);
    error if dn mod P eq 0, "denominator hits P";
    TlF := ChangeRing(Matrix(Q,n,n,[dn*x : x in Eltseq(TlL)]), F) * (F!dn)^-1;
    W := W meet Kernel(TlF - ScalarMatrix(F, n, F!(1+l)));
    printf "T%o done t=%o dimW %o\n", l, Cputime(tl), Dimension(W);
    if Dimension(W) eq 0 then break; end if;
end for;
if Dimension(W) gt 0 then
    for q in badq do
        // try both classical Eisenstein U_q eigenvalues: 1 and q
        tl := Cputime();
        Tq := HeckeOperator(S, q);
        TqL := Cm * ChangeRing(Tq, Q) * Cm^-1;
        dn := LCM([Denominator(x) : x in Eltseq(TqL)]);
        error if dn mod P eq 0, "denominator hits P";
        TqF := ChangeRing(Matrix(Q,n,n,[dn*x : x in Eltseq(TqL)]), F) * (F!dn)^-1;
        W1 := W meet Kernel(TqF - ScalarMatrix(F, n, F!1));
        Wq := W meet Kernel(TqF - ScalarMatrix(F, n, F!q));
        printf "U%o done t=%o dims: eig1 %o eigq %o\n", q, Cputime(tl), Dimension(W1), Dimension(Wq);
        // keep the union info; continue with the larger (both tracked would fork; keep bigger)
        W := Dimension(W1) ge Dimension(Wq) select W1 else Wq;
        printf "U%o keep dim %o\n", q, Dimension(W);
        if Dimension(W) eq 0 then break; end if;
    end for;
end if;
printf "EIS_KERNEL_FINAL_DIM %o\n", Dimension(W);
if Dimension(W) gt 0 then
    // identify the Hecke field of the piece: act with T_l for one good l on W-lift?
    // (mod-P kernel only; full identification needs char-0 work -- report basis)
    for b in Basis(W) do printf "KERNEL_VEC %o\n", b; end for;
end if;
printf "EIS_FULLSPACE_DONE t=%o\n", Cputime(tt0);
quit;
