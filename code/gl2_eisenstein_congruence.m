// Eisenstein congruence certification for f = 2190.2.a.v mod 37
// (and 1830.2.a.q mod 31): proves reducibility of the mod-p_frak Galois
// representation, hence (standard isogeny-shift argument, cf. Katz 1981)
// the existence of a rational point of order 37 (resp. 31) on SOME abelian
// surface isogenous to A_f -- exactly the form of the AC prediction.
//
// Method: on the f-isotypic modular symbols piece Af (dim 4), with its
// integral lattice L, compute W := the intersection over ell <= Sturm
// bound (and p | N) of ker(T_ell - (1+ell)  mod P) resp.
// ker(T_p - ap_target mod P) on L/P L, P = 37 or 31.
// W != 0  ==>  the Eisenstein maximal ideal is in the support of the
// f-part of H_1 mod P  ==>  a_ell(f) == 1 + ell mod some prime above P
// for all ell (Deligne-Serre lifting / Brauer-Nesbitt)  ==>  rhobar_f
// reducible with semisimplification 1 + chi_cyc.
// Eisenstein a_p targets at p | N (matching f's AL signs):
//   2190: (a2,a3,a5,a73) = (1,1,1,73)  [73 == -1 mod 37 matches w73=+1]
//   1830: (a2,a3,a5,a61) = (1,1,1,61)  [61 == -1 mod 31 matches w61=+1]
// Sturm bound: k/12 * [SL2(Z):Gamma0(N)] = index/6.
//
// Run: magma -b Lv:=2190 code/gl2_eisenstein_congruence.m > results/gl2_eis_2190.log

SetColumns(0);
SetSeed(1);
if not assigned Lv then Lv := 2190; elif Type(Lv) eq MonStgElt then Lv := StringToInteger(Lv); end if;
if not assigned MemGB then MemGB := 16; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

if Lv eq 2190 then
    trtargets := [<7,2>, <11,-6>, <13,-2>, <17,6>, <19,10>];
    P := 37;
    badtargets := [<2,1>, <3,1>, <5,1>, <73,73>];
elif Lv eq 1830 then
    trtargets := [<7,0>, <11,8>, <13,-2>, <17,6>, <19,-6>];
    P := 31;
    badtargets := [<2,1>, <3,1>, <5,1>, <61,61>];
elif Lv eq 23 then
    trtargets := [];   // validation: single piece, Eisenstein prime 11
    P := 11;
    badtargets := [<23,1>];
else
    error "unknown level";
end if;

idx := Index(Gamma0(Lv));
sturm := Ceiling(idx/6);
printf "LEVEL %o P %o STURM %o\n", Lv, P, sturm;

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
printf "PIECE %o\n", target;

// integral lattice of the Af piece: intersect integral cuspidal homology
// with the Af subspace; work in the Af basis
BAf := BasisMatrix(VectorSpace(Af));
LatS := Lattice(Af);
printf "LATTICE_RANK %o\n", Rank(LatS);
Lb := BasisMatrix(LatS);   // rows: lattice basis in ambient coords
// Hecke matrices on the lattice basis, reduced mod P
F := GF(P);
W := VectorSpace(F, 4);
dimW := 4;
ells := [ l : l in PrimesUpTo(sturm) | Lv mod l ne 0 ];
count := 0;
for l in ells cat [bt[1] : bt in badtargets] do
    if Lv mod l eq 0 then
        av := [ bt[2] : bt in badtargets | bt[1] eq l ][1];
    else
        av := 1 + l;
    end if;
    Tl := HeckeOperator(Af, l);       // matrix on Basis(Af) (rational)
    // express on lattice basis: Tl_L := Lb * Tl_amb * Lb^-1 -- but
    // HeckeOperator(Af, l) is ALREADY w.r.t. Basis(Af); move to lattice:
    // rows of Lb are in ambient coords; convert: C := coordinates of
    // lattice basis in Basis(Af): C := Lb * pseudoinv(BAf)
    C := Solution(ChangeRing(BAf, Rationals()), ChangeRing(Lb, Rationals()));  // C * BAf = Lb
    TlL := C * Tl * C^-1;
    den := LCM([Denominator(x) : x in Eltseq(TlL)]);
    error if den mod P eq 0, "denominator hits P -- lattice not Hecke-stable?";
    TlF := Matrix(F, 4, 4, [ F!(Numerator(x)*InverseMod(Denominator(x), P)) : x in Eltseq(TlL) ]);
    W := W meet Kernel(TlF - ScalarMatrix(F, 4, F!av));
    count +:= 1;
    if Dimension(W) eq 0 then
        printf "KERNEL_EMPTY after %o operators (last l=%o)\n", count, l;
        break;
    end if;
    if count mod 25 eq 0 then printf "PROGRESS ops=%o dimW=%o\n", count, Dimension(W); end if;
end for;
printf "FINAL_EISENSTEIN_KERNEL_DIM %o (after %o Hecke operators to Sturm=%o)\n",
       Dimension(W), count, sturm;
if Dimension(W) gt 0 then
    printf "THEOREM: a_ell(f) == 1+ell mod a prime above %o for all ell (Eisenstein congruence proven to Sturm bound); rhobar reducible; some abelian surface isogenous to A_f has a rational point of order %o\n", P, P;
else
    printf "NO_EISENSTEIN_KERNEL: congruence fails -- AC prediction refuted for this form?!\n";
end if;
printf "GL2_EIS_DONE level=%o\n", Lv;
quit;
