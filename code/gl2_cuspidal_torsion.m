// Cuspidal-torsion certification for GL2-type abelian surfaces A_f.
// Target (Filip, 2026-08-03): does A_f for f = 2190.2.a.v have a rational
// 37-torsion point?  (Alessandri-Coppola arXiv:2602.21047 predict order 37
// = their gcd bound; existence unproven.)  Also f = 1830.2.a.q (order 31).
//
// Method: level N squarefree => all cusps of X_0(N) rational => the image
// of the rational cuspidal subgroup in the optimal quotient A_f is rational
// torsion.  Magma: TorsionLowerBound(A) computes exactly this lower bound;
// TorsionMultiple(A) the upper multiple.  37 | lower bound => THEOREM.
//
// The newform piece is located inside the Atkin-Lehner eigenspace (cheap)
// and identified by Hecke traces against LMFDB:
//   2190.2.a.v: dim 2, K = Q(sqrt3), AL(2,3,5,73) = (-,-,-,+),
//               tr a_7, a_11, a_13, a_17, a_19 = 2, -6, -2, 6, 10
//   1830.2.a.q: dim 2, K = Q(sqrt2), AL(2,3,5,61) = (-,-,-,+),
//               tr a_7, a_11, a_13, a_17, a_19 = 0, 8, -2, 6, -6
// (sign-0 spaces double the trace).
//
// Run: magma -b Lv:=2190 code/gl2_cuspidal_torsion.m > results/gl2_cusp_2190.log

SetColumns(0);
SetSeed(1);
if not assigned Lv then Lv := 2190; elif Type(Lv) eq MonStgElt then Lv := StringToInteger(Lv); end if;
if not assigned MemGB then MemGB := 24; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

if Lv eq 2190 then
    ALsigns := [<2,-1>, <3,-1>, <5,-1>, <73,1>];
    trtargets := [<7,2>, <11,-6>, <13,-2>, <17,6>, <19,10>];
    targetprime := 37;
elif Lv eq 1830 then
    ALsigns := [<2,-1>, <3,-1>, <5,-1>, <61,1>];
    trtargets := [<7,0>, <11,8>, <13,-2>, <17,6>, <19,-6>];
    targetprime := 31;
else
    error "unknown level";
end if;

printf "LEVEL %o TARGETPRIME %o\n", Lv, targetprime;
M := ModularSymbols(Lv, 2, 0);
printf "AMBIENT_DIM %o\n", Dimension(M);
S := CuspidalSubspace(M);
NS := NewSubspace(S);
printf "NEW_CUSPIDAL_DIM %o\n", Dimension(NS);

D := NewformDecomposition(NS);
printf "PIECES %o dims=%o\n", #D, [Dimension(d) : d in D];

found := 0;
for i in [1..#D] do
    if Dimension(D[i]) ne 4 then continue; end if;
    ok := true;
    for tt in trtargets do
        tr := Trace(HeckeOperator(D[i], tt[1]));
        if tr ne 2*tt[2] then ok := false; break; end if;
    end for;
    if not ok then continue; end if;
    found +:= 1;
    als := [ Integers()!(Trace(AtkinLehnerOperator(D[i], al[1]))/4) : al in ALsigns ];
    printf "MATCHED_PIECE %o (all 5 Hecke traces agree) AL_signs=%o (expect %o)\n",
           i, als, [al[2] : al in ALsigns];
    A := ModularAbelianVariety(D[i]);
    printf "MODABVAR dim=%o\n", Dimension(A);
    tl := TorsionLowerBound(A);
    printf "TORSION_LOWER_BOUND %o\n", tl;
    tm := TorsionMultiple(A);
    printf "TORSION_MULTIPLE %o\n", tm;
    if tl mod targetprime eq 0 then
        printf "THEOREM: A_f has a rational point of order %o (cuspidal)\n", targetprime;
    else
        printf "CUSPIDAL_BOUND_INSUFFICIENT (lower=%o; %o-part unproven this way)\n", tl, targetprime;
    end if;
end for;
printf "MATCHED_TOTAL %o\n", found;
printf "GL2_CUSPIDAL_DONE level=%o\n", Lv;
quit;
