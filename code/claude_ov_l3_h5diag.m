// ===========================================================================
// Lane 3 : DIAGNOSE the H5-vs-census disagreement found by claude_ov_l3_h5validate.
//
// h5validate cross-tabulated the derived Humbert-5 equation H5 against the
// Frobenius real-subfield-disc census on 1000 LMFDB geom_end_alg='RM' curves:
//     H5=0 & census={5}    103
//     H5=0 & census scatters 242      <-- must be explained
//     H5<>0 & census={5}    30        <-- must be explained
//     H5<>0 & census={d}, d<>5  254
//     H5<>0 & census scatters 371     <-- must be explained (input is ALL RM!)
// Since every input curve has geometric RM, a SCATTERING census is by itself a
// falsification of the census screen, independent of H5.  So diagnose the
// census first.
//
// (A) per-prime detail on a few scattering rows: a1,a2,n,core, L-poly
//     factorisation type over Q, p-rank (ordinary / p-rank 1 / supersingular).
// (B) negative control: random sextics must have H5 <> 0.
// (C) forward control: random (g,h) on the Elkies-Kumar parametrisation ->
//     Mestre curve -> census must be {5}.  (This is step (c) of
//     claude_ov_l3_humbert5.m, which crashed there on a FldQuadElt coercion.)
//
// usage: code/claude_magma_slot.sh -b code/claude_ov_l3_h5diag.m
// ===========================================================================
SetColumns(0);
SetMemoryLimit(4*10^9);
Z := Integers();  Q := Rationals();
Px<x> := PolynomialRing(Q);
PT<T> := PolynomialRing(Q);
R<J2,J4,J6,J10> := PolynomialRing(Q, 4);
H5 := R ! eval Read("results/claude_ov_l3_humbert5_eqn.txt");
printf "H5 loaded: #terms=%o\n", #Terms(H5);

function H5Val(f)
    C := HyperellipticCurve(f);
    I := IgusaClebschInvariants(C);
    return Evaluate(H5, [Q!I[1], Q!I[2], Q!I[3], Q!I[4]]);
end function;

// ---------------------------------------------------------------- (A) detail
// rows taken from results/claude_ov_l3_h5validate.log; eqn from LMFDB
// (f-coeffs ; h-coeffs), model y^2 + h y = f, so g = 4f + h^2.
rows := [
 <"21168.b3",  [-18,-27,12,60,76,44,17], [Z|]>,      // H5=0, census scattered (10)
 <"112500.a2", [Z|], [Z|]>,                          // placeholder, filled below
 <"21168.b8",  [3,30,-36,-654,-240,288,45], [Z|]>,   // H5=0, scattered
 <"960400.hk2",[-72,196,-175,-35,175,-84,-7], [Z|1]>,
 <"1012500.ee2",[-24,-145,0,337,0,-135,-34],[Z|1,0,0,1]>,
 <"162000.ff6",[0,15,-165,390,210,15], [Z|]>,
 <"155236.a2", [128,32,-40,-9,14,4,-2], [Z|0,1,1]>,
 <"824464.a1", [2,-1,4,3,-3,2,1], [Z|]>,
 <"152881.b1", [-2,4,-3,-2,2], [Z|1,1,0,1]>
];
rows := [r : r in rows | #r[2] gt 0];

printf "\n==== (A) per-prime census detail ====\n";
for r in rows do
    f := Px ! r[2];  hh := Px ! r[3];
    g := 4*f + hh^2;
    v := H5Val(g);
    dsc := Z ! Discriminant(g);
    printf "\nCURVE %o  H5zero=%o  deg=%o\n", r[1], v eq 0, Degree(g);
    discs := {};  nord := 0; nnonord := 0; nsplit := 0; nbad := 0;
    detail := [];
    for p in PrimesInInterval(3, 200) do
        if dsc mod p eq 0 then nbad +:= 1; continue; end if;
        Fp := GF(p);
        fp := PolynomialRing(Fp) ! g;
        if Degree(fp) lt 5 or not IsSquarefree(fp) then nbad +:= 1; continue; end if;
        Cp := HyperellipticCurve(fp);
        L := LPolynomial(Cp);
        a1 := Z ! (-Coefficient(L,1));  a2 := Z ! Coefficient(L,2);
        n := a1^2 - 4*(a2 - 2*p);
        d0 := (n eq 0) select 0 else SquarefreeFactorization(n);
        // p-rank from the Newton polygon of L: ordinary <=> p does not divide a2... use
        // the standard test: ordinary iff a2 is coprime to p and a1 coprime to p
        prk := 2;
        if a1 mod p eq 0 then prk := (a2 mod p eq 0) select 0 else 1; end if;
        LQ := PT ! L;
        fct := [<Degree(t[1]), t[2]> : t in Factorization(LQ)];
        Include(~discs, d0);
        Append(~detail, <p, a1, a2, n, d0, prk, fct>);
    end for;
    printf "  discs = %o\n", Sort(Setseq(discs));
    for t in detail do
        if t[5] ne 5 then
            printf "  p=%3o a1=%o a2=%o n=%o core=%o prank=%o Lfact=%o\n",
                   t[1], t[2], t[3], t[4], t[5], t[6], t[7];
        end if;
    end for;
    n5 := #[t : t in detail | t[5] eq 5];
    printf "  SUMMARY %o : primes=%o core5=%o other=%o  (nonordinary among 'other': %o)\n",
        r[1], #detail, n5, #detail - n5, #[t : t in detail | t[5] ne 5 and t[6] lt 2];
end for;

// ------------------------------------------------------- (B) negative control
printf "\n==== (B) negative control: random sextics ====\n";
nz := 0; nt := 0;
for trial in [1..300] do
    f := x^6 + &+[Random([-20..20])*x^i : i in [0..5]];
    if Discriminant(f) eq 0 then continue; end if;
    nt +:= 1;
    if H5Val(f) ne 0 then nz +:= 1; end if;
end for;
printf "RANDOM_SEXTICS tested=%o  H5nonzero=%o\n", nt, nz;

// -------------------------------------------------------- (C) forward control
printf "\n==== (C) forward control: random (g,h) on the EK parametrisation ====\n";
function DiscCensusRat(f, B)
    m := LCM([Denominator(c) : c in Coefficients(f)] cat [Z|1]);
    fI := Px ! [ Q ! (Coefficient(f,i)*m^(6-i)) : i in [0..Degree(f)] ];
    fI := Px ! [ Z ! c : c in Coefficients(fI) ];
    dsc := Z ! Numerator(Discriminant(fI));
    discs := {};  np := 0;
    for p in PrimesInInterval(3, B) do
        if dsc mod p eq 0 then continue; end if;
        fp := PolynomialRing(GF(p)) ! fI;
        if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
        L := LPolynomial(HyperellipticCurve(fp));
        a1 := Z!(-Coefficient(L,1)); a2 := Z!Coefficient(L,2);
        n := a1^2 - 4*(a2 - 2*p);
        if n ne 0 then Include(~discs, SquarefreeFactorization(n)); end if;
        np +:= 1;
    end for;
    return discs, np;
end function;

cnt := 0;
for trial in [1..200] do
    gv := Random([-8..8]) + Random([0..6])/Random([1..6]);
    hv := Random([-8..8]) + Random([0..6])/Random([1..6]);
    IC := [Q| 24*gv+6, 9*gv^2, 81*gv^3+18*gv^2+36*hv, 4*hv^2 ];
    if IC[4] eq 0 or IC[2] eq 0 then continue; end if;
    ok := true;  f := Px ! 0;
    try
        C := HyperellipticCurveFromIgusaClebsch(IC);
        if BaseRing(C) cmpne Q then ok := false; else
            f := Px ! HyperellipticPolynomials(SimplifiedModel(C));
        end if;
    catch e
        ok := false;
    end try;
    if not ok or Degree(f) lt 5 then continue; end if;
    hv5 := H5Val(f);
    dc, np := DiscCensusRat(f, 120);
    cnt +:= 1;
    printf "  FWD (g,h)=(%o,%o) H5zero=%o nprimes=%o census=%o\n",
        gv, hv, hv5 eq 0, np, Sort(Setseq(dc));
    if cnt ge 12 then break; end if;
end for;

printf "\nDIAG_DONE\n";
quit;
