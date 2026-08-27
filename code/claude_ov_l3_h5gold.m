// ===========================================================================
// Lane 3 : GOLD-STANDARD validation of the derived Humbert-5 equation, and of
// the repaired RM screen, against LMFDB's certified endomorphism data.
//
// Ground truth: the 116 curves of the LMFDB production genus-2 database with
// geom_end_alg = 'RM' that carry a g2c_endomorphisms row.  For each we have
//   * factorsQQ_geom[0][0] = the LMFDB number-field label of the RM FIELD
//     ('2.2.5.1' = Q(sqrt 5), '2.2.8.1' = Q(sqrt 2), '2.2.12.1' = Q(sqrt 3),
//      '2.2.13.1', '2.2.17.1');
//   * ring_geom[0] = the INDEX of the RM order in the maximal order.
// H_5 is the Humbert surface of DISCRIMINANT 5, i.e. RM by the MAXIMAL order
// of Q(sqrt 5).  So the two exact predictions are
//   H5 = 0    <=>   field = 2.2.5.1  AND  index = 1
//   repaired census core  =  the squarefree core of the field discriminant
//                            (5, 2, 3, 13, 17 for the labels above),
//                            for EVERY curve including index > 1.
//
// usage: code/claude_magma_slot.sh -b code/claude_ov_l3_h5gold.m
// ===========================================================================
SetColumns(0);
SetMemoryLimit(4*10^9);
if not assigned IN then IN := "/tmp/claude-1000/-home-claude-torsion-jac/1c22f408-33bb-486e-870f-c5fca84bb522/scratchpad/l3_gold.txt"; end if;
if not assigned PB then PB := 400; elif Type(PB) eq MonStgElt then PB := StringToInteger(PB); end if;
if not assigned OUT then OUT := "results/claude_ov_l3_h5gold.log"; end if;

Z := Integers();  Q := Rationals();
Px<x> := PolynomialRing(Q);
R<J2,J4,J6,J10> := PolynomialRing(Q, 4);
H5 := R ! eval Read("results/claude_ov_l3_humbert5_eqn.txt");
printf "H5 loaded: #terms=%o\n", #Terms(H5);

// squarefree core of the discriminant of the field with LMFDB label 2.2.D.1
function CoreOfLabel(lab)
    s := Split(lab, ".");
    return SquarefreeFactorization(StringToInteger(s[3]));
end function;

function Census(g, B)
    dsc := Z ! Discriminant(g);
    dn := {};  dr := {};  np := 0; na1z := 0; nsq := 0;
    for p in PrimesInInterval(3, B) do
        if dsc mod p eq 0 then continue; end if;
        fp := PolynomialRing(GF(p)) ! g;
        if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
        L := LPolynomial(HyperellipticCurve(fp));
        a1 := Z!(-Coefficient(L,1));  a2 := Z!Coefficient(L,2);
        n := a1^2 - 4*(a2 - 2*p);
        np +:= 1;
        if n eq 0 then continue; end if;
        d := SquarefreeFactorization(n);
        if d eq 1 then nsq +:= 1; continue; end if;
        Include(~dn, d);
        if a1 eq 0 then na1z +:= 1; continue; end if;
        Include(~dr, d);
    end for;
    return dn, dr, np, na1z, nsq;
end function;

lines := [l : l in Split(Read(IN), "\n") | #l gt 5];
printf "GOLD_CURVES %o  PB=%o\n", #lines, PB;
G := Open(OUT, "w");
Puts(G, Sprintf("GOLD_CURVES %o PB=%o", #lines, PB));

nh5ok := 0; nh5bad := 0; nnaive := 0; nrep := 0; ntot := 0;
for l in lines do
    s := Split(l, ";");
    // fields: label ; f ; h ; nflabel ; index ; gl2      (h may be empty ->
    // Magma's Split drops it, so realign from the end)
    lab := s[1];
    nfl := s[#s-2];  idx := StringToInteger(s[#s-1]);  gl2 := s[#s];
    fc := [StringToInteger(t) : t in Split(s[2], ",")];
    hc := [Z|];
    if #s ge 6 then hc := [StringToInteger(t) : t in Split(s[3], ",")]; end if;
    f := Px ! fc;  hh := Px ! hc;
    g := 4*f + hh^2;
    if Degree(g) lt 5 or Discriminant(g) eq 0 then
        Puts(G, "ROW " cat lab cat " BAD"); continue;
    end if;
    C := HyperellipticCurve(g);
    I := IgusaClebschInvariants(C);
    h5z := Evaluate(H5, [Q!I[1],Q!I[2],Q!I[3],Q!I[4]]) eq 0;
    dn, dr, np, na1z, nsq := Census(g, PB);
    core := CoreOfLabel(nfl);
    pred := (nfl eq "2.2.5.1") and (idx eq 1);
    okH5 := (h5z eq pred);
    okRep := (dr eq {core});
    okNaive := (dn eq {core});
    ntot +:= 1;
    if okH5 then nh5ok +:= 1; else nh5bad +:= 1; end if;
    if okNaive then nnaive +:= 1; end if;
    if okRep then nrep +:= 1; end if;
    line := Sprintf("ROW %o nf=%o idx=%o gl2=%o | H5zero=%o pred=%o %o | naive=%o %o | rep=%o %o | np=%o na1z=%o nsq=%o",
        lab, nfl, idx, gl2, h5z, pred, okH5 select "OK" else "MISMATCH",
        Sort(Setseq(dn)), okNaive select "OK" else "SCATTER",
        Sort(Setseq(dr)), okRep select "OK" else "SCATTER", np, na1z, nsq);
    printf "%o\n", line;  Puts(G, line);
end for;
res := Sprintf("GOLD_RESULT curves=%o  H5 correct=%o mismatch=%o | naive census correct=%o (%o%%) | repaired census correct=%o (%o%%)",
    ntot, nh5ok, nh5bad, nnaive, RealField(4)!(100*nnaive/ntot), nrep, RealField(4)!(100*nrep/ntot));
printf "%o\n", res;  Puts(G, res);
Flush(G); delete G;
printf "GOLD_DONE\n";
quit;
