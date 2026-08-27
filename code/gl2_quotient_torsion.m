// Cuspidal torsion in the OPTIMAL QUOTIENT A_f = J_0(N)/B for
// f = 2190.2.a.v (target 37) and f = 1830.2.a.q (target 31).
// The subvariety-side TorsionLowerBound gave 1 (results/gl2_cusp_*.log):
// the cuspidal group meets the SUBvariety trivially, as expected; the
// canonical home of cuspidal torsion is the quotient.  Here: build
// J = JZero(N) in the ModAbVar category, identify the factor attached to
// the newform (Hecke traces), form B = sum of the OTHER factors, and
// compute TorsionLowerBound / TorsionMultiple of Q = J/B.
//
// Run: magma -b Lv:=2190 code/gl2_quotient_torsion.m > results/gl2_quot_2190.log

SetColumns(0);
SetSeed(1);
if not assigned Lv then Lv := 2190; elif Type(Lv) eq MonStgElt then Lv := StringToInteger(Lv); end if;
if not assigned MemGB then MemGB := 24; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

if Lv eq 2190 then
    trtargets := [<7,2>, <11,-6>, <13,-2>, <17,6>, <19,10>];
    targetprime := 37;
elif Lv eq 1830 then
    trtargets := [<7,0>, <11,8>, <13,-2>, <17,6>, <19,-6>];
    targetprime := 31;
else
    error "unknown level";
end if;

printf "LEVEL %o TARGET %o\n", Lv, targetprime;
J := JZero(Lv);
D := Decomposition(J);
printf "FACTORS %o dims=%o\n", #D, [Dimension(d) : d in D];

// identify the factor: dimension 2 (ModAbVar dims are true dims), new of
// level Lv, matching Hecke traces on its modular symbols
target := 0;
for i in [1..#D] do
    if Dimension(D[i]) ne 2 then continue; end if;
    ms := ModularSymbols(D[i])[1];
    if Level(ms) ne Lv then continue; end if;
    ok := true;
    for tt in trtargets do
        if Trace(HeckeOperator(ms, tt[1])) ne 2*tt[2] then ok := false; break; end if;
    end for;
    if ok then target := i; break; end if;
end for;
error if target eq 0, "newform factor not found";
printf "TARGET_FACTOR %o\n", target;

B := &+[ D[i] : i in [1..#D] | i ne target ];
printf "COMPLEMENT_DIM %o\n", Dimension(B);
Q := J/B;
printf "QUOTIENT_DIM %o\n", Dimension(Q);
tl := TorsionLowerBound(Q);
printf "Q_TORSION_LOWER_BOUND %o\n", tl;
tm := TorsionMultiple(Q);
printf "Q_TORSION_MULTIPLE %o\n", tm;
if tl mod targetprime eq 0 then
    printf "THEOREM: the optimal quotient A_f has a rational point of order %o\n", targetprime;
else
    printf "QUOTIENT_CUSPIDAL_INSUFFICIENT lower=%o\n", tl;
end if;
printf "GL2_QUOT_DONE level=%o\n", Lv;
quit;
