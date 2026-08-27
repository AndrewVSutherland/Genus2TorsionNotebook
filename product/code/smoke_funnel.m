// smoke_funnel.m — validate glue + funnel end-to-end on 210.e6 x 30.a6
// ([2,8] x [2,6], both full rational 2-torsion; lane8's [2,24] anchor came
// from this conductor pair), with exact stage forced on one output.
SetColumns(0);
SetMemoryLimit(3*10^9);
load "split_lab.m";  // run from product/code/

E1 := EllipticCurve([1,0,0,-1070,7812]);   // 210.e6  [2,8]
E2 := EllipticCurve([1,0,1,-19,26]);       // 30.a6   [2,6]
printf "E1 %o  E2 %o\n", Invariants(TorsionSubgroup(E1)), Invariants(TorsionSubgroup(E2));
L := Genus2Elliptic2(E1,E2);
printf "glued: %o curves\n", #L;
for k in [1..#L] do
    C := L[k];
    ok, g := IntegralSextic(C);
    printf "-- curve %o: g = %o\n", k, g;
    invsList, ps := ProfileCurve(g, 8);
    printf "   primes %o\n", ps;
    for i in [1..#invsList] do printf "   p=%o invs=%o\n", ps[i], invsList[i]; end for;
    compat := CompatibleGroups(invsList);
    printf "   compat = %o\n", compat;
    unknown := [ v : v in compat | not v in KNOWN ];
    printf "   unknown = %o\n", unknown;
    t0 := Cputime();
    einv := ExactTorsion(g);
    printf "   EXACT = %o  (%o s)\n", einv, Cputime()-t0;
end for;
quit;
