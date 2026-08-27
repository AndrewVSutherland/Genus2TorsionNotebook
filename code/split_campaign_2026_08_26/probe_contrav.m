SetColumns(0);
SetMemoryLimit(8*10^9);
AttachSpec("/home/claude/Magma/magma.spec");
Attach("/home/claude/torsion_jac/genus2.m");
Q := Rationals();
// E10(2)
tv := Q!2; de := tv^2-3*tv+1;
bv := tv^3*(tv-1)*(2*tv-1)/de^2; cv := -tv*(tv-1)*(2*tv-1)/de;
E1 := EllipticCurve([1-cv, -bv, -bv, 0, 0]);
printf "E1 tors %o\n", Invariants(TorsionSubgroup(E1));
U := GenusOneModel(3, E1);
printf "model U ok, degree %o\n", Degree(U);
H := Hessian(U);
printf "Hessian ok\n";
P, QQ := Contravariants(U);
printf "Contravariants ok\n";
for lam in [Q| 0, 1, -1, 2, 1/2, 5, -3 ] do
    // contravariant pencil member: P + lam*Q  (and mu=0 handled by lam='infty' below)
    ok := true; Em := 0;
    try
        M := P + lam*QQ;
        Em := Jacobian(M);
    catch e;
        ok := false;
    end try;
    if not ok then printf "lam=%o: degenerate\n", lam; continue; end if;
    // congruence test mod 3
    bad := 3*Conductor(E1)*Conductor(Em);
    mm := 0; nt := 0;
    for p in PrimesInInterval(5, 200) do
        if bad mod p eq 0 then continue; end if;
        nt +:= 1;
        if (TraceOfFrobenius(E1,p) - TraceOfFrobenius(Em,p)) mod 3 ne 0 then mm +:= 1; end if;
    end for;
    iso := IsIsomorphic(E1, Em);
    g3 := [];
    if mm eq 0 and not iso then
        try g3 := Genus2Elliptic3(E1, Em); catch e; end try;
    end if;
    printf "lam=%o: cond(Em)=%o iso=%o trace-mismatch=%o/%o glue3=%o\n",
        lam, Conductor(Em), iso, mm, nt, #g3;
end for;
// also the Hesse (covariant) pencil for comparison
for lam in [Q| 1, -1, 2 ] do
    ok := true; Em := 0;
    try
        M := U + lam*H;
        Em := Jacobian(M);
    catch e;
        ok := false;
    end try;
    if not ok then printf "hesse lam=%o: degenerate\n", lam; continue; end if;
    bad := 3*Conductor(E1)*Conductor(Em);
    mm := 0; nt := 0;
    for p in PrimesInInterval(5, 200) do
        if bad mod p eq 0 then continue; end if;
        nt +:= 1;
        if (TraceOfFrobenius(E1,p) - TraceOfFrobenius(Em,p)) mod 3 ne 0 then mm +:= 1; end if;
    end for;
    iso := IsIsomorphic(E1, Em);
    g3 := [];
    if mm eq 0 and not iso then
        try g3 := Genus2Elliptic3(E1, Em); catch e; end try;
    end if;
    printf "hesse lam=%o: cond(Em)=%o iso=%o trace-mismatch=%o/%o glue3=%o\n",
        lam, Conductor(Em), iso, mm, nt, #g3;
end for;
quit;
