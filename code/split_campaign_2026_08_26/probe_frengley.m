// probe: Z(3,2) point -> (E1,E2) via ZNrModuli; verify anti-3-congruence + glue
SetColumns(0);
SetMemoryLimit(4*10^9);
AttachSpec("/home/claude/Magma/magma.spec");
Attach("/home/claude/torsion_jac/product/N-congruences/ZNr-equations.m");
Attach("/home/claude/torsion_jac/genus2.m");

Q := Rationals();
c1f := func<v | -432*v^3 - 216*v^2 + 576*v - 184>;
c0f := func<v | 4*v^6 - 20*v^5 + 41*v^4 - 44*v^3 + 26*v^2 - 8*v + 1>;

pts := [];
for vv in [Q| 2, 3, 1/2, -1, 5] do
    B := c0f(vv) - c1f(vv)^2/46656;
    if B eq 0 then continue; end if;
    for d in [Q| 1, 2, 3, 1/2, -1, 6] do
        z := (d + B/d)/2; w := (B/d - d)/2;
        u := (w - c1f(vv)/216)/108;
        fval := 11664*u^2 + c1f(vv)*u + c0f(vv);
        assert fval eq z^2;
        Append(~pts, [u, vv, z]);
    end for;
end for;
printf "constructed %o surface points\n", #pts;

nok := 0;
for i in [1..#pts] do
    if nok ge 4 then break; end if;
    pt := pts[i];
    printf "== point %o\n", pt;
    ok := true; E1 := 0; E2 := 0;
    try
        E1, E2 := ZNrModuli(3, 2, pt);
    catch e
        printf "   ZNrModuli FAILED: %o\n", e`Object;
        ok := false;
    end try;
    if not ok then continue; end if;
    printf "   E1 cond %o tors %o\n", Conductor(E1), Invariants(TorsionSubgroup(E1));
    printf "   E2 cond %o tors %o\n", Conductor(E2), Invariants(TorsionSubgroup(E2));
    if IsIsomorphic(E1, E2) then printf "   (isomorphic pair -- skip)\n"; continue; end if;
    bad := 3*Conductor(E1)*Conductor(E2);
    mm := 0; tested := 0;
    for p in PrimesInInterval(5, 300) do
        if bad mod p eq 0 then continue; end if;
        tested +:= 1;
        if (TraceOfFrobenius(E1, p) - TraceOfFrobenius(E2, p)) mod 3 ne 0 then mm +:= 1; end if;
    end for;
    printf "   trace mod-3 mismatches: %o / %o\n", mm, tested;
    g3 := [];
    try
        g3 := Genus2Elliptic3(E1, E2);
    catch e
        printf "   Genus2Elliptic3 error: %o\n", e`Object;
    end try;
    printf "   Genus2Elliptic3 returned %o curve(s)\n", #g3;
    if #g3 gt 0 then
        C := g3[1];
        J := Jacobian(C);
        printf "   glued torsion (first curve): %o\n", Invariants(TorsionSubgroup(J));
        nok +:= 1;
    end if;
end for;
printf "PROBE_DONE ok_glues=%o\n", nok;
quit;
