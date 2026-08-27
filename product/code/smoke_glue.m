// smoke_glue.m — validate the gluing + funnel pipeline on the 5-isogenous
// [10]+[10] pair 66.c3 / 66.c4 (candidate for split [10,10] or [5,10]).
// NOTE (2026-08-25, public-snapshot review): this pair does not appear to
// have compatible 2-torsion fields (Q(sqrt 33) vs Q(sqrt -2)), so the
// Genus2Elliptic2 call below fails its same-2-torsion-field require -- it
// always did, in the original genus2.m as well as in gluing.m.  The 66-class
// gluing that succeeded went through the analytic route instead (see
// product/logs/glue66.log).
SetColumns(0);
SetMemoryLimit(3*10^9);
spec := GetEnv("TORSION_MAGMA_SPEC"); if spec eq "" then spec := "/home/claude/Magma/magma.spec"; end if; AttachSpec(spec);
Attach("gluing.m");

E1 := EllipticCurve([1,0,0,-45,81]); E2 := EllipticCurve([1,0,0,115,561]);  // 66.c3, 66.c4
printf "E1 tors %o  E2 tors %o\n", Invariants(TorsionSubgroup(E1)), Invariants(TorsionSubgroup(E2));

L := Genus2Elliptic2(E1,E2);
printf "num glued curves: %o\n", #L;

for i in [1..#L] do
    C := L[i];
    f0,h0 := HyperellipticPolynomials(C);
    g := 4*f0 + h0^2;   // y^2 = g(x)
    printf "GLUED %o: g = %o\n", i, g;
    disc := Integers()!Discriminant(g);
    // profile at 8 good odd primes
    p := 2; cnt := 0;
    ords := [];
    while cnt lt 8 do
        p := NextPrime(p);
        if disc mod p ne 0 then
            cnt +:= 1;
            gp := PolynomialRing(GF(p))!g;
            Jp := Jacobian(HyperellipticCurve(gp));
            Ap := AbelianGroup(Jp);
            printf "  p=%o  invs=%o\n", p, Invariants(Ap);
            Append(~ords, #Ap);
        end if;
    end while;
    printf "  gcd of orders: %o\n", GCD(ords);
    t0 := Cputime();
    T := TorsionSubgroup(Jacobian(C));
    printf "  EXACT torsion: %o  (%.1o s)\n", Invariants(T), Cputime()-t0;
end for;
quit;
