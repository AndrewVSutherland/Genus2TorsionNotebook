// claude_z31_ek8_ref.m  (Task B4, validation)
// Reference-vector generator for the Elkies-Kumar disc-8 31-torsion sieve.
// For deterministic pseudo-random (r0,s0) in F_P^2 it builds the genus-2
// curve with the EK Igusa-Clebsch invariants via Magma's Mestre
// implementation and dumps the twist-CLASS data {chi(1), chi(-1)} (multiset;
// the individual values depend on the random twist Magma picks, the multiset
// does not).  TSV columns: P r0 s0 flag chi1 chim1 s1 s2
//   flag OK    - curve built, invariants cross-checked
//   flag DEGEN - I10 = 0 or Mestre conic det = 0  (C side must agree)
//   flag ERR   - Magma error (C side skips)
//   flag MISM  - built curve failed the invariant cross-check (C side skips)
// Run: magma -b claude_z31_ek8_ref.m   (writes refvec_ek8_P<P>.tsv per prime)
SetColumns(0);
SetMemoryLimit(3*10^9);

function EKIG(r, s)
    A1 := 2*r*s^2;
    A  := -(9*r*s + 4*r^2 + 4*r + 1)/3;
    B1 := r*s^2*(3*s + 8*r - 2)/3;
    B  := -(54*r^2*s + 81*r*s - 16*r^3 - 24*r^2 - 12*r - 2)/27;
    B2 := r^2;
    return [-24*B1/A1, -12*A, 96*(A/A1)*B1 - 36*B, -4*A1*B2];
end function;

function ConicDet(ic)
    AP := ic[1]; BP := ic[2]; CP := ic[3]; DP := ic[4];
    A := -AP/120;
    B := (BP + 720*A^2)/6750;
    C := (CP - 8640*A^3 + 108000*A*B)/202500;
    D := (DP + 62208*A^5 - 972000*A^3*B - 1620000*A^2*C
          + 3037500*A*B^2 + 6075000*B*C)/(-4556250);
    A11 := 2*C + A*B/3;
    A22 := D;
    A33 := B*D/2 + 2*C*(B^2 + A*C)/9;
    A23 := B*(B^2 + A*C)/3 + C*(2*C + A*B/3)/3;
    A31 := D;
    A12 := 2*(B^2 + A*C)/3;
    return A11*(A22*A33 - A23^2) - A12*(A12*A33 - A23*A31)
           + A31*(A12*A23 - A22*A31);
end function;

// I10-normalized weighted-projective equality of two IC quadruples
function ICEqual(u, v)
    return u[1]^5*v[4] eq v[1]^5*u[4]
       and u[2]^5*v[4]^2 eq v[2]^5*u[4]^2
       and u[3]^5*v[4]^3 eq v[3]^5*u[4]^3;
end function;

for pr in [<103, 250>, <211, 250>, <1019, 120>, <1031, 120>] do
    P := pr[1]; N := pr[2];
    Fp := GF(P);
    out := Sprintf("refvec_ek8_P%o.tsv", P);
    System("rm -f " cat out);
    nok := 0; ndeg := 0; nerr := 0; nmism := 0;
    T0 := Cputime();
    for k in [1..N] do
        r0 := Fp!(37*k^2 + 11*k + 5);
        s0 := Fp!(53*k^2 + 29*k + 7);
        if r0 eq 0 or s0 eq 0 then continue; end if;   // I10 = 0, trivial
        ig := EKIG(r0, s0);
        if ig[4] eq 0 or ConicDet(ig) eq 0 then
            fprintf out, "%o\t%o\t%o\tDEGEN\t0\t0\t0\t0\n", P, r0, s0;
            ndeg +:= 1;
            continue;
        end if;
        try
            C := HyperellipticCurveFromIgusaClebsch(ig);
            icc := IgusaClebschInvariants(C);
            if not ICEqual(ig, icc) then
                fprintf out, "%o\t%o\t%o\tMISM\t0\t0\t0\t0\n", P, r0, s0;
                nmism +:= 1;
            else
                L := LPolynomial(C);
                c1 := Evaluate(L, 1); cm1 := Evaluate(L, -1);
                s1 := -Coefficient(L, 1); s2 := Coefficient(L, 2);
                fprintf out, "%o\t%o\t%o\tOK\t%o\t%o\t%o\t%o\n",
                    P, r0, s0, c1, cm1, s1, s2;
                nok +:= 1;
            end if;
        catch e
            fprintf out, "%o\t%o\t%o\tERR\t0\t0\t0\t0\n", P, r0, s0;
            nerr +:= 1;
        end try;
    end for;
    printf "REFVEC P=%o ok=%o degen=%o err=%o mism=%o (%.1o s)\n",
        P, nok, ndeg, nerr, nmism, Cputime() - T0;
end for;
printf "REFVEC_DONE\n";
quit;
