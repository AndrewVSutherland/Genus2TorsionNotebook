/* claude_z31_kummer_probe.m
 * Dump ground-truth Kummer vectors from Magma's Jacobian arithmetic for
 * validating the Flynn duplication/biquadratic transcription used by
 * code/claude_z31_kummer_sieve.c.
 * For random curves y^2 = f (deg 6) over GF(p) and random A,B in J:
 *   VEC p f0..f6 kA(4) kB(4) kApB(4) kAmB(4) k2A(4)
 * All vectors are Magma's Eltseq of K!point (projective reps).
 */
SetColumns(0);
SetMemoryLimit(3*10^9);

for p in [101, 211, 1009, 2003, 4001, 10007] do
    SetSeed(555000 + p);
    Fp := GF(p);
    R<x> := PolynomialRing(Fp);
    n := 0;
    tries := 0;
    while n lt 10 and tries lt 500 do
        tries +:= 1;
        f := R![Random(Fp) : i in [0..6]];
        if Degree(f) ne 6 then continue; end if;
        ok := true;
        kA := []; kB := []; kApB := []; kAmB := []; k2A := [];
        try
            C := HyperellipticCurve(f);
            Jc := Jacobian(C);
            K := KummerSurface(Jc);
            A := Random(Jc);
            B := Random(Jc);
            kA := Eltseq(K!A);
            kB := Eltseq(K!B);
            kApB := Eltseq(K!(A+B));
            kAmB := Eltseq(K!(A-B));
            k2A := Eltseq(K!(A+A));
        catch e
            ok := false;
        end try;
        if not ok then continue; end if;
        n +:= 1;
        co := Coefficients(f);
        printf "VEC %o", p;
        for i in [1..7] do printf " %o", Integers()!co[i]; end for;
        for s in [kA, kB, kApB, kAmB, k2A] do
            for i in [1..4] do printf " %o", Integers()!s[i]; end for;
        end for;
        printf "\n";
    end while;
    printf "PRIMEDONE %o n=%o tries=%o\n", p, n, tries;
end for;
print "ALLDONE";
quit;
