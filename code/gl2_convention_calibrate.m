// 37-hunt: convention calibration on the KNOWN case J_0(23) = Jac(X_0(23)).
// The pipeline (modular symbols -> lattice -> intersection pairing ->
// periods -> symplectic basis -> big period matrix -> ReconstructCurve)
// must reproduce X_0(23): y^2 = x^6-8x^5+2x^4+2x^3-11x^2+10x-7.
// All convention switches are iterated:
//   frame in {basis, lattice=Cm*PM}, block in {[A|B],[B|A]}, sign in {+,-}
// The switch combination that passes IsBigPeriodMatrix AND reproduces the
// right G2 invariants is THE convention for the 2190 driver.
//
// Run: magma -b code/gl2_convention_calibrate.m > results/gl2_convention_cal.log

SetColumns(0);
SetSeed(1);
SetMemoryLimit(8*10^9);

Attach("~/.claude/jobs/a1db5dd4/tmp/polredabs_shim.m");
AttachSpec("~/.claude/jobs/a1db5dd4/tmp/endomorphisms/endomorphisms/magma/spec");
AttachSpec("~/.claude/jobs/a1db5dd4/tmp/quartic/magma/spec");
AttachSpec("~/.claude/jobs/a1db5dd4/tmp/curve_reconstruction/magma/spec");

Prec := 80;
SetDefaultRealField(RealField(Prec + 10));
FX := RationalsExtra(Prec);
CC := FX`CC;

Q := Rationals();
Px<x> := PolynomialRing(Q);
X23 := HyperellipticCurve(x^6-8*x^5+2*x^4+2*x^3-11*x^2+10*x-7);
g2target := G2Invariants(X23);

M := ModularSymbols(23, 2, 0);
S := CuspidalSubspace(M);
Af := S;   // whole thing: J0(23), dim 4 sign-0
BAf := ChangeRing(BasisMatrix(VectorSpace(Af)), Q);
Lb := ChangeRing(BasisMatrix(Lattice(Af)), Q);
Cm := Solution(BAf, Lb);
IP := ChangeRing(IntersectionPairing(Af), Q);
E := Cm * IP * Transpose(Cm);
den := LCM([Denominator(v) : v in Eltseq(E)]);
EZ0 := Matrix(Integers(), 4, 4, [Integers()!(den*v) : v in Eltseq(E)]);
g0 := GCD([v : v in Eltseq(EZ0) | v ne 0]);
EZ := EZ0 div g0;
printf "TYPE %o\n", ElementaryDivisors(EZ);

nterms := Ceiling(23 * Prec * Log(10)/(4*Pi(RealField(20)))) + 200;
Pv := Periods(Af, nterms);
PMb := Matrix(CC, 4, 2, &cat[ [CC!Pv[i][1], CC!Pv[i][2]] : i in [1..4] ]);

for frame in [1, 2] do
    PM := frame eq 1 select PMb else ChangeRing(Cm, CC)*PMb;
    for sgn in [1, -1] do
        EpN := sgn*EZ;
        ok, T := true, IdentityMatrix(Integers(),4);
        try
            _, T := FrobeniusFormAlternating(EpN);
        catch e ok := false; end try;
        if not ok then printf "frame=%o sgn=%o : Frobenius failed\n", frame, sgn; continue; end if;
        Bs := ChangeRing(ChangeRing(T, Q), CC);
        Om := Bs * PM;
        for blk in [1, 2] do
            if blk eq 1 then
                bigP := Matrix(CC, 2, 4, [ Om[i][j] : i in [1..4], j in [1..2] ]);
            else
                ord := [3,4,1,2];
                bigP := Matrix(CC, 2, 4, [ Om[ord[i]][j] : i in [1..4], j in [1..2] ]);
            end if;
            isb := false;
            try isb := IsBigPeriodMatrix(bigP); catch e ; end try;
            printf "frame=%o sgn=%o blk=%o : IsBigPeriodMatrix=%o\n", frame, sgn, blk, isb;
            if not isb then continue; end if;
            try
                Crec := ReconstructCurve(bigP, FX : Base := true);
                gg := G2Invariants(ChangeRing(Crec, Q));
                match := gg eq g2target;
                printf "  RECONSTRUCTED over %o ; G2_MATCH_X023 = %o\n", BaseRing(Crec), match;
                if match then
                    printf "CONVENTION_FOUND frame=%o sgn=%o blk=%o\n", frame, sgn, blk;
                end if;
            catch e
                msg := Sprint(e`Object);
                printf "  RECFAIL %o\n", #msg gt 150 select Substring(msg,1,150) else msg;
            end try;
        end for;
    end for;
end for;
printf "CALIBRATION_DONE\n";
quit;
