// claude_z31_ek8_witness.m  (Task B4, step 4)
// Locate the RM witness (1830.2.a.q curve) on the Elkies-Kumar disc-8
// chart (r,s); validate the simplified Igusa-Clebsch formulas; find the
// quadratic twist d linking the EK-reconstructed curve to the witness;
// smoke-test the mod-P Mestre pipeline at the intended sieve primes.
// Source of the IG map: data/claude_z31_ek8_source/igusa8.txt
//   (arXiv:1209.3527 ancillary tarball, directory 8/).
SetColumns(0);
SetMemoryLimit(3*10^9);
Q := Rationals();

// ---------------------------------------------------------------- EK map
function EKIG(r, s)   // works over any field of char != 2,3
    A1 := 2*r*s^2;
    A  := -(9*r*s + 4*r^2 + 4*r + 1)/3;
    B1 := r*s^2*(3*s + 8*r - 2)/3;
    B  := -(54*r^2*s + 81*r*s - 16*r^3 - 24*r^2 - 12*r - 2)/27;
    B2 := r^2;
    return [-24*B1/A1, -12*A, 96*(A/A1)*B1 - 36*B, -4*A1*B2];
end function;

// simplified closed forms; assert equality as rational functions
K2<rr,ss> := RationalFunctionField(Q, 2);
IGraw := EKIG(rr, ss);
IGsimp := [
    -4*(3*ss + 8*rr - 2),
    4*(9*rr*ss + 4*rr^2 + 4*rr + 1),
    -(16/3)*(9*rr*ss + 4*rr^2 + 4*rr + 1)*(3*ss + 8*rr - 2)
        + (4/3)*(54*rr^2*ss + 81*rr*ss - 16*rr^3 - 24*rr^2 - 12*rr - 2),
    -8*rr^3*ss^2 ];
for i in [1..4] do assert IGraw[i] eq IGsimp[i]; end for;
printf "IG_SIMPLIFICATION_OK\n";

// --------------------------------------------------- Mestre conic det tool
function ClebschABCD(ic)
    AP := ic[1]; BP := ic[2]; CP := ic[3]; DP := ic[4];
    A := -AP/120;
    B := (BP + 720*A^2)/6750;
    C := (CP - 8640*A^3 + 108000*A*B)/202500;
    D := (DP + 62208*A^5 - 972000*A^3*B - 1620000*A^2*C
          + 3037500*A*B^2 + 6075000*B*C)/(-4556250);
    return A, B, C, D;
end function;

function ConicDet(ic)
    A, B, C, D := ClebschABCD(ic);
    A11 := 2*C + A*B/3;
    A22 := D;
    A33 := B*D/2 + 2*C*(B^2 + A*C)/9;
    A23 := B*(B^2 + A*C)/3 + C*(2*C + A*B/3)/3;
    A31 := D;
    A12 := 2*(B^2 + A*C)/3;
    return A11*(A22*A33 - A23^2) - A12*(A12*A33 - A23*A31)
           + A31*(A12*A23 - A22*A31);
end function;

// ------------------------------------------------------------- witness IC
P<x> := PolynomialRing(Q);
F := -3356*x^6 + 11364*x^5 - 18347*x^4 + 17202*x^3 - 9863*x^2 + 3264*x - 504;
CW := HyperellipticCurve(F);
ic := IgusaClebschInvariants(CW);
printf "WITNESS_IC %o %o %o %o\n", ic[1], ic[2], ic[3], ic[4];

// ------------------------------------------- solve for (r,s) on the chart
R3<r, s, w> := PolynomialRing(Q, 3);
I2 := -4*(3*s + 8*r - 2);
I4 := 4*(9*r*s + 4*r^2 + 4*r + 1);
I6 := -(16/3)*(9*r*s + 4*r^2 + 4*r + 1)*(3*s + 8*r - 2)
      + (4/3)*(54*r^2*s + 81*r*s - 16*r^3 - 24*r^2 - 12*r - 2);
I10 := -8*r^3*s^2;
// weighted-projective match with witness invariants (ic[1] != 0 checked)
assert ic[1] ne 0;
e1 := I4*ic[1]^2 - ic[2]*I2^2;
e2 := 27*(I6*ic[1]^3 - ic[3]*I2^3);      // clear /27
e3 := I10*ic[1]^5 - ic[4]*I2^5;
// saturate away r=0, s=0 (product locus) via Rabinowitsch
I := ideal< R3 | e1, e2, e3, w*r*s - 1 >;
printf "IDEAL_DIM %o\n", Dimension(I);
V := Variety(I);
printf "N_PREIMAGES %o\n", #V;
goodrs := [];
for v in V do
    r0 := v[1]; s0 := v[2];
    printf "PREIMAGE r=%o s=%o\n", r0, s0;
    // exact weighted-projective check of ALL four invariants
    ig := EKIG(r0, s0);
    lam2 := ig[1]/ic[1];              // lambda^2
    ok := ig[2] eq ic[2]*lam2^2 and 27*ig[3] eq 27*ic[3]*lam2^3
          and ig[4] eq ic[4]*lam2^5;
    printf "PREIMAGE_WEIGHTED_MATCH %o (lam2=%o)\n", ok, lam2;
    if ok then Append(~goodrs, [r0, s0]); end if;
    printf "PREIMAGE_CONIC_DET %o\n", ConicDet(ig);
end for;

// -------------------------------- rebuild over Q + find the twist d
for rs in goodrs do
    r0 := rs[1]; s0 := rs[2];
    C0 := HyperellipticCurveFromIgusaClebsch(EKIG(r0, s0));
    if BaseRing(C0) ne Q then
        printf "REBUILD r=%o s=%o NOT_OVER_Q\n", r0, s0;
        continue;
    end if;
    icc := IgusaClebschInvariants(C0);
    lam2 := icc[1]/ic[1];
    printf "REBUILD_IC_MATCH %o\n",
        icc[2] eq ic[2]*lam2^2 and icc[3] eq ic[3]*lam2^3 and icc[4] eq ic[4]*lam2^5;
    // find quadratic twist d with C0_d isomorphic to witness
    found := false;
    for dd in [1,-1,2,-2,3,-3,5,-5,6,-6,10,-10,15,-15,30,-30,61,-61,122,-122,
               183,-183,305,-305,366,-366,610,-610,915,-915,1830,-1830] do
        Cd := QuadraticTwist(C0, Q!dd);
        if IsIsomorphic(Cd, CW) then
            printf "TWIST_FOUND d=%o\n", dd;
            found := true;
            break;
        end if;
    end for;
    if not found then printf "TWIST_NOT_FOUND in candidate list\n"; end if;
end for;

// ----------------------- smoke: mod-P build + L-poly at the sieve primes
if #goodrs ge 1 then
    r0 := goodrs[1][1]; s0 := goodrs[1][2];
    for p in [103, 211, 1019, 1031, 1039, 1051, 1063, 1087] do
        Fp := GF(p);
        num := Numerator(r0); den := Denominator(r0);
        nums := Numerator(s0); dens := Denominator(s0);
        if IsDivisibleBy(den, p) or IsDivisibleBy(dens, p)
           or Fp!num eq 0 or Fp!nums eq 0 then
            printf "SMOKE p=%o BAD_REDUCTION\n", p;
            continue;
        end if;
        rp := Fp!num/Fp!den; sp := Fp!nums/Fp!dens;
        igp := EKIG(rp, sp);
        dt := ConicDet(igp);
        if igp[4] eq 0 or dt eq 0 then
            printf "SMOKE p=%o DEGEN\n", p;
            continue;
        end if;
        t0 := Cputime();
        Cp := HyperellipticCurveFromIgusaClebsch(igp);
        L := LPolynomial(Cp);
        t1 := Cputime();
        c1 := Evaluate(L, 1); cm1 := Evaluate(L, -1);
        printf "SMOKE p=%o rp=%o sp=%o chi1=%o chim1=%o mask31=[%o,%o] (%.2o s)\n",
            p, rp, sp, c1, cm1,
            c1 mod 31 eq 0 select 1 else 0,
            cm1 mod 31 eq 0 select 1 else 0, t1 - t0;
    end for;
end if;
printf "WITNESS_SCRIPT_DONE\n";
quit;
