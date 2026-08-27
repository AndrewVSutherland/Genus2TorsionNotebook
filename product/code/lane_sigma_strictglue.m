// lane_sigma_strictglue.m — GlueScan + STRICT rationality on sieve survivors
// (2026-08-14 plan, section 2, step 2).  For each candidate pair in CANDS and
// each N in NLIST: run GlueScan (weak recognition = prefilter only), then
// RECOMPUTE the complex absolute invariants for every returned anti-isometry
// matrix and accept them as rational only under the strict test
//   err < 10^-(2*height(q)+15),  heights <= (prec-30)/2
// (glue_window.m idiom; analytic_glue.m's RatApprox is the weak |err|<eps
// test that accepts essentially every real number — 08-13 lesson).  Verified
// triples with heights <= MaxMestreH digits get one capped Mestre attempt
// (ReduceIC on the I2=1-normalized tuple, then PARI hyperellred; NEVER
// ReducedMinimalWeierstrassModel on Mestre output); the invariant triple is
// logged REGARDLESS of Mestre success (a certified Q-surface with torsion
// injection is reportable without the curve model).
// CANDS ships with one placeholder: the [70] control pair x7(2/11) x
// x10(-1/2) (3-congruent, non-isogenous; algebraic control in x7x10_glue.m).
// Smoke-test it with NLIST:="3"; real [2,6]x[2,6] survivors from
// lane_cong_sieve.m Fam:="26" (NONISO lines) replace it with NLIST:="5,7".
// Usage: cd product/code && magma -b lane_sigma_strictglue.m > ../logs/strictglue.log
//   optional: Prec:=<int> (default 200), NLIST:="5,7", MaxMestreH:=<int> (60),
//             MemGB:=<int> (default 6)
SetColumns(0);
if not assigned Prec then Prec := 200; elif Type(Prec) eq MonStgElt then Prec := StringToInteger(Prec); end if;
if not assigned NLIST then NLIST := "5,7"; end if;
if not assigned MaxMestreH then MaxMestreH := 60; elif Type(MaxMestreH) eq MonStgElt then MaxMestreH := StringToInteger(MaxMestreH); end if;
if not assigned MemGB then MemGB := 6; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
load "split_lab.m";
load "analytic_glue.m";

RQ := Rationals(); ZZ := Integers();
NLISTI := [ StringToInteger(x) : x in Split(NLIST, ",") ];

// candidate pairs as <aInvariants(EA), aInvariants(EB)> (any model; minimal
// models are taken at runtime).  Placeholder = the [70] control pair.
u7 := RQ!2/11; u10 := RQ!-1/2;
b7 := u7^3 - u7^2;  c7 := u7^2 - u7;
b10 := u10^3*(u10-1)*(2*u10-1)/(u10^2-3*u10+1)^2;
c10 := -u10*(u10-1)*(2*u10-1)/(u10^2-3*u10+1);
CANDS := [*
    < [ RQ | 1-c7, -b7, -b7, 0, 0 ], [ RQ | 1-c10, -b10, -b10, 0, 0 ] >
*];

function StrictRat(z, prec)
    // strict rational recognition (glue_window.m idiom): laddered
    // BestApproximation; accept only if err < 10^-(2*height+15), so the
    // required precision genuinely certifies the height reached.
    if Abs(Im(z)) gt RealField(20)!10.0^(-(prec div 3)) then return false, RQ!0, 0; end if;
    r := Re(z);
    for hb in [20, 40, 60, 90] do
        if 2*hb + 20 gt prec then break; end if;
        q := BestApproximation(r, 10^hb);
        hq := Max(Ilog(10, 1+Abs(Numerator(q))), Ilog(10, 1+Denominator(q)));
        if Abs(r - q) lt RealField(20)!10.0^(-(2*hq + 15)) then
            return true, q, hq;
        end if;
    end for;
    return false, RQ!0, 0;
end function;

function InvariantsForM(E, F, N, Mv, prec)
    // recompute the COMPLEX absolute invariants for one anti-isometry matrix
    // (GlueScan discards the reals after its weak recognition; the strict
    // test needs them).  Same lattice pipeline as GlueScan/GlueCandidates.
    CC := ComplexField(prec);
    wE := NormPeriods(E, CC); wF := NormPeriods(F, CC);
    QQ := Rationals();
    m11 := Mv[1]; m12 := Mv[2]; m21 := Mv[3]; m22 := Mv[4];
    rows := [ [QQ|1,0,0,0], [QQ|0,1,0,0], [QQ|0,0,1,0], [QQ|0,0,0,1],
              [QQ|1/N, 0, m11/N, m12/N], [QQ|0, 1/N, m21/N, m22/N] ];
    H := HermiteForm(Matrix(ZZ, 6, 4, [ [ZZ| N*c : c in r ] : r in rows ]));
    B := [ [ QQ | H[i][j]/N : j in [1..4] ] : i in [1..4] ];
    J := Matrix(QQ, 4,4, [ [ N*EProdPair(B[i], B[j]) : j in [1..4] ] : i in [1..4] ]);
    if not forall{ <i,j> : i,j in [1..4] | Denominator(J[i][j]) eq 1 } then return false, [CC|]; end if;
    JZ := Matrix(ZZ, 4,4, [ [ ZZ!J[i][j] : j in [1..4] ] : i in [1..4] ]);
    if Abs(Determinant(JZ)) ne 1 then return false, [CC|]; end if;
    F0, T := FrobeniusFormAlternating(JZ);
    if F0[1][3] ne 1 or F0[2][4] ne 1 then return false, [CC|]; end if;
    Cb := [ [ &+[ QQ | T[i][j]*B[j][k] : j in [1..4] ] : k in [1..4] ] : i in [1..4] ];
    cols := [ [ CC | v[1]*wE[1] + v[2]*wE[2], v[3]*wF[1] + v[4]*wF[2] ] where v := Cb[i] : i in [1..4] ];
    PA := Matrix(CC, 2,2, [ cols[1][1], cols[2][1], cols[1][2], cols[2][2] ]);
    PB := Matrix(CC, 2,2, [ cols[3][1], cols[4][1], cols[3][2], cols[4][2] ]);
    okt := true; tau := 0;
    try tau := PA^-1*PB; catch e okt := false; end try;
    if not okt then return false, [CC|]; end if;
    function imposdef(t)
        a := Im(t[1][1]); d := Im(t[2][2]); b := Im(t[1][2]);
        return a gt 0 and a*d - b^2 gt 0;
    end function;
    if not imposdef(tau) then
        try tau := PB^-1*PA; catch e okt := false; end try;
        if not okt or not imposdef(tau) then return false, [CC|]; end if;
    end if;
    okr := true; ros := [];
    try ros := RosenhainInvariants(tau); catch e okr := false; end try;
    if not okr then return false, [CC|]; end if;
    if exists{ r : r in ros | Abs(r) lt 10.0^(-8) or Abs(r-1) lt 10.0^(-8) } then return false, [CC|]; end if;
    PCx<xx> := PolynomialRing(CC);
    g := xx*(xx-1)*&*[ xx - r : r in ros ];
    IC := IgusaClebschInvariants(HyperellipticCurve(g));
    if Abs(IC[4]) lt 10.0^(-prec div 2) then return false, [CC|]; end if;
    return true, [ CC | IC[1]^5/IC[4], IC[1]^3*IC[2]/IC[4], IC[1]^2*IC[3]/IC[4] ];
end function;

printf "STRICTGLUE_START pairs=%o Ns=%o prec=%o\n", #CANDS, NLISTI, Prec;
nver := 0;
t0 := Cputime();
for idx in [1..#CANDS] do
    cd := CANDS[idx];
    EA := MinimalModel(EllipticCurve([ RQ | x : x in cd[1] ]));
    EB := MinimalModel(EllipticCurve([ RQ | x : x in cd[2] ]));
    T1 := Invariants(TorsionSubgroup(EA)); T2 := Invariants(TorsionSubgroup(EB));
    printf "PAIR idx=%o %o x %o torsions %o %o\n", idx, aInvariants(EA), aInvariants(EB), T1, T2;
    for N in NLISTI do
        scan := [* *];
        try scan := GlueScan(EA, EB, N : prec := Prec); catch e printf "SCANERR idx=%o N=%o\n", idx, N; end try;
        printf "SCAN idx=%o N=%o candidates=%o\n", idx, N, #scan;
        for s in scan do
            okI, js := InvariantsForM(EA, EB, N, s[1], Prec);
            if not okI then printf "STRICTRAT idx=%o N=%o M=%o RECOMPUTE-FAILED\n", idx, N, s[1]; continue; end if;
            qs := [ RQ | ]; hs := [ ZZ | ]; okall := true;
            for j in js do
                okq, q, hq := StrictRat(j, Prec);
                if not okq then okall := false; break; end if;
                Append(~qs, q); Append(~hs, hq);
            end for;
            if not okall then
                printf "STRICTRAT idx=%o N=%o M=%o REJECTED (weak-only)\n", idx, N, s[1];
                continue;
            end if;
            nver +:= 1;
            printf "STRICTRAT idx=%o N=%o M=%o VERIFIED heights=%o\n", idx, N, s[1], hs;
            printf "STRICTINV idx=%o N=%o q1=%o q2=%o q3=%o\n", idx, N, qs[1], qs[2], qs[3];
            hsum := [ Ilog(10,1+Abs(Numerator(q))) + Ilog(10,1+Denominator(q)) : q in qs ];
            if qs[1] eq 0 then printf "MESTRESKIP idx=%o I2=0 (needs other normalization)\n", idx; continue; end if;
            if Max(hsum) gt MaxMestreH then
                printf "MESTRESKIP idx=%o N=%o heights=%o > %o (invariants recorded above)\n", idx, N, hsum, MaxMestreH;
                continue;
            end if;
            ICq := ReduceIC([ RQ | 1, qs[2]/qs[1], qs[3]/qs[1], 1/qs[1] ]);
            printf "MESTREIN idx=%o reduced IC heights %o\n", idx,
                [ Ilog(10, 1+Abs(Numerator(x))) + Ilog(10, 1+Denominator(x)) : x in ICq ];
            okc := true; C0 := 0;
            tM := Cputime();
            try C0 := HyperellipticCurveFromIgusaClebsch(ICq); catch e okc := false; end try;
            if not okc then
                printf "MESTREFAIL idx=%o N=%o M=%o (invariants recorded above)\n", idx, N, s[1];
                continue;
            end if;
            try C0 := hyperellred(C0); catch e; end try;
            printf "MESTREOK idx=%o N=%o M=%o %o s C: %o\n", idx, N, s[1], Cputime()-tM, HyperellipticPolynomials(C0);
        end for;
    end for;
    System(Sprintf("echo 'strictglue pair %o/%o done %o verified' >> ../logs/strictglue.progress", idx, #CANDS, nver));
end for;
printf "STRICTGLUE_DONE n=%o verified=%o %o s\n", #CANDS, nver, Cputime()-t0;
quit;
