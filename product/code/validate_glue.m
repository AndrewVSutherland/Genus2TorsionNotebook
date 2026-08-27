// validate_glue.m — validate the analytic gluer at N=2 against
// Genus2Elliptic2 on the [2,2,24]-witness pair E26(241/81) x E28(1/3).
// Two independent gates (rewritten 2026-08-15, codex round 2):
//   (1) RECOGNITION-FREE: the six analytic anti-isometries' complex
//       (j1,j2,j3) triples must match the six algebraic gluings' EXACT
//       invariant triples bijectively to ~working accuracy (< 10^-150 at
//       prec 200).  This is the actual soundness proof of the lattice ->
//       tau -> Rosenhain -> Igusa-Clebsch pipeline.
//   (2) STRICT RECOGNITION: GlueScan at prec 600 (the true invariant
//       heights here are 100-182 digits; StrictRatApprox needs
//       prec >= 2*height + 35) must return EXACTLY the algebraic rational
//       triples.
// History: the 2026-08-13 run (weak RatApprox, prec 120) printed a
// height-47 "rational invariant" for M=[0,1;1,0]; that value was a
// best-approximant artifact of the weak test — the true invariant of that
// gluing is the height-182 rational of algebraic gluing #6 (10^-80 away).
// The old log's GLUE line is superseded by this validation.
// No Mestre reconstruction here: invariant-triple identity IS the
// validation, and Mestre stalls at these heights (documented 08-13).
// Usage: cd product/code && magma -b validate_glue.m > ../logs/validate_glue.log
SetColumns(0);
SetMemoryLimit(6*10^9);
load "split_lab.m";
load "analytic_glue.m";

RQ := Rationals();

function E26(t)
    r1 := (-2*t+10)/((t+3)*(t-3));
    r2 := (-t^3+7*t^2-11*t+5)/(4*(t+3)*(t-3)^2);
    r3 := (-2*t^2+4*t-2)/((t+3)^2*(t-3));
    f := (RQx.1 - r1)*(RQx.1 - r2)*(RQx.1 - r3);
    return EllipticCurve(RQx!f);
end function;
function E28(u)
    r1 := (16*u^3+12*u^2+2*u)/(8*u^2-1)^2;
    r2 := (32*u^3+24*u^2+8*u+1)/(16*u^2*(8*u^2-1));
    r3 := (-32*u^4-32*u^3-12*u^2-2*u)/((4*u+1)^2*(8*u^2-1));
    f := (RQx.1 - r1)*(RQx.1 - r2)*(RQx.1 - r3);
    return EllipticCurve(RQx!f);
end function;

EA := E26(241/81); EB := E28(1/3);
ref := Genus2Elliptic2(EA, EB);
printf "reference: %o algebraic gluings\n", #ref;
algJ := [];
for k in [1..#ref] do
    IC := IgusaClebschInvariants(ref[k]);
    Append(~algJ, [ IC[1]^5/IC[4], IC[1]^3*IC[2]/IC[4], IC[1]^2*IC[3]/IC[4] ]);
    printf "alg %o: invariant heights %o %o %o\n", k,
        Ilog(10, 1+Abs(Numerator(algJ[k][1]))),
        Ilog(10, 1+Abs(Numerator(algJ[k][2]))),
        Ilog(10, 1+Abs(Numerator(algJ[k][3])));
end for;

// ---- gate (1): recognition-free bijective match at prec 200 ----
prec1 := 200;
CC := ComplexField(prec1);
ZZ := Integers(); QQ := Rationals();
wE := NormPeriods(EA, CC); wF := NormPeriods(EB, CC);
N := 2;
lg := func<x | x eq 0 select -9999 else Round(Log(10, x))>;
matched := {};
nanal := 0;
for m11, m12, m21, m22 in [0..N-1] do
    if (m11*m22 - m12*m21 + 1) mod N ne 0 then continue; end if;
    M := Matrix(Integers(N), 2, 2, [m11,m12,m21,m22]);
    if not IsInvertible(M) then continue; end if;
    rows := [ [QQ|1,0,0,0], [QQ|0,1,0,0], [QQ|0,0,1,0], [QQ|0,0,0,1],
              [QQ|1/N, 0, m11/N, m12/N], [QQ|0, 1/N, m21/N, m22/N] ];
    H := HermiteForm(Matrix(ZZ, 6, 4, [ [ZZ| N*c : c in r ] : r in rows ]));
    B := [ [ QQ | H[i][j]/N : j in [1..4] ] : i in [1..4] ];
    J := Matrix(QQ, 4,4, [ [ N*EProdPair(B[i], B[j]) : j in [1..4] ] : i in [1..4] ]);
    if not forall{ <i,j> : i,j in [1..4] | Denominator(J[i][j]) eq 1 } then continue; end if;
    JZ := Matrix(ZZ, 4,4, [ [ ZZ!J[i][j] : j in [1..4] ] : i in [1..4] ]);
    if Abs(Determinant(JZ)) ne 1 then continue; end if;
    F0, T := FrobeniusFormAlternating(JZ);
    if F0[1][3] ne 1 or F0[2][4] ne 1 then continue; end if;
    Cb := [ [ &+[ QQ | T[i][j]*B[j][k] : j in [1..4] ] : k in [1..4] ] : i in [1..4] ];
    cols := [ [ CC | v[1]*wE[1] + v[2]*wE[2], v[3]*wF[1] + v[4]*wF[2] ]
              where v := Cb[i] : i in [1..4] ];
    PA := Matrix(CC, 2,2, [ cols[1][1], cols[2][1], cols[1][2], cols[2][2] ]);
    PB := Matrix(CC, 2,2, [ cols[3][1], cols[4][1], cols[3][2], cols[4][2] ]);
    okt := true; tau := 0;
    try tau := PA^-1*PB; catch e okt := false; end try;
    if not okt then continue; end if;
    a := Im(tau[1][1]); d := Im(tau[2][2]); b := Im(tau[1][2]);
    if not (a gt 0 and a*d - b^2 gt 0) then
        try tau := PB^-1*PA; catch e okt := false; end try;
        if not okt then continue; end if;
        a := Im(tau[1][1]); d := Im(tau[2][2]); b := Im(tau[1][2]);
        if not (a gt 0 and a*d - b^2 gt 0) then continue; end if;
    end if;
    okr := true; ros := [];
    try ros := RosenhainInvariants(tau); catch e okr := false; end try;
    if not okr then continue; end if;
    PCx<xx> := PolynomialRing(CC);
    g := xx*(xx-1)*&*[ xx - r : r in ros ];
    IC := IgusaClebschInvariants(HyperellipticCurve(g));
    if Abs(IC[4]) lt 10.0^(-prec1 div 2) then continue; end if;
    j1 := IC[1]^5/IC[4]; j2 := IC[1]^3*IC[2]/IC[4]; j3 := IC[1]^2*IC[3]/IC[4];
    nanal +:= 1;
    best := Re(CC!10)^100; bestk := 0;
    for k in [1..#algJ] do
        d1 := Abs(j1 - CC!algJ[k][1]) + Abs(j2 - CC!algJ[k][2]) + Abs(j3 - CC!algJ[k][3]);
        if Re(d1) lt best then best := Re(d1); bestk := k; end if;
    end for;
    ok := best lt Re(CC!10)^(-150);
    printf "MATCH M=[%o,%o;%o,%o] -> alg %o at distance 10^%o %o\n",
        m11, m12, m21, m22, bestk, lg(best), ok select "PASS" else "FAIL";
    if ok then Include(~matched, bestk); end if;
end for;
printf "VALIDATE_GLUE_MATCH %o/%o analytic gluings matched, %o/%o algebraic gluings covered\n",
    #matched, nanal, #matched, #algJ;
error if #matched ne #algJ, "recognition-free match FAILED";

// ---- gate (2): strict recognition at prec 600 returns the exact rationals ----
scan := GlueScan(EA, EB, 2 : prec := 600, Verbose := false);
printf "strict GlueScan at prec 600: %o rational triples\n", #scan;
nexact := 0;
for s in scan do
    hit := exists(k){ k : k in [1..#algJ] | algJ[k][1] eq s[2] and algJ[k][2] eq s[3] and algJ[k][3] eq s[4] };
    if hit then nexact +:= 1;
    else printf "STRICT M=%o returned a triple NOT among the algebraic invariants\n", s[1];
    end if;
end for;
printf "VALIDATE_GLUE_STRICT %o/%o strict recognitions exactly equal algebraic invariants\n", nexact, #scan;
error if nexact ne #algJ or #scan ne #algJ, "strict recognition validation FAILED";
printf "VALIDATE_GLUE_DONE\n";
quit;
