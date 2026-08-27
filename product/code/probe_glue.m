// probe_glue.m — measure the true numerical error of the analytic-gluing
// j-invariants against the KNOWN rational values for the [2,2,24]-witness
// pair (M=[0,1;1,0]), at several working precisions.  Diagnoses the
// strict-recognition regression: the strict test needs err < 10^-(2*hq+15)
// with hq = 48, i.e. err < 10^-111 — this prints what the pipeline delivers.
SetColumns(0);
SetMemoryLimit(4*10^9);
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
tru1 := 144833980842901263766603102638411016587229898409/4486263029391955874429830415322275260360;

ZZ := Integers(); QQ := Rationals();
for prec in [120, 200, 300, 400] do
    CC := ComplexField(prec);
    wE := NormPeriods(EA, CC); wF := NormPeriods(EB, CC);
    N := 2; m11 := 0; m12 := 1; m21 := 1; m22 := 0;
    rows := [ [QQ|1,0,0,0], [QQ|0,1,0,0], [QQ|0,0,1,0], [QQ|0,0,0,1],
              [QQ|1/N, 0, m11/N, m12/N], [QQ|0, 1/N, m21/N, m22/N] ];
    H := HermiteForm(Matrix(ZZ, 6, 4, [ [ZZ| N*c : c in r ] : r in rows ]));
    B := [ [ QQ | H[i][j]/N : j in [1..4] ] : i in [1..4] ];
    J := Matrix(QQ, 4,4, [ [ N*EProdPair(B[i], B[j]) : j in [1..4] ] : i in [1..4] ]);
    error if not forall{ <i,j> : i,j in [1..4] | Denominator(J[i][j]) eq 1 }, "not integral";
    JZ := Matrix(ZZ, 4,4, [ [ ZZ!J[i][j] : j in [1..4] ] : i in [1..4] ]);
    error if Abs(Determinant(JZ)) ne 1, "not unimodular";
    F0, T := FrobeniusFormAlternating(JZ);
    Cb := [ [ &+[ QQ | T[i][j]*B[j][k] : j in [1..4] ] : k in [1..4] ] : i in [1..4] ];
    function colv(v, wE, wF, CC)
        return [ CC | v[1]*wE[1] + v[2]*wE[2], v[3]*wF[1] + v[4]*wF[2] ];
    end function;
    cols := [ colv(Cb[i], wE, wF, CC) : i in [1..4] ];
    PA := Matrix(CC, 2,2, [ cols[1][1], cols[2][1], cols[1][2], cols[2][2] ]);
    PB := Matrix(CC, 2,2, [ cols[3][1], cols[4][1], cols[3][2], cols[4][2] ]);
    tau := PA^-1*PB;
    a := Im(tau[1][1]); d := Im(tau[2][2]); b := Im(tau[1][2]);
    if not (a gt 0 and a*d - b^2 gt 0) then tau := PB^-1*PA; end if;
    ros := RosenhainInvariants(tau);
    PCx<xx> := PolynomialRing(CC);
    g := xx*(xx-1)*&*[ xx - r : r in ros ];
    IC := IgusaClebschInvariants(HyperellipticCurve(g));
    j1 := IC[1]^5/IC[4];
    err := Abs(Re(j1) - tru1);
    printf "prec %o: Im(j1) = 10^%o, |Re(j1) - truth| = 10^%o (strict needs 10^-111)\n",
        prec,
        Im(j1) eq 0 select -9999 else Round(Log(10, Abs(Im(j1)))),
        err eq 0 select -9999 else Round(Log(10, err));
    // what BestApproximation delivers at each ladder height
    for hb in [20, 40, 60, 90, 130] do
        if 2*hb + 20 gt prec then break; end if;
        q := BestApproximation(Re(j1), 10^hb);
        hq := Max(Ilog(10, 1+Abs(Numerator(q))), Ilog(10, 1+Abs(Denominator(q))));
        e2 := Abs(Re(j1) - q);
        printf "   hb=%o: hq=%o, |r-q|=10^%o, thr=10^-%o, correct-q=%o\n",
            hb, hq, e2 eq 0 select -9999 else Round(Log(10, e2)), 2*hq+15, q eq tru1;
    end for;
end for;
printf "PROBE_GLUE_DONE\n";
quit;
