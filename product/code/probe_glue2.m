// probe_glue2.m — locate the 10^-80 accuracy ceiling: compare Rosenhain
// invariants and periods of the [2,2,24]-witness gluing computed at prec
// 200 vs 500; print the parent precisions and the cross-precision drift of
// each stage (periods, tau, ros, IC, j1).
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

ZZ := Integers(); QQ := Rationals();
function Stage(prec, EA, EB)
    CC := ComplexField(prec);
    wE := NormPeriods(EA, CC); wF := NormPeriods(EB, CC);
    N := 2; m11 := 0; m12 := 1; m21 := 1; m22 := 0;
    rows := [ [QQ|1,0,0,0], [QQ|0,1,0,0], [QQ|0,0,1,0], [QQ|0,0,0,1],
              [QQ|1/N, 0, m11/N, m12/N], [QQ|0, 1/N, m21/N, m22/N] ];
    H := HermiteForm(Matrix(ZZ, 6, 4, [ [ZZ| N*c : c in r ] : r in rows ]));
    B := [ [ QQ | H[i][j]/N : j in [1..4] ] : i in [1..4] ];
    J := Matrix(QQ, 4,4, [ [ N*EProdPair(B[i], B[j]) : j in [1..4] ] : i in [1..4] ]);
    JZ := Matrix(ZZ, 4,4, [ [ ZZ!J[i][j] : j in [1..4] ] : i in [1..4] ]);
    F0, T := FrobeniusFormAlternating(JZ);
    Cb := [ [ &+[ QQ | T[i][j]*B[j][k] : j in [1..4] ] : k in [1..4] ] : i in [1..4] ];
    cols := [ [ CC | v[1]*wE[1] + v[2]*wE[2], v[3]*wF[1] + v[4]*wF[2] ]
              where v := Cb[i] : i in [1..4] ];
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
    return wE, tau, ros, j1;
end function;

wE2, tau2, ros2, j2 := Stage(200, EA, EB);
wE5, tau5, ros5, j5 := Stage(500, EA, EB);
r2s := [ x : x in ros2 ]; r5s := [ x : x in ros5 ];
printf "parent precisions: w %o, tau %o, j %o\n",
    Precision(Parent(wE2[1])), Precision(Parent(tau2[1][1])), Precision(Parent(j2));
C5 := ComplexField(500);
dW := Abs(C5!wE5[1] - C5!wE2[1]);
dT := Abs(C5!tau5[1][1] - C5!tau2[1][1]);
// ros orderings may differ between runs; match each prec-200 value to its
// nearest prec-500 partner
dR := Max([ Min([ Abs(C5!a - C5!b) : b in r5s ]) : a in r2s ]);
dJ := Abs(C5!j5 - C5!j2);
lg := func<x | x eq 0 select -9999 else Round(Log(10, x))>;
printf "cross-precision drift: dW=10^%o dTau=10^%o dRos=10^%o dJ=10^%o\n", lg(dW), lg(dT), lg(dR), lg(dJ);
tru1 := 144833980842901263766603102638411016587229898409/4486263029391955874429830415322275260360;
printf "|j(500) - truth| = 10^%o\n", lg(Abs(C5!j5 - C5!tru1));
printf "PROBE_GLUE2_DONE\n";
quit;
