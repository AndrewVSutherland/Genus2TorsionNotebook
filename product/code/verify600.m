SetColumns(0);
SetMemoryLimit(8*10^9);
// run from product/code (house convention); set TORSION_JAC_ROOT to
// run from elsewhere
root := GetEnv("TORSION_JAC_ROOT");
if root ne "" then ChangeDirectory(root cat "/product/code"); end if;
load "split_lab.m";
load "analytic_glue.m";
D := CremonaDatabase();
cls := EllipticCurves(D, 90, 3);
Ts := [ Invariants(TorsionSubgroup(E)) : E in cls ];
i12 := [ k : k in [1..#cls] | Ts[k] eq [12] ][1];
i26 := [ k : k in [1..#cls] | Ts[k] eq [2,6] ][1];
EA := cls[i12]; EB := cls[i26];
// recompute the M=[0,4;4,2] gluing at two precisions and compare
for prec in [300, 700] do
    CC := ComplexField(prec);
    wE := NormPeriods(EA, CC); wF := NormPeriods(EB, CC);
    ZZ := Integers(); QQ := Rationals();
    N := 5; m11 := 0; m12 := 4; m21 := 4; m22 := 2;
    rows := [ [QQ|1,0,0,0], [QQ|0,1,0,0], [QQ|0,0,1,0], [QQ|0,0,0,1],
              [QQ|1/N, 0, m11/N, m12/N], [QQ|0, 1/N, m21/N, m22/N] ];
    H := HermiteForm(Matrix(ZZ, 6, 4, [ [ZZ| N*c : c in r ] : r in rows ]));
    B := [ [ QQ | H[i][j]/N : j in [1..4] ] : i in [1..4] ];
    JZ := Matrix(ZZ, 4,4, [ [ ZZ!(N*EProdPair(B[i], B[j])) : j in [1..4] ] : i in [1..4] ]);
    F0, T := FrobeniusFormAlternating(JZ);
    Cb := [ [ &+[ QQ | T[i][j]*B[j][k] : j in [1..4] ] : k in [1..4] ] : i in [1..4] ];
    cols := [ [ CC | v[1]*wE[1] + v[2]*wE[2], v[3]*wF[1] + v[4]*wF[2] ] where v := Cb[i] : i in [1..4] ];
    PA := Matrix(CC, 2,2, [ cols[1][1], cols[2][1], cols[1][2], cols[2][2] ]);
    PB := Matrix(CC, 2,2, [ cols[3][1], cols[4][1], cols[3][2], cols[4][2] ]);
    tau := PA^-1*PB;
    if Im(tau[1][1]) lt 0 then tau := PB^-1*PA; end if;
    ros := RosenhainInvariants(tau);
    PCx<xx> := PolynomialRing(CC);
    g := xx*(xx-1)*&*[ xx - r : r in ros ];
    IC := IgusaClebschInvariants(HyperellipticCurve(g));
    j1 := IC[1]^5/IC[4]; j2 := IC[1]^3*IC[2]/IC[4]; j3 := IC[1]^2*IC[3]/IC[4];
    // proper verified rational recognition: err < 10^-(2*hq+10)
    for pair in [* <"j1", j1>, <"j2", j2>, <"j3", j3> *] do
        r := Re(pair[2]);
        q := BestApproximation(r, 10^((prec-30) div 2));
        hq := Ilog(10, 1+Abs(Numerator(q))) + Ilog(10, 1+Denominator(q));
        err := Abs(r - q);
        okv := err lt RealField(20)!10.0^(-(2*Ilog(10,1+Denominator(q)) + 10));
        printf "prec=%o %o: height %o digits, err 10^%o, VERIFIED=%o\n",
            prec, pair[1], hq, err eq 0 select -prec else Ilog(10, Ceiling(1/err)) * -1, okv;
        if okv then printf "   q = %o\n", q; end if;
    end for;
end for;
printf "VERIFY600_DONE\n";
quit;
