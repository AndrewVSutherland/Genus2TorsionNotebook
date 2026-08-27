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
ok, m := IsIsogenous(EA, EB);
printf "connecting degree data: %o %o\n", ok, m;
MS := [ [0,1,1,0],[0,2,3,0],[0,3,2,0],[0,4,4,0],[1,0,0,4],[2,0,0,2],[3,0,0,3],[4,0,0,1] ];
N := 5;
prec := 500;
CC := ComplexField(prec);
wE := NormPeriods(EA, CC); wF := NormPeriods(EB, CC);
ZZ := Integers(); QQ := Rationals();
for Mv in MS do
    m11 := Mv[1]; m12 := Mv[2]; m21 := Mv[3]; m22 := Mv[4];
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
    okr := true; ros := [];
    try ros := RosenhainInvariants(tau); catch e okr := false; end try;
    if not okr then printf "M=%o rosenhain failed\n", Mv; continue; end if;
    PCx<xx> := PolynomialRing(CC);
    g := xx*(xx-1)*&*[ xx - r : r in ros ];
    IC := IgusaClebschInvariants(HyperellipticCurve(g));
    if Abs(IC[4]) lt 10.0^(-60) then printf "M=%o DEGENERATE (I10~0)\n", Mv; continue; end if;
    j1 := IC[1]^5/IC[4];
    r := Re(j1);
    // strict: try denominator bounds well below prec/2 and demand err << 1/q^2
    found := false;
    for hb in [30, 60, 100, 160] do
        q := BestApproximation(r, 10^hb);
        hq := Max(Ilog(10, 1+Abs(Numerator(q))), Ilog(10, 1+Denominator(q)));
        if Abs(r - q) lt RealField(20)!10.0^(-(2*hq + 20)) then
            printf "M=%o RATIONAL j1 (height %o digits): STABLE\n", Mv,
                Ilog(10,1+Abs(Numerator(q))) + hq;
            found := true;
            break;
        end if;
    end for;
    if not found then printf "M=%o irrational (real) j1\n", Mv; end if;
end for;
printf "PIN90_DONE\n";
quit;
