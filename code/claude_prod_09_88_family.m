// (8,8) lane, step 1.0: validation gate for the rebuilt Nicholls Lambda_334
// family (definitions in code/claude_prod_09_88_defs.m; the original
// scratchpad scripts are lost).  This script must reproduce
//   (a) the recorded (a,b,c) of member (s,t,v)=(2,3,1)  [cert.m record]
//   (b) all 10 rows of data/claude_prod_09_88_stage1_members.txt
//       (J2 = [2,2,2,4], order-4 point doubling to class {0,a}, J1 = [4,4])
//   (c) all 3 rows of data/claude_prod_09_88_double_stage1.txt
//       (J2 = [2,2,4,4], order-4 points over BOTH {0,a} and {c,oo}, J1=[4,4])
//   (d) the (s,t,v)=(3,1/2,1) spot check J1 = [2,4,4]  [family.m record]
// Any mismatch is a transcription bug -- do not build on the defs file until
// FAMILY_REBUILD_VALIDATED prints.
//
// Run from repo root:
//   nohup magma -b code/claude_prod_09_88_family.m > results/claude_prod_09_88_family_rebuild.log 2>&1 &

SetColumns(0);
SetSeed(1);
SetMemoryLimit(3*10^9);

load "code/claude_prod_09_88_defs.m";

Q := Q_88; P := P88;

////////////////////////////////////////////////////////////////////////
// (a) formula gate: member (s,t,v) = (2,3,1)
////////////////////////////////////////////////////////////////////////
_, _, a0, b0, c0 := Lambda334(2, 3, 1);
printf "GATE_A (2,3,1): a=%o b=%o c=%o\n", a0, b0, c0;
assert a0 eq 17/2 and b0 eq 670761/69169 and c0 eq 9;
printf "GATE_A_PASS\n";

////////////////////////////////////////////////////////////////////////
// (b) stage-1 members: J2=[2,2,2,4], 4-pt over {0,a}, J1=[4,4]
////////////////////////////////////////////////////////////////////////
stage1 := [
    [2, 1, 1], [2, 1, -1], [2, 1, 2], [2, 3, 1], [3, 1, 1], [3, 2, -1],
    [1/2, 2, 1], [2/3, 1, 1], [2, 5, 1/2], [3, 1/3, 1]
];
for row in stage1 do
    m := row[1]; n := row[2]; v := row[3];
    s, t := StageOneST(m, n);
    h, g1, a, b, c := Lambda334(s, t, v);
    hi := IntSextic(h);
    J2 := Jacobian(HyperellipticCurve(hi));
    T0a := J2![x^2 - a*x, P!0];
    ok0a, inv2 := HasOrder4Over(J2, T0a);
    inv1 := ExactTorsion(g1);
    printf "STAGE1 m=%o n=%o v=%o J2=%o over0a=%o J1=%o\n",
           m, n, v, inv2, ok0a, inv1;
    assert inv2 eq [2,2,2,4] and ok0a and inv1 eq [4,4];
end for;
printf "GATE_B_PASS (10/10 stage-1 rows reproduced)\n";

////////////////////////////////////////////////////////////////////////
// (c) double-stage-1 members at base (m,n)=(3,1/3)
////////////////////////////////////////////////////////////////////////
dvals := [-729/17500, 26244/7975, 729/38425];
s, t := StageOneST(3, 1/3);
assert s eq 265/54 and t eq 5/3;
for v in dvals do
    h, g1, a, b, c := Lambda334(s, t, v);
    hi := IntSextic(h);
    J2 := Jacobian(HyperellipticCurve(hi));
    T0a := J2![x^2 - a*x, P!0];
    Tcoo := J2![x - c, P!0];
    ok0a, inv2 := HasOrder4Over(J2, T0a);
    okcoo, _ := HasOrder4Over(J2, Tcoo);
    inv1 := ExactTorsion(g1);
    printf "DOUBLE m=3 n=1/3 v=%o J2=%o over0a=%o overcoo=%o J1=%o\n",
           v, inv2, ok0a, okcoo, inv1;
    assert inv2 eq [2,2,4,4] and ok0a and okcoo and inv1 eq [4,4];
end for;
printf "GATE_C_PASS (3/3 double-stage-1 rows reproduced)\n";

////////////////////////////////////////////////////////////////////////
// (d) spot check: (s,t,v) = (3,1/2,1) has J1 = [2,4,4]
////////////////////////////////////////////////////////////////////////
_, g1d := Lambda334(3, 1/2, 1);
invd := ExactTorsion(g1d);
printf "GATE_D (3,1/2,1): J1=%o\n", invd;
assert invd eq [2,4,4];
printf "GATE_D_PASS\n";

printf "FAMILY_REBUILD_VALIDATED\n";
quit;
