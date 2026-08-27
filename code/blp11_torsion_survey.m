// Exact torsion survey of the 19 Bernard-Leprevost-Pohst order-11 curves
// (Exp. Math. 18 (2009) 65-70): C_{a,b,c,d}: y^2 = R^2 - 4c^2 S^2,
// R = x^3 - x^2 + a x + b, S = x^2 + d.  D_inf = (+inf) - (-inf) has order
// 11 on each.  Question (Filip): does any have 2-rank 2, i.e. torsion
// containing [2,22]?  The paper says nothing about 2-torsion; several
// curves may lie outside the LMFDB/alpha discriminant range, so this is
// computed directly.  Row C5 = X0(23) (RM by sqrt(5), the modular case).
// Row tC6 as printed duplicates tC3 (c = 15/9 = 5/3, same a,b,d; suspected
// typo in the paper) -- included and flagged.
//
// Run: magma -b code/blp11_torsion_survey.m > results/blp11_torsion_survey.log

SetColumns(0);
SetSeed(1);
SetMemoryLimit(3*10^9);

Q := Rationals();
P<x> := PolynomialRing(Q);

rows := [
    // <label, a, b, c, d>
    <"C1",  -101/48,   -61/48,     1/4,   -5/12>,
    <"C2",  473/147,   -4013/343,  6/7,   207/49>,
    <"C3",  8/49,      -134/49,    3/7,   47/49>,
    <"C4",  1159/81,   261607/2187, 40/9, 13/27>,
    <"C5",  -1/13,     -191/2197,  8/13,  15/169>,   // = X0(23), RM sqrt5
    <"C6",  -28/169,   103/2197,   3/13,  -4/169>,
    <"C7",  594/1805,  13348/34295, 8/19, -64/361>,
    <"C8",  208/867,   1338/4913,  5/17,  -39/289>,
    <"C9",  415/1089,  -2207/1089, 8/33,  119/121>,
    <"C10", 4989/2500, -13599/12500, 27/50, -81/250>,
    <"tC1", -3,        59,         4,     -7>,
    <"tC2", -163/1215, -367/3645,  2/3,   13/243>,
    <"tC3", -13/18,    71/6,       5/3,   -13/3>,
    <"tC4", -2287/27,  -1171/9,    10/3,  -323/3>,
    <"tC5", 121/147,   -141/343,   2/7,   15/49>,
    <"tC6", -13/18,    71/6,       15/9,  -13/3>,    // as printed; = tC3?
    <"tC7", -1494/847, 19480/9317, 2/11,  -256/121>,
    <"tC8", 125/121,   -223/1331,  6/11,  29/121>,
    <"tC9", 187/361,   -649/6859,  6/19,  23/361>
];

for row in rows do
    lbl := row[1]; a := row[2]; b := row[3]; c := row[4]; d := row[5];
    R := x^3 - x^2 + a*x + b;
    S := x^2 + d;
    f := R^2 - 4*c^2*S^2;
    if Degree(f) ne 6 or Discriminant(f) eq 0 then
        printf "ROW %o DEGENERATE\n", lbl;
        continue;
    end if;
    dtype := {* Degree(fp[1])^^fp[2] : fp in Factorization(f) *};
    C0 := HyperellipticCurve(f);
    Cc := C0;
    try
        Cm := ReducedMinimalWeierstrassModel(C0);
        Cc := SimplifiedModel(Cm);
    catch e ; end try;
    T := TorsionSubgroup(Jacobian(Cc));
    invs := Invariants(T);
    printf "ROW %o a=%o b=%o c=%o d=%o factortype=%o TORSION=%o\n",
           lbl, a, b, c, d, dtype, invs;
    if #invs ge 2 and invs[#invs] mod 11 eq 0 and invs[#invs-1] mod 2 eq 0 then
        printf "HIT22 %o TORSION=%o f=%o\n", lbl, invs, f;
    end if;
end for;
printf "BLP11_SURVEY_DONE\n";
quit;
