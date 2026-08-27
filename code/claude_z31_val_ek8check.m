// Validator task 4 (B4): independent checks of the EK disc-8 chart lane.
// (1) (r,s) = (6000/96721, 3557520/96721) -> Igusa-Clebsch via the archived
//     Elkies-Kumar formulas (data/claude_z31_ek8_source/igusa8.txt, transcribed
//     here BY THE VALIDATOR from the archived file) must match the reduced
//     witness model y^2 = -240x^6-120x^5-839x^4-310x^3-1103x^2-240x-504
//     weighted-projectively (weights 1,2,3,5).
// (2) TorsionSubgroup of that model = Z/31.
// (3) 6 random OK rows of refvec_ek8_P103/P1019.tsv recomputed from scratch:
//     EK IC mod p -> Mestre curve -> L-polynomial -> {chi(1),chi(-1)} multiset.
SetColumns(0);
SetMemoryLimit(8*10^9);

Q := Rationals();
r := Q!6000/96721; s := Q!3557520/96721;
A1 := 2*r*s^2;
A  := -(9*r*s + 4*r^2 + 4*r + 1)/3;
B1 := r*s^2*(3*s + 8*r - 2)/3;
B  := -(54*r^2*s + 81*r*s - 16*r^3 - 24*r^2 - 12*r - 2)/27;
B2 := r^2;
IC := [ -24*B1/A1, -12*A, 96*(A/A1)*B1 - 36*B, -4*A1*B2 ];
printf "EK_IC %o %o %o %o\n", IC[1], IC[2], IC[3], IC[4];

bad := r*s*(16*r*s^2+32*r^2*s-40*r*s-s+16*r^3+24*r^2+12*r+2)
          *(27*r*s^2+36*r^2*s+18*r*s-72*s+16*r^3+48*r^2+48*r+16);
printf "BADLOCUS_ZERO %o\n", bad eq 0;

R<x> := PolynomialRing(Q);
fw := -240*x^6 - 120*x^5 - 839*x^4 - 310*x^3 - 1103*x^2 - 240*x - 504;
Cw := HyperellipticCurve(fw);
ICw := IgusaClebschInvariants(Cw);
printf "MODEL_IC %o %o %o %o\n", ICw[1], ICw[2], ICw[3], ICw[4];
if IC[1] ne 0 then
    lam := ICw[1]/IC[1];
    okic := (ICw[2] eq lam^2*IC[2]) and (ICw[3] eq lam^3*IC[3]) and (ICw[4] eq lam^5*IC[4]);
    printf "IC_WEIGHTED_MATCH %o lambda=%o\n", okic, lam;
else
    print "IC_WEIGHTED_MATCH SKIP I2=0";
end if;

J := Jacobian(Cw);
T := TorsionSubgroup(J);
printf "TORSION_INVARIANTS %o ORDER %o\n", Invariants(T), #T;

// (3) independent refvec row recomputation
rows := [
<103,46,8,10524,10524>,
<103,32,102,7567,14639>,
<103,29,61,10771,10771>,
<1019,486,950,1076722,1003282>,
<1019,238,889,1047118,1030798>,
<1019,543,371,1037706,1037706>
];
for row in rows do
    p := row[1]; Fp := GF(p);
    rr := Fp!row[2]; ss := Fp!row[3];
    a1 := 2*rr*ss^2;
    aa := -(9*rr*ss + 4*rr^2 + 4*rr + 1)/3;
    b1 := rr*ss^2*(3*ss + 8*rr - 2)/3;
    bb := -(54*rr^2*ss + 81*rr*ss - 16*rr^3 - 24*rr^2 - 12*rr - 2)/27;
    b2 := rr^2;
    ICp := [ -24*b1/a1, -12*aa, 96*(aa/a1)*b1 - 36*bb, -4*a1*b2 ];
    okrow := false;
    try
        Cp := HyperellipticCurveFromIgusaClebsch(ICp);
        Lp := LPolynomial(Cp);
        c1 := Evaluate(Lp, 1); cm1 := Evaluate(Lp, -1);
        okrow := {* c1, cm1 *} eq {* row[4], row[5] *};
        printf "REFROW p=%o r0=%o s0=%o chi={%o,%o} want={%o,%o} %o\n",
            p, row[2], row[3], c1, cm1, row[4], row[5], okrow select "PASS" else "FAIL";
    catch e
        printf "REFROW p=%o r0=%o s0=%o ERROR %o\n", p, row[2], row[3], e`Object;
    end try;
end for;
print "ALLDONE";
quit;
