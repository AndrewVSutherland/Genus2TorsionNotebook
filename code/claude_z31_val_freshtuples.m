// Validator task 2: fresh random chart-U' tuples -> exact order of D = [div(u,v)-D_inf]
// Truth source for cross-check against B2's Kummer ladder (claude_z31_kummer_sieve selftest).
// Tuples generated independently by the validator (python random, seed 20260730).
SetColumns(0);
SetMemoryLimit(3*10^9);

tuples := [
[101,89,76,22,32,68,21,4,93,72],
[211,102,143,151,15,44,138,84,145,151],
[499,256,239,418,182,36,392,265,126,387],
[1009,302,971,368,201,231,413,723,286,250],
[2003,141,1422,133,977,781,804,471,855,1984],
[101,50,86,90,62,68,84,56,5,97],
[211,88,25,179,31,32,69,62,66,90],
[499,299,250,62,94,139,1,181,307,324],
[1009,739,689,646,73,496,378,811,843,679],
[2003,1438,887,171,725,642,1490,941,791,1694],
[101,42,28,74,79,67,51,5,97,77],
[211,200,155,205,176,68,71,100,197,51],
[499,429,203,256,204,115,307,491,379,454],
[1009,970,797,902,290,447,829,640,644,953],
[2003,513,792,189,1529,285,1628,16,953,1302],
[101,12,29,49,64,14,98,14,61,88],
[211,35,155,151,121,185,0,68,113,154],
[499,140,335,54,195,226,232,164,417,403],
[1009,197,806,711,325,284,441,237,897,957],
[2003,1574,2001,925,770,1244,1765,1719,616,1501]
];

for t in tuples do
    P := t[1];
    Fp := GF(P); R<x> := PolynomialRing(Fp);
    u1 := Fp!t[2]; u0 := Fp!t[3]; v1 := Fp!t[4]; v0 := Fp!t[5];
    w4 := Fp!t[6]; w3 := Fp!t[7]; w2 := Fp!t[8]; w1 := Fp!t[9]; w0 := Fp!t[10];
    u := x^2 + u1*x + u0;
    v := v1*x + v0;
    w := w4*x^4 + w3*x^3 + w2*x^2 + w1*x + w0;
    f := v^2 + u*w;
    if Degree(f) ne 6 or Discriminant(f) eq 0 then
        printf "SKIP degenf P=%o t=%o\n", P, t; continue;
    end if;
    if u1^2 - 4*u0 eq 0 then
        printf "SKIP den0 P=%o t=%o\n", P, t; continue;
    end if;
    ok := false; meth := 0;
    C := HyperellipticCurve(f); J := Jacobian(C);
    DD := J!0;
    try
        DD := J![u, v]; ok := true; meth := 1;
    catch e
        ok := false;
    end try;
    if not ok then
        try
            DD := elt<J | u, v, 2>; ok := true; meth := 2;
        catch e
            ok := false;
        end try;
    end if;
    if not ok then
        printf "SKIP nodd P=%o t=%o\n", P, t; continue;
    end if;
    N := Order(DD);
    flag := (N mod 31 eq 0) select 1 else 0;
    printf "ROW %o %o %o %o %o %o %o %o %o %o %o %o meth=%o\n",
        P, t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], t[10], N, flag, meth;
end for;
print "ALLDONE";
quit;
