// [2,2,10] probe: Elkies chart with FOUR rational Weierstrass points.
// det of interpolation matrix = 0; fix (t1,t2)=(2,3), study the curve in (t3,t4).
SetColumns(0);
SetMemoryLimit(6*10^9);
Q := Rationals();
K<t3,t4> := RationalFunctionField(Q, 2);

function WPdata(t)
    q := t^2/(t^2-1); m := t/(t^2-1);
    return (q-1)*(-(2*q-1)+2*m)/q, q;
end function;

t1 := K!2; t2 := K!3;
ts := [t1, t2, t3, t4];
lams := []; qs := [];
for t in ts do
    lam, q := WPdata(t);
    Append(~lams, lam); Append(~qs, q);
end for;
M := Matrix(K, 4, 4, [ [lams[i]^2, lams[i], 1, qs[i]*lams[i]] : i in [1..4] ]);
d := Determinant(M);
NR := PolynomialRing(Q, 2);
num := NR!Numerator(d);
printf "det numerator: total degree %o\n", TotalDegree(num);
fac := Factorization(num);
for tt in fac do
    F := tt[1];
    printf "FACTOR deg %o mult %o: %o\n", TotalDegree(F), tt[2], F;
end for;
printf "ROUND7_DONE\n";
quit;
