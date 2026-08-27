// [2,10]-locus: collect (a,b) with quartic [2,2]-split, wide box; fit curve.
SetColumns(0);
SetMemoryLimit(6*10^9);
Q := Rationals();
P<x> := PolynomialRing(Q);
pts := [];
for aN in [-40..40], aD in [1..4], bN in [-40..40], bD in [1..4] do
    if GCD(Abs(aN),aD) ne 1 or GCD(Abs(bN),bD) ne 1 then continue; end if;
    a := Q!aN/aD; b := Q!bN/bD;
    if 1+a+b eq 0 or b eq 0 then continue; end if;
    f := (1+a*x+b*x^2)^2 - (1+a+b)^2*x^5;
    if Degree(f) lt 5 or not IsSquarefree(f) then continue; end if;
    quart := f div (x-1);
    if Degree(quart) ne 4 then continue; end if;
    degs := Sort([Degree(t[1]) : t in Factorization(quart)]);
    if degs eq [2,2] then
        Append(~pts, [a,b]);
        printf "PT %o %o\n", a, b;
    end if;
end for;
printf "NPTS %o\n", #pts;
quit;
