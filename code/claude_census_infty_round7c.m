// Antidiagonal [2,2,10] slice: t3,t4 = +-s (common q), condition
// R = (q1-q) lam1 Cq(lam2) - (q2-q) lam2 Cq(lam1) = 0, kappa eliminated.
// Fix t2; study the (t1,s)-curve.  Known point (t1,t2,s) = (2,3,6).
SetColumns(0);
SetMemoryLimit(6*10^9);
Q := Rationals();

function WPq(t)
    return t^2/(t^2-1), t/(t^2-1);
end function;

K<t1,s> := RationalFunctionField(Q, 2);
for t2v in [3, 2, 5] do
    q1, m1 := WPq(t1);
    lam1 := (q1-1)*(-(2*q1-1)+2*m1)/q1;
    q2, m2 := WPq(K!t2v);
    lam2 := (q2-1)*(-(2*q2-1)+2*m2)/q2;
    q, m := WPq(s);
    Cq := func< xx | xx^2*q^2 + 2*xx*q*(2*q-1)*(q-1) + (q-1)^2 >;
    R := (q1-q)*lam1*Cq(lam2) - (q2-q)*lam2*Cq(lam1);
    NR<a,b> := PolynomialRing(Q, 2);
    num := NR!Numerator(R);
    printf "t2=%o: deg %o; factors:\n", t2v, TotalDegree(num);
    for tt in Factorization(num) do
        F := tt[1];
        d := TotalDegree(F);
        if d le 1 then printf "  lin: %o\n", F; continue; end if;
        A2 := AffineSpace(NR);
        Cu := Curve(A2, F);
        printf "  FACTOR deg %o genus %o : %o\n", d, Genus(Cu), F;
    end for;
end for;
printf "ROUND7C_DONE\n";
quit;
