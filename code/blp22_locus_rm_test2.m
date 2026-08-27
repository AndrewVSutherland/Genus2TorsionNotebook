// BLP [2,22] engine, Humbert-sqrt5 containment test, v2.
// Fixes of v1 (results/blp22_locus_rm_test.log):
//  (1) SAMPLING: the locus is codim 2 in (r,q,s,d); random (r,s,d) misses
//      it (~1/p).  Now: fix random (r,s), treat (q,d) as unknowns, solve
//      T0 = T1 = 0 via Res_q over F_p(d), root-find in d, then gcd in q.
//  (2) TWIST-AWARENESS: a mod-p reduction of an RM-sqrt5 surface whose RM
//      is defined only over F_{p^2} can fail D = 5 v^2 for pi while its
//      SQUARE satisfies it.  Signature now: RM5-compatible iff
//         D1 := a^2 - 4b + 8p in 5*squares   (untwisted)
//      or D2 := a2^2 - 4 b2 + 8 p^2 in 5*squares, where
//         a2 := a^2 - 2b,  b2 := b^2 - 2 p a^2 + 2 p^2   (pi^2 coeffs).
// Verdict: all/most locus points RM5-compatible => locus likely inside
// Humbert-5 (engine cannot give non-RM [2,22]); substantial fraction
// incompatible => locus NOT contained, keep hunting rational points.
//
// Run: magma -b code/blp22_locus_rm_test2.m > results/blp22_locus_rm_test2.log

SetColumns(0);
SetSeed(1);
SetMemoryLimit(3*10^9);

function DeriveT(a, b, c, d)
    K := Parent(a);
    P := PolynomialRing(K); x := P.1;
    R := x^3 - x^2 + a*x + b;
    S := x^2 + d;
    F := R^2 - 4*c^2*S^2;
    if Degree(F) ne 6 then return K!0, K!0, false; end if;
    fc := [Coefficient(F, i) : i in [0..6]];
    s := [K| 1];
    for k in [1..9] do
        Ftk := k le 6 select fc[6-k+1] else K!0;
        acc := Ftk;
        for i in [1..k-1] do
            acc -:= s[i+1]*s[k-i+1];
        end for;
        Append(~s, acc/2);
    end for;
    M := Matrix(K, 3, 3, [ [ s[3+m+1], s[4+m+1], s[5+m+1] ] : m in [1..3] ]);
    rhs := Vector(K, [ -s[6+m+1] : m in [1..3] ]);
    if Determinant(M) eq 0 then return K!0, K!0, false; end if;
    u := Solution(Transpose(M), rhs);
    U := x^3 + u[3]*x^2 + u[2]*x + u[1];
    uco := [u[1], u[2], u[3], K!1];
    V := P!0;
    for n in [0..6] do
        vn := K!0;
        for i in [Max(0, n-3)..3] do
            j := i + 3 - n;
            if j ge 0 and j le 9 then vn +:= uco[i+1]*s[j+1]; end if;
        end for;
        V +:= vn*x^n;
    end for;
    G := F*U^2 - V^2;
    if Degree(G) gt 5 then return K!0, K!0, false; end if;
    t := Coefficient(G, 2);
    return Coefficient(G, 0) - t*d, Coefficient(G, 1), true;
end function;

function Is5Sq(D)
    return D ge 0 and D mod 5 eq 0 and IsSquare(D div 5);
end function;

p := 499;
Fp := GF(p);
K2<qv, dv> := RationalFunctionField(Fp, 2);
P2 := PolynomialRing(Fp, 2);
PF<X> := PolynomialRing(Fp);
PZ := PolynomialRing(Integers());

ncompat := 0; nincompat := 0; nsplit := 0; npts := 0; trials := 0;
while npts lt 40 and trials lt 200 do
    trials +:= 1;
    r := Fp!Random(Fp); sv := Fp!Random(Fp);
    if sv eq r then continue; end if;
    rK := K2!r; sK := K2!sv;
    den := 2*(sK^2 + dv) - sK*(sK - rK);
    if den eq 0 then continue; end if;
    pf := (2*(rK-1)*(sK^2+dv) + (sK-rK)*(sK^2+qv)) / den;
    cf := (rK - 1 - pf)/2;
    aq := qv - rK*pf;
    bq := -rK*qv + 2*cf*dv;
    t0, t1, ok := DeriveT(aq, bq, cf, K2!dv);
    if not ok or t0 eq 0 or t1 eq 0 then continue; end if;
    N0 := P2!Numerator(t0); N1 := P2!Numerator(t1);
    if Degree(N0, 1) lt 1 or Degree(N1, 1) lt 1 then continue; end if;
    Rd := Resultant(N0, N1, P2.1);    // eliminate q -> poly in d
    if Rd eq 0 then continue; end if;
    Rdn := UnivariatePolynomial(Rd);
    for rtd in Roots(Rdn) do
        d0 := rtd[1];
        if sv^2 + d0 eq 0 then continue; end if;
        N0d := UnivariatePolynomial(Evaluate(N0, 2, d0));
        N1d := UnivariatePolynomial(Evaluate(N1, 2, d0));
        if N0d eq 0 or N1d eq 0 then continue; end if;
        gg := GCD(PF!N0d, PF!N1d);
        if Degree(gg) lt 1 then continue; end if;
        for rtq in Roots(gg) do
            q0 := rtq[1];
            den0 := 2*(sv^2 + d0) - sv*(sv - r);
            if den0 eq 0 then continue; end if;
            p0 := (2*(r-1)*(sv^2+d0) + (sv-r)*(sv^2+q0)) / den0;
            c0 := (r - 1 - p0)/2;
            if c0 eq 0 then continue; end if;
            g0 := (X - r)*(X^2 + p0*X + q0);
            f := g0*(g0 + 4*c0*(X^2 + d0));
            if Degree(f) ne 6 or Degree(GCD(f, Derivative(f))) gt 0 then continue; end if;
            C := HyperellipticCurve(f);
            chi := PZ!Reverse(Coefficients(LPolynomial(C)));
            npts +:= 1;
            if not IsIrreducible(chi) then
                nsplit +:= 1;
            else
                av := -Coefficient(chi, 3); bv := Coefficient(chi, 2);
                D1 := av^2 - 4*bv + 8*p;
                a2 := av^2 - 2*bv; b2 := bv^2 - 2*p*av^2 + 2*p^2;
                D2 := a2^2 - 4*b2 + 8*p^2;
                if Is5Sq(D1) or Is5Sq(D2) then
                    ncompat +:= 1;
                else
                    nincompat +:= 1;
                    printf "INCOMPAT r=%o s=%o d=%o q=%o D1=%o D2=%o\n",
                           r, sv, d0, q0, D1, D2;
                end if;
            end if;
            break;
        end for;
        if npts ge 40 then break; end if;
    end for;
end while;
printf "P=%o points=%o RM5_compatible=%o incompatible=%o split=%o\n",
       p, npts, ncompat, nincompat, nsplit;
printf "LOCUS_RM_TEST2_DONE\n";
quit;
