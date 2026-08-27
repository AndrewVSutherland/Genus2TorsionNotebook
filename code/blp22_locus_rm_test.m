// BLP [2,22] engine: is the [1,1,2,2]+order-11 locus contained in the
// Humbert-sqrt5 (RM by sqrt 5) surface?
//
// Mod-p test: points of the locus over F_p are plentiful (q solves the
// gcd of the two closure numerators over F_p).  If the GLOBAL locus lies
// in Humbert-5, then for every mod-p locus point the Frobenius char poly
// chi = x^4 - a x^3 + b x^2 - p a x + p^2 must be the Q(sqrt5)-norm of a
// quadratic, which happens iff
//     D := a^2 - 4b + 8p  =  5 v^2   for some integer v
// (alpha = u + v(1+sqrt5)/2 gives Tr(alpha)^2 - 4N(alpha) = 5 v^2).
// A generic non-RM point has D essentially random in [0, 16p]: the
// probability of 5*square is O(1/sqrt p).  So:
//   fraction(D = 5 v^2) ~ 1  =>  locus inside Humbert-5 (engine can NEVER
//                                yield a non-RM [2,22]);
//   fraction well below 1   =>  locus NOT contained; keep hunting.
// Split/reducible-chi points are tallied separately.
//
// Run: magma -b code/blp22_locus_rm_test.m > results/blp22_locus_rm_test.log

SetColumns(0);
SetSeed(1);
SetMemoryLimit(3*10^9);

// ---- DeriveT (verbatim from code/blp22_engine_derive.m) ----
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

for p in [997, 1009] do
    Fp := GF(p);
    Kq<qv> := RationalFunctionField(Fp);
    Pq := PolynomialRing(Fp);
    PF<X> := PolynomialRing(Fp);
    n5 := 0; nother := 0; nsplit := 0; ntried := 0;
    trials := 0;
    while ntried lt 40 and trials lt 500 do
        trials +:= 1;
        r := Random(Fp); sv := Random(Fp); d := Random(Fp);
        if sv eq r or sv^2 + d eq 0 then continue; end if;
        den := 2*(sv^2 + d) - sv*(sv - r);
        if den eq 0 then continue; end if;
        pfun := (2*(r-1)*(sv^2+d) + (sv-r)*(sv^2+qv)) / den;
        cfun := (r - 1 - pfun)/2;
        aq := qv - r*pfun;
        bq := -r*qv + 2*cfun*d;
        t0, t1, ok := DeriveT(aq, bq, cfun, d);
        if not ok or (t0 eq 0 and t1 eq 0) then continue; end if;
        N0 := Pq!Numerator(t0); N1 := Pq!Numerator(t1);
        if N0 eq 0 or N1 eq 0 then continue; end if;
        g := GCD(N0, N1);
        if Degree(g) lt 1 then continue; end if;
        for rt in Roots(g) do
            q0 := rt[1];
            p0 := Evaluate(pfun, q0); c0 := Evaluate(cfun, q0);
            if c0 eq 0 then continue; end if;
            g0 := (X - r)*(X^2 + p0*X + q0);
            f := g0*(g0 + 4*c0*(X^2 + d));
            if Degree(f) ne 6 or Degree(GCD(f, Derivative(f))) gt 0 then continue; end if;
            C := HyperellipticCurve(f);
            chi := PolynomialRing(Integers())!Reverse(Coefficients(LPolynomial(C)));
            ntried +:= 1;
            if not IsIrreducible(chi) then nsplit +:= 1; continue; end if;
            av := -Coefficient(chi, 3); bv := Coefficient(chi, 2);
            D := av^2 - 4*bv + 8*p;
            is5sq := D ge 0 and D mod 5 eq 0 and IsSquare(D div 5);
            if is5sq then n5 +:= 1; else nother +:= 1;
                printf "NON5 p=%o r=%o s=%o d=%o q=%o D=%o\n", p, r, sv, d, q0, D;
            end if;
            break;   // one locus point per (r,s,d) triple
        end for;
    end while;
    printf "P=%o locus_points=%o D_in_5squares=%o other=%o split_chi=%o\n",
           p, ntried, n5, nother, nsplit;
end for;
printf "LOCUS_RM_TEST_DONE\n";
quit;
