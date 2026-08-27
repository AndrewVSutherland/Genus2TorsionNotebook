//////////////////////////////////////////////////////////////////////
//  claude_ov_2216_delta.m      -- shared x-T (delta) module
//
//  Odd quintic model  y^2 = f(X),  deg f = 5, factor type [1,1,1,2]:
//      f = c * (X-r1)(X-r2)(X-r3) * qtm(X),   qtm monic irreducible.
//  L = Q[X]/(f) = Q x Q x Q x K,  K = Q[X]/(qtm).
//
//  CONVENTION.  We compute a representative delta' of the x-T image in
//  L*/(L*)^2 Q*  (i.e. only up to a common rational scalar), by the
//  classical formulas
//      coprime class D :  delta'_j = u_D(a_j)
//      2-torsion T_S   :  delta'_j = uS(a_j)            (a_j notin S)
//                         delta'_j = (f/uS)(a_j)        (a_j in S)
//  with S an even subset of the FINITE Weierstrass points (use the
//  complementary subset when the class is [(r_i,0) - infty]).
//
//  The honest image in L*/(L*)^2 is then  delta = lambda0 * delta'  where
//      lambda0 = squarefree part of N_{L/Q}(delta'),
//  because N(delta) must be a square and N(lambda*delta') = lambda^5 N(delta').
//
//  Then, x-T being injective on J(Q)/2J(Q) for odd degree,
//      D in 2 J(Q)  <==>  lambda0*delta'_j square in Q (j=1,2,3)
//                          and lambda0*delta'_K square in K.
//////////////////////////////////////////////////////////////////////

Qq := Rationals();
PP<XX> := PolynomialRing(Qq);

function CO_SqClass(a)
    n := Numerator(a); d := Denominator(a);
    m := n*d; s := Sign(m); m := AbsoluteValue(m); r := 1;
    for pe in Factorization(m) do
        if IsOdd(pe[2]) then r *:= pe[1]; end if;
    end for;
    return s*r;
end function;

function CO_EvalK(g, KK, th)
    val := KK!0; cs := Coefficients(g); pw := KK!1;
    for i in [1..#cs] do val +:= (KK!cs[i])*pw; pw *:= th; end for;
    return val;
end function;

// delta' for a class with u_D coprime to f
function CO_DeltaPrimeCoprime(uD, roots, KK, th)
    d := [];
    for r in roots do
        val := Evaluate(uD, r);
        if val eq 0 then return false, _; end if;
        Append(~d, val);
    end for;
    valK := CO_EvalK(uD, KK, th);
    if valK eq 0 then return false, _; end if;
    return true, <d[1], d[2], d[3], valK>;
end function;

// delta' for the 2-torsion class of an even subset S of the finite
// Weierstrass points, encoded by the monic divisor uS | f.
//
// f has ODD degree 5, so infinity is the sixth Weierstrass point and it
// is NOT in S (S consists of finite points).  In the projective form
// F(x,z) = z^6 f(x/z) the "difference" of a finite root a and the root
// [1:0] at infinity is  a*0 - 1*1 = -1.  Hence the factor for a_j in S
// picks up exactly one extra  (-1)  from the infinite root:
//      a_j notin S :  delta'_j = uS(a_j)
//      a_j in S    :  delta'_j = -(f/uS)(a_j)
// (the leading coefficient of f/uS is a GLOBAL constant and drops mod Q*).
function CO_DeltaPrimeTwo(uS, f, roots, KK, th)
    cofac := f div uS;
    if uS*cofac ne f then return false, _; end if;
    d := [];
    for r in roots do
        if Evaluate(uS, r) eq 0 then val := -Evaluate(cofac, r);
        else val := Evaluate(uS, r); end if;
        if val eq 0 then return false, _; end if;
        Append(~d, val);
    end for;
    vS := CO_EvalK(uS, KK, th);
    if vS eq 0 then valK := -CO_EvalK(cofac, KK, th); else valK := vS; end if;
    if valK eq 0 then return false, _; end if;
    return true, <d[1], d[2], d[3], valK>;
end function;

function CO_MulPrime(a, b)
    return <a[1]*b[1], a[2]*b[2], a[3]*b[3], a[4]*b[4]>;
end function;

// normalize a delta' to the honest delta and test triviality
function CO_IsTrivial(dp)
    nrm := dp[1]*dp[2]*dp[3]*Norm(dp[4]);
    if nrm eq 0 then return false; end if;
    lam := CO_SqClass(nrm);
    if not IsSquare(lam*dp[1]) then return false; end if;
    if not IsSquare(lam*dp[2]) then return false; end if;
    if not IsSquare(lam*dp[3]) then return false; end if;
    return IsSquare(lam*dp[4]);
end function;

// compact printable signature of delta' (square classes + lambda0)
function CO_Signature(dp)
    nrm := dp[1]*dp[2]*dp[3]*Norm(dp[4]);
    lam := CO_SqClass(nrm);
    return [ CO_SqClass(lam*dp[1]), CO_SqClass(lam*dp[2]), CO_SqClass(lam*dp[3]) ],
           IsSquare(lam*dp[4]), lam;
end function;
