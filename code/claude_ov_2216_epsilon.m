//////////////////////////////////////////////////////////////////////
//  claude_ov_2216_epsilon.m
//
//  Pin down the normalisation of the x-T map for ODD-degree genus-2
//  models.  Write  Delta_n(D)_i = u_D(a_i)  (the "naive" value, defined
//  when gcd(u_D,f)=1) and posit
//
//        delta(D) = eps^{deg u_D} * Delta_n(D)   in L*/(L*)^2 .
//
//  Norm computation:  N(Delta_n(D)) = prod_i u_D(a_i)
//                   = Res(f,u_D)/c^{m} = (-1)^m * (prod_k y_k)^2 / c^m
//                   ~ (-c)^m ,
//  and the honest x-T image has square norm, so N(eps) ~ -c.  If eps is
//  a CONSTANT class it must therefore be eps = -c = -lc(f).
//
//  TEST: find random odd-degree curves carrying a class D with
//  deg u_D = 1 that IS divisible by 2 in J(Q); then delta(D) is trivial,
//  so eps = Delta_n(D).  Compare with -lc(f).
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(8*10^9);
Qq := Rationals();
P<X> := PolynomialRing(Qq);

function SqClass(a)
    n := Numerator(a); d := Denominator(a);
    m := n*d; s := Sign(m); m := AbsoluteValue(m); r := 1;
    for pe in Factorization(m) do
        if IsOdd(pe[2]) then r *:= pe[1]; end if;
    end for;
    return s*r;
end function;

SetSeed(11223344);
ntest := 0; nagree := 0; ndis := 0;

for iter in [1..4000] do
    if ntest ge 25 then break; end if;
    // random monic-ish quintic with 5 distinct rational roots so that all
    // components of L are Q and square classes are easy to read off.
    rts := [];
    while #rts lt 5 do
        r := Random(-9,9);
        if r notin rts then Append(~rts, r); end if;
    end while;
    c := Random([-6,-3,-2,-1,1,2,3,5]);
    f := c*&*[ X - r : r in rts ];
    if Discriminant(f) eq 0 then continue; end if;
    C := HyperellipticCurve(f);
    J := Jacobian(C);

    // build a few classes and double them, looking for deg u = 1
    pts := Points(C : Bound := 40);
    cls := [];
    for pp in pts do
        try D := J!(Divisor(pp) - Divisor(PointsAtInfinity(C)[1]));
             Append(~cls, D);
        catch e ; end try;
    end for;
    if #cls lt 2 then continue; end if;
    cand := [];
    for i in [1..#cls] do
        for j in [i..#cls] do
            Append(~cand, cls[i]+cls[j]);
        end for;
    end for;
    for D in cand do
        E := 2*D;
        if E eq J!0 then continue; end if;
        uE := E[1];
        if Degree(uE) ne 1 then continue; end if;
        if GCD(uE, f) ne 1 then continue; end if;
        // delta(E) must be trivial; so eps = Delta_n(E)
        vals := [ SqClass(Evaluate(uE, Qq!r)) : r in rts ];
        epsguess := SqClass(-c);
        ok := &and[ w eq epsguess : w in vals ];
        ntest +:= 1;
        if ok then nagree +:= 1;
        else ndis +:= 1;
            printf "MISMATCH f=%o  lc=%o  -lc_sqclass=%o  Delta_n(2D)=%o\n",
                   f, c, epsguess, vals;
        end if;
        break;
    end for;
end for;

printf "eps-test: cases=%o  eps==-lc(f): agree=%o disagree=%o\n", ntest, nagree, ndis;
print "SEARCH_DONE";
quit;
