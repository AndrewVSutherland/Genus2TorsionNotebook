//////////////////////////////////////////////////////////////////////
//  Agent Z/48 probe: simple A(16) plus rational 3-torsion.
//
//  If J(Q) contains a point P of order 16 and a nonzero rational
//  3-torsion point Q, then P+Q has order 48.  Thus, for an already
//  certified A(16) tuple, the only additional condition is
//  J(Q)[3] != 0.
//
//  This script applies the standard good-reduction filter:
//      J(Q)[3] != 0  =>  3 divides #J(F_p)
//  for every good prime p != 3.  A single good prime with
//  #J(F_p) not divisible by 3 rules out the Z/48 route for that curve.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned prime_bound then
    prime_bound := 97;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;

if not assigned do_exact_torsion then
    do_exact_torsion := true;
elif Type(do_exact_torsion) eq MonStgElt then
    do_exact_torsion := do_exact_torsion in {"true", "True", "1", "yes", "Yes"};
end if;

Qq := Rationals();
Z := Integers();
Px<x> := PolynomialRing(Qq);
Pz<xz> := PolynomialRing(Z);

function NormalizeInvariants(inv)
    return Sort([Z!n : n in inv | Z!n ne 1]);
end function;

function ProductOfInvariants(inv)
    ninv := NormalizeInvariants(inv);
    if #ninv eq 0 then
        return 1;
    end if;
    return &*ninv;
end function;

function FactorDegrees(f)
    return [<Degree(item[1]), item[2]> : item in Factorization(f)];
end function;

function SquareIntegralScale(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return L^2;
end function;

function SquareIntegralPolynomial(f)
    return Pz!(SquareIntegralScale(f)*f);
end function;

function A8Data(rv, pv, tv)
    r := Qq!rv;
    p := Qq!pv;
    t := Qq!tv;

    e := t^2 - 2*p*t/r;
    s := p - r^2;
    d := e + 2*p - r^2;
    lambda := r/t;
    u := p + r*t - 2*r;
    v := e + r^2 - r*p - r^2*t + 3*p*t - r*t^2;

    a := r^2 - lambda;
    b := 2*r*p - 2*lambda*(p + r*t) + 2*r*lambda;
    c := p^2 + 2*p*r^2 - r^4 - r^3*t - r*p^2/t
         - lambda*(r^2 + e)
         + 2*lambda*(r*p + r^2*t - 3*p*t + r*t^2);

    Qpoly := x^2 + d;
    q := a*x^2 + b*x + c;
    L := r*x + s;
    g8 := x^2 + u*x + v;
    f := q*(Qpoly^2 + q);
    ell8base := -(q + Qpoly*L);
    return f, g8, ell8base;
end function;

function A8CandidateData(rv, tv, pv, muv, yv, Nv, zv)
    r := Qq!rv;
    t := Qq!tv;
    p := Qq!pv;
    mu := Qq!muv;
    y := Qq!yv;
    N := Qq!Nv;
    z := Qq!zv;

    f, g8, ellBase := A8Data(r, p, t);
    ell := ellBase + g8*(mu*x + N);
    W := x^2 + y*x + z;
    scale := mu^2 - 2*r*mu + r/t;

    squareOK := false;
    if (ell^2 - f) mod g8 eq 0 then
        S := ExactQuotient(ell^2 - f, g8);
        squareOK := S eq scale*W^2;
    end if;

    return f, g8, ellBase, ell, W, squareOK, scale;
end function;

function VerifyA16Tuple(rv, tv, pv, muv, yv, Nv, zv)
    f, g8, ellBase, ell, W, squareOK, scale := A8CandidateData(
        rv, tv, pv, muv, yv, Nv, zv);
    if not squareOK or Degree(f) lt 5 or Discriminant(f) eq 0 then
        return false, 0, squareOK, scale, f;
    end if;

    try
        C := HyperellipticCurve(f);
        J := Jacobian(C);
        ZJ := J!0;
        D8 := J![g8, -ellBase mod g8];
        D16minus := J![W, -ell mod W];
        D16plus := J![W, ell mod W];

        if 8*D8 ne ZJ or 4*D8 eq ZJ then
            return false, 0, squareOK, scale, f;
        end if;
        if 2*D16minus eq D8 and 16*D16minus eq ZJ and 8*D16minus ne ZJ then
            return true, -1, squareOK, scale, f;
        elif 2*D16plus eq D8 and 16*D16plus eq ZJ and 8*D16plus ne ZJ then
            return true, 1, squareOK, scale, f;
        end if;
    catch err
        return false, 0, squareOK, scale, f;
    end try;

    return false, 0, squareOK, scale, f;
end function;

function GoodReductionPolynomial(f, p)
    F := GF(p);
    PF<xp> := PolynomialRing(F);
    try
        fp := PF![F!Coefficient(f, i) : i in [0..Degree(f)]];
    catch err
        return false, PF!0;
    end try;
    if Degree(fp) lt 5 or Discriminant(fp) eq 0 then
        return false, fp;
    end if;
    return true, fp;
end function;

function PointCountGCD3Filter(f, prime_bound)
    primes := [p : p in PrimesUpTo(prime_bound) | p notin {2,3}];
    goodData := [];
    gcdN := 0;
    firstKill := <0, 0>;

    for p in primes do
        good, fp := GoodReductionPolynomial(f, p);
        if not good then
            continue;
        end if;

        C := HyperellipticCurve(fp);
        Np := Z!Evaluate(LPolynomial(C), 1);
        Append(~goodData, <p, Np, Np mod 3>);
        if gcdN eq 0 then
            gcdN := Np;
        else
            gcdN := GCD(gcdN, Np);
        end if;

        if firstKill[1] eq 0 and Np mod 3 ne 0 then
            firstKill := <p, Np>;
        end if;
    end for;

    return goodData, gcdN, firstKill;
end function;

TupleFormat := recformat<label, r, t, p, mu, y, N, z>;

// These are the certified A(8)->A(16) tuples recorded in
// a8_order16_cover_status.md and torsion_goal_log.md.  The first two are
// distinct halves on the same underlying curve; keeping both checks
// tuple-level reconstruction.
a16Tuples := [
    rec<TupleFormat | label := "rt_3_1over3_p2_mu9",
        r := Qq!3, t := Qq!1/3, p := Qq!2,
        mu := Qq!9, y := Qq!-4, N := Qq!-12, z := Qq!8/3>,
    rec<TupleFormat | label := "rt_3_1over3_p2_mu27over11",
        r := Qq!3, t := Qq!1/3, p := Qq!2,
        mu := Qq!27/11, y := Qq!-8, N := Qq!-60/11, z := Qq!32/3>,
    rec<TupleFormat | label := "rt_3_1over3_p34over9",
        r := Qq!3, t := Qq!1/3, p := Qq!34/9,
        mu := Qq!9/17, y := Qq!-4/21, N := Qq!-12/17,
        z := Qq!-1544/567>,
    rec<TupleFormat | label := "rt_minus1_1over3_p1over3",
        r := Qq!-1, t := Qq!1/3, p := Qq!1/3,
        mu := Qq!3/7, y := Qq!4/3, N := Qq!8/7,
        z := Qq!8/3>,
    rec<TupleFormat | label := "rt_minus1_1over2_pminus35over6",
        r := Qq!-1, t := Qq!1/2, p := Qq!-35/6,
        mu := Qq!1, y := Qq!-7, N := Qq!-13/2,
        z := Qq!77/4>,
    rec<TupleFormat | label := "rt_3_1over2_p17over6",
        r := Qq!3, t := Qq!1/2, p := Qq!17/6,
        mu := Qq!1, y := Qq!-1/3, N := Qq!-7/6,
        z := Qq!-227/36>
];

print "AGENT_Z48_A16_PLUS3_PROBE";
print "prime_bound", prime_bound, "do_exact_torsion", do_exact_torsion;
print "candidate_count", #a16Tuples;

seenModels := {};

for tup in a16Tuples do
    print "";
    print "CANDIDATE", tup`label;
    printf "  tuple r=%o t=%o p=%o mu=%o y=%o N=%o z=%o\n",
        tup`r, tup`t, tup`p, tup`mu, tup`y, tup`N, tup`z;

    a16OK, sign, squareOK, scale, f := VerifyA16Tuple(
        tup`r, tup`t, tup`p, tup`mu, tup`y, tup`N, tup`z);
    fInt := SquareIntegralPolynomial(f);
    modelKey := Sprint([Coefficient(fInt, i) : i in [0..Degree(fInt)]]);
    duplicateModel := modelKey in seenModels;
    Include(~seenModels, modelKey);

    print "  square_relation", squareOK, "scale", scale;
    print "  a16_verified", a16OK, "sign", sign;
    print "  factor_degrees", FactorDegrees(f);
    print "  square_integral_model", fInt;
    print "  duplicate_integral_model", duplicateModel;

    if do_exact_torsion then
        try
            Cq := HyperellipticCurve(fInt);
            Jq := Jacobian(Cq);
            A, mp := TorsionSubgroup(Jq);
            inv := NormalizeInvariants(Invariants(A));
            print "  exact_torsion", inv, "order", ProductOfInvariants(inv);
        catch err
            print "  exact_torsion_failed", err`Object;
        end try;
    end if;

    goodData, gcdN, firstKill := PointCountGCD3Filter(f, prime_bound);
    print "  good_primes_used", #goodData;
    print "  point_count_gcd", gcdN, "gcd_mod_3", gcdN mod 3;
    if firstKill[1] ne 0 then
        print "  THREE_TORSION_KILLED", "p", firstKill[1],
              "J_order", firstKill[2], "mod3", firstKill[2] mod 3;
    else
        print "  THREE_TORSION_SURVIVES_FILTER";
    end if;
    print "  point_counts", goodData;
end for;

print "";
print "DONE";
