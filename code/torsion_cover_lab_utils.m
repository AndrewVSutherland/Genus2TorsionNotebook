//////////////////////////////////////////////////////////////////////
//  Lightweight shared helpers for torsion-cover experiments.
//
//  This file is intentionally small and additive.  Existing scripts in this
//  directory mostly carry their own helper functions; new scripts can load
//  this file when they only need standard height enumeration, factor-degree
//  summaries, and point-count gates.
//////////////////////////////////////////////////////////////////////

function TC_SumInts(xs)
    if #xs eq 0 then
        return 0;
    end if;
    return &+xs;
end function;

function TC_MakeMonic(g)
    if g eq 0 then
        return g;
    end if;
    return g/LeadingCoefficient(g);
end function;

function TC_HeightRationals(H)
    QQ := Rationals();
    vals := [QQ!0];
    for den in [1..H] do
        for num in [-H..H] do
            if GCD(num, den) eq 1 then
                Append(~vals, QQ!num/QQ!den);
            end if;
        end for;
    end for;
    return Sort(Setseq(Seqset(vals)));
end function;

function TC_FactorDegreeMults(g)
    if g eq 0 or Degree(g) lt 1 then
        return [];
    end if;
    return [<Degree(fe[1]), fe[2]> : fe in Factorization(g)];
end function;

function TC_ContainsPoint(seq, T)
    for S in seq do
        if S eq T then
            return true;
        end if;
    end for;
    return false;
end function;

function TC_GoodReductionPolynomial(f, ell)
    F := GF(ell);
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

function TC_PointCountGate(f, primes, minGood, gateMod)
    Z := Integers();
    used := [];
    gcdN := 0;
    good := 0;

    for ell in primes do
        ok, fp := TC_GoodReductionPolynomial(f, ell);
        if not ok then
            continue;
        end if;

        C := HyperellipticCurve(fp);
        Np := Z!#Jacobian(C);
        if gcdN eq 0 then
            gcdN := Np;
        else
            gcdN := GCD(gcdN, Np);
        end if;
        good +:= 1;
        Append(~used, <ell, Np, Np mod gateMod, gcdN>);

        if Np mod gateMod ne 0 then
            return false, "killed", ell, Np, good, gcdN, used;
        end if;
        if good ge minGood then
            return true, "passed", 0, 0, good, gcdN, used;
        end if;
    end for;

    return false, "insufficient_good", 0, 0, good, gcdN, used;
end function;

function TC_NormalizeTorsionInvariants(inv)
    Z := Integers();
    return Sort([Z!n : n in inv | Z!n ne 1]);
end function;
