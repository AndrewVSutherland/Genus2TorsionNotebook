//////////////////////////////////////////////////////////////////////
//  Quick Z/48 probe from known A(16) examples.
//
//  A rational point of order 48 would in particular force 3 | #J(F_p)
//  at every good prime p != 3.  This script tests the known A(8)->A(16)
//  hits and the split sanity checks against that necessary condition.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
P<x> := PolynomialRing(Q);

function IntegralSquareModel(f)
    L := 1;
    for c in Coefficients(f) do
        L := LCM(L, Denominator(c));
    end for;
    return P!(L^2*f);
end function;

function A8f(rv, pv, tv)
    e := tv^2 - 2*pv*tv/rv;
    d := e + 2*pv - rv^2;
    lambda := rv/tv;
    a := rv^2 - lambda;
    b := 2*rv*pv - 2*lambda*(pv + rv*tv) + 2*rv*lambda;
    c := pv^2 + 2*pv*rv^2 - rv^4 - rv^3*tv - rv*pv^2/tv
         - lambda*(rv^2 + e) + 2*lambda*(rv*pv + rv^2*tv - 3*pv*tv + rv*tv^2);
    Qpoly := x^2 + d;
    q := a*x^2 + b*x + c;
    return q*(Qpoly^2 + q);
end function;

function IsGoodAt(f, p)
    F := GF(p);
    PF<xp> := PolynomialRing(F);
    fp := PF![F!Coefficient(f, i) : i in [0..Degree(f)]];
    if Degree(fp) lt 5 then
        return false, fp;
    end if;
    if Discriminant(fp) eq 0 then
        return false, fp;
    end if;
    return true, fp;
end function;

function FirstPrimeKilling3(f)
    for p in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97] do
        ok, fp := IsGoodAt(f, p);
        if not ok then
            continue;
        end if;
        C := HyperellipticCurve(fp);
        n := Integers()!#Jacobian(C);
        if n mod 3 ne 0 then
            return true, p, n;
        end if;
    end for;
    return false, 0, 0;
end function;

examples := [
    <"simple_A16_rt_3_1over3_p2", IntegralSquareModel(A8f(Q!3, Q!2, Q!1/3))>,
    <"simple_A16_rt_3_1over3_p34over9", IntegralSquareModel(A8f(Q!3, Q!34/9, Q!1/3))>,
    <"split_A8_plain_minus8_target_416", IntegralSquareModel(A8f(Q!-8, Q!-8, Q!-8))>,
    <"split_A8_2rank2_2216_even",
      P!(184320000000000000*x^6 - 369930240000000000*x^4
         + 183582028800000000*x^2 - 3945417984000000)>
];

print "Z48_KNOWN_A16_3FILTER";
for rec in examples do
    label := rec[1];
    f := rec[2];
    print "EXAMPLE", label;
    print "  degree", Degree(f), "disc_zero", Discriminant(f) eq 0;
    tors := [];
    try
        tors := Invariants(TorsionSubgroup(Jacobian(HyperellipticCurve(f))));
    catch e
        tors := ["torsion_failed"];
    end try;
    print "  torsion", tors;
    killed, p, n := FirstPrimeKilling3(f);
    if killed then
        print "  Z48_3PART_KILLED", "p", p, "#J", n, "mod3", n mod 3;
    else
        print "  Z48_3PART_SURVIVES_TEST_PRIMES";
    end if;
end for;

print "DONE";

