//////////////////////////////////////////////////////////////////////
//  Symbolic HPL check for simultaneous Fi/Fj genus-one fibers.
//
//  On the D-square + F0-square surface, write sigma = q^2*rho/tau.
//  For fixed (rho,tau), each target T in {F1,F2,F4} has the form
//
//      y^2 = A_T*(1 + B_T*q^2).
//
//  Given a point (q0,y0) on the first conic, a slope m gives
//
//      q(m) = (q0*m^2 - 2*y0*m + q0*K)/(m^2 - K),
//      K = A_1*B_1.
//
//  Imposing a second conic gives the reciprocal quartic
//
//      W^2 = A_2*((m^2-K)^2
//                 + B_2*(q0*m^2 - 2*y0*m + q0*K)^2).
//
//  Equivalently, with z = m + K/m and W = m*Z, it is the fiber
//  product of two conics:
//
//      Z^2 = A_2*((z^2 - 4*K) + B_2*(q0*z - 2*y0)^2),
//      u^2 = z^2 - 4*K.
//
//  This script verifies this structure at the HPL point and optionally
//  constructs elliptic models for the three pair fibers.
//////////////////////////////////////////////////////////////////////

Q := Rationals();
R<m> := PolynomialRing(Q);

if not assigned do_elliptic then
    do_elliptic := true;
elif Type(do_elliptic) eq MonStgElt then
    do_elliptic := do_elliptic eq "true" or do_elliptic eq "1";
end if;

if not assigned do_rank then
    do_rank := false;
elif Type(do_rank) eq MonStgElt then
    do_rank := do_rank eq "true" or do_rank eq "1";
end if;

function IsNonzeroSquareQ(q)
    if q eq 0 then
        return false, Q!0;
    end if;
    num := Integers()!Numerator(q);
    den := Integers()!Denominator(q);
    if num lt 0 then
        return false, Q!0;
    end if;
    okN, rtN := IsSquare(num);
    if not okN then return false, Q!0; end if;
    okD, rtD := IsSquare(den);
    if not okD then return false, Q!0; end if;
    return true, Q!rtN/Q!rtD;
end function;

function TargetConicAB(rho, tau, tgt)
    if tgt eq "F1" then
        return true, (1 + rho)*(1 + tau), rho/tau;
    elif tgt eq "F2" then
        return true, (1 + rho)*(rho + tau), Q!1/tau;
    elif tgt eq "F4" then
        return true, (1 + tau)*(rho + tau), rho/(tau^2);
    end if;
    return false, Q!0, Q!0;
end function;

function PairQuartic(A1, B1, A2, B2, q0, y0)
    K := A1*B1;
    den := m^2 - K;
    num := q0*m^2 - 2*y0*m + q0*K;
    return A2*(den^2 + B2*num^2);
end function;

function SquarefreePart(f)
    if f eq 0 then
        return R!0;
    end if;
    fac := Factorization(f);
    sf := R!1;
    for ff in fac do
        if ff[2] mod 2 eq 1 then
            sf *:= ff[1];
        end if;
    end for;
    return sf;
end function;

function GenusFromSquarefreeDegree(deg)
    if deg le 0 then
        return -1;
    end if;
    return Floor((deg - 1)/2);
end function;

function CoeffOrZero(f, i)
    if i gt Degree(f) then
        return Q!0;
    end if;
    return Coefficient(f, i);
end function;

function CoverValues(rho, sigma, tau)
    F0 := rho*sigma*tau;
    F1 := (1 + rho)*(1 + sigma)*(1 + tau);
    F2 := rho*(1 + rho)*(rho + sigma)*(rho + tau);
    F3 := sigma*(1 + sigma)*(rho + sigma)*(sigma + tau);
    F4 := tau*(1 + tau)*(rho + tau)*(sigma + tau);
    return F0, F1, F2, F3, F4;
end function;

pairs := [
    <"F1","F2">,
    <"F1","F4">,
    <"F2","F4">
];

rhoH := Q!58466134224 / Q!53109477625;
sigmaH := Q!719363573659505664 / Q!749082246897952705;
tauH := Q!307598400 / Q!352612321;
qH := tauH;

print "M2248 HPL simultaneous pair-fiber symbolic check";
print "do_elliptic", do_elliptic, "do_rank", do_rank;
print "rho", rhoH;
print "sigma", sigmaH;
print "tau", tauH;
print "q", qH;

F0H, F1H, F2H, F3H, F4H := CoverValues(rhoH, sigmaH, tauH);
assert IsNonzeroSquareQ(F0H);
assert IsNonzeroSquareQ(F1H);
assert IsNonzeroSquareQ(F2H);
assert IsNonzeroSquareQ(F3H);
assert IsNonzeroSquareQ(F4H);

for pair in pairs do
    first := pair[1];
    second := pair[2];
    ok1, A1, B1 := TargetConicAB(rhoH, tauH, first);
    ok2, A2, B2 := TargetConicAB(rhoH, tauH, second);
    assert ok1 and ok2;

    okY1, y1 := IsNonzeroSquareQ(A1*(1 + B1*qH^2));
    okY2, y2 := IsNonzeroSquareQ(A2*(1 + B2*qH^2));
    assert okY1 and okY2;

    K := A1*B1;
    f := PairQuartic(A1, B1, A2, B2, qH, y1);
    sf := SquarefreePart(f);
    deg := Degree(sf);
    genus := GenusFromSquarefreeDegree(deg);

    a4 := CoeffOrZero(f, 4);
    a3 := CoeffOrZero(f, 3);
    a2 := CoeffOrZero(f, 2);
    a1 := CoeffOrZero(f, 1);
    a0 := CoeffOrZero(f, 0);

    print "PAIR", first, second;
    print "A1", A1;
    print "B1", B1;
    print "A2", A2;
    print "B2", B2;
    print "K", K;
    print "quartic_degree", Degree(f), "squarefree_degree", deg, "genus", genus;
    print "reciprocal_check_a0_eq_a4K2", a0 eq a4*K^2;
    print "reciprocal_check_a1_eq_a3K", a1 eq a3*K;
    print "z_conic_coefficients",
          "z2", A2*(1 + B2*qH^2),
          "z1", -4*A2*B2*qH*y1,
          "z0", A2*(4*B2*y1^2 - 4*K);

    if do_elliptic then
        tangent_m := K*qH/y1;
        den := tangent_m^2 - K;
        W := y2*den;
        C := HyperellipticCurve(f);
        P0 := C![tangent_m, W, 1];
        assert P0 in C;
        E, phi := EllipticCurve(C, P0);
        print "base_slope", tangent_m;
        print "elliptic_model", E;
        T, tors_map := TorsionSubgroup(E);
        print "torsion_invariants", Invariants(T);
        if do_rank then
            lo, hi := RankBounds(E);
            print "rank_bounds", lo, hi;
        end if;
    end if;
end for;

print "DONE";
