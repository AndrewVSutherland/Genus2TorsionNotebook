
//////////////////////////////////////////////////////////////////////
//  Z/2 x Z/24 via exact 2-descent on the rational A(12) chart:
//  formulation and validation.
//
//  Chart (a12_parameterization.tex / search_A2_24_from_A12.m), params
//  (p, z, r):
//    s = (z^2-4p^2+1)/(2z),  t = (z^2+4p^2-1)^2/(8p^2 z),
//    mu = ((s^2-1)(2pr+1) - p^2(2st-4))/(4p^3),
//    lambda = (4-mu^2)p^2/(s^2-1),
//    T1 = p*x + r,  R = (T1^2+x-1)/lambda,  ell = s*x+t,
//    Qp = 2*T1 + mu*R,  F = R*x^2 + 4*(R+x-1)*(R-1) = Qp^2 + R*ell^2,
//    curve  C: Y^2 = f := R*F   (degree 6),
//    P4 = [Qp^, (R*ell) mod Qp^],  P6 = [(R+x-1)^, (x*R) mod .],
//    P12 = P4 + P6 of order 12.
//
//  Target [2,24]: some 2-torsion translate D = P12 + T divisible by 2,
//  plus 2-rank >= 2.
//
//  Descent (x - T map, even degree 6):
//    delta([u,v]) = u_monic(T) in A*/(A*^2 Q*),  A = Q[T]/f_monic.
//    D in 2J(Q) ==> delta(D) trivial ==> NECESSARY squareclass
//    conditions; per component (K_R quadratic, K_F quartic):
//      N1 := Norm_{K_R}(delta(D)) = square,
//      N2 := Norm_{K_F}(delta(D)) = square,
//    computable as resultants -- symbolic, factorable rational functions
//    of (p,z,r).  Translates multiply delta by delta(T_i) (computable),
//    so "some translate halvable" = coset condition, no exact solves.
//
//  This script:
//   1. builds the chart symbolically, verifies F = Qp^2 + R*ell^2;
//   2. computes N1, N2 for delta(P12) = u_{P4}(T)*u_{P6}(T) and factors
//      them over Q(p,z,r);
//   3. VALIDATES the necessary conditions numerically: for samples,
//      N1,N2 square-tests vs Magma's IsDivisibleBy(P12 + T, 2) over the
//      translate subgroup (criterion must never reject a halvable
//      sample -- the conditions are necessary for the trivial-translate
//      component; translate handling checked separately);
//   4. reports everything needed to build the sieve.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

// ---- symbolic chart over Q(p,z,r) ----
K3<p,z,r> := RationalFunctionField(Q, 3);
PX<x> := PolynomialRing(K3);

s := (z^2 - 4*p^2 + 1)/(2*z);
t := (z^2 + 4*p^2 - 1)^2/(8*p^2*z);
mu := ((s^2 - 1)*(2*p*r + 1) - p^2*(2*s*t - 4))/(4*p^3);
lambda := (4 - mu^2)*p^2/(s^2 - 1);
T1 := p*x + r;
Rpol := (T1^2 + x - 1)/lambda;
ell := s*x + t;
Qpol := 2*T1 + mu*Rpol;
F := Rpol*x^2 + 4*(Rpol + x - 1)*(Rpol - 1);
assert F eq Qpol^2 + Rpol*ell^2;
print "chart identity F = Q^2 + R*ell^2 verified over Q(p,z,r)";
f := Rpol*F;
printf "deg f = %o, lc(f) = (rational function)\n", Degree(f);

// Mumford u-polys of P4 and P6 (monic)
u4 := Qpol/LeadingCoefficient(Qpol);
u6 := (Rpol + x - 1)/LeadingCoefficient(Rpol + x - 1);

// components of A: R-part (quadratic) and F-part (quartic), both monic
Rm := Rpol/LeadingCoefficient(Rpol);
Fm := F/LeadingCoefficient(F);

// delta(P12) = u4(T)*u6(T).  Component norms via resultants:
//   Norm_{K_g}(u(theta)) = Res_x(g_monic, u) for monic g.
// (For u monic of deg 2 and g monic: Res(g,u) = prod u(theta_i).)
N1_4 := Resultant(Rm, u4);
N1_6 := Resultant(Rm, u6);
N2_4 := Resultant(Fm, u4);
N2_6 := Resultant(Fm, u6);
N1 := N1_4*N1_6;   // Norm_{K_R}(delta(P12))
N2 := N2_4*N2_6;   // Norm_{K_F}(delta(P12))

procedure ShowFactored(name, val)
    num := Numerator(val); den := Denominator(val);
    printf "%o:\n", name;
    printf "  num: %o\n", Factorization(num);
    printf "  den: %o\n", Factorization(den);
end procedure;

print "\n== norm conditions for delta(P12) (necessary for P12 in 2J) ==";
ShowFactored("N1 = Norm_{K_R} delta(P12)", N1);
ShowFactored("N2 = Norm_{K_F} delta(P12)", N2);
// global product consistency: N1*N2 = Norm_A(delta) should be related to
// a square times powers of lc(f) etc.
ShowFactored("N1*N2", N1*N2);

// ---- numeric validation ----
P2<xq> := PolynomialRing(Q);

function IsSquareQ(qv)
    qv := Q!qv;
    if qv le 0 then return false; end if;
    return IsSquare(Numerator(qv)) and IsSquare(Denominator(qv));
end function;

function ChartData(pv, zv, rv)
    if pv eq 0 or zv eq 0 then return false, _, _, _, _, _; end if;
    sv := (zv^2 - 4*pv^2 + 1)/(2*zv);
    if sv^2 eq 1 then return false, _, _, _, _, _; end if;
    tv := (zv^2 + 4*pv^2 - 1)^2/(8*pv^2*zv);
    muv := ((sv^2 - 1)*(2*pv*rv + 1) - pv^2*(2*sv*tv - 4))/(4*pv^3);
    lav := (4 - muv^2)*pv^2/(sv^2 - 1);
    if lav eq 0 then return false, _, _, _, _, _; end if;
    T1v := pv*xq + rv;
    Rv := (T1v^2 + xq - 1)/lav;
    ellv := sv*xq + tv;
    Qv := 2*T1v + muv*Rv;
    Fv := Rv*xq^2 + 4*(Rv + xq - 1)*(Rv - 1);
    fv := Rv*Fv;
    if Degree(fv) ne 6 or Discriminant(fv) eq 0 then return false, _, _, _, _, _; end if;
    return true, fv, Rv, Qv, ellv, Fv;
end function;

// numeric norms via the same resultants
function NormsAt(pv, zv, rv)
    ok, fv, Rv, Qv, ellv, Fv := ChartData(pv, zv, rv);
    if not ok then return false, 0, 0; end if;
    Rmv := Rv/LeadingCoefficient(Rv);
    Fmv := Fv/LeadingCoefficient(Fv);
    u4v := Qv/LeadingCoefficient(Qv);
    u6v := (Rv + xq - 1)/LeadingCoefficient(Rv + xq - 1);
    n1 := Resultant(Rmv, u4v)*Resultant(Rmv, u6v);
    n2 := Resultant(Fmv, u4v)*Resultant(Fmv, u6v);
    return true, n1, n2;
end function;

print "\n== validation vs exact Jacobian divisibility ==";
// samples: check that whenever P12 itself is divisible by 2 (no translate),
// N1 and N2 are squares; and record the joint distribution.
nsamp := 0; nhalvable := 0; nviol := 0;
agreeN := 0;
for pv0 in [1,2,-1,3] do for zv0 in [1,2,3,-2] do for rv0 in [0,1,-1,2,-3] do
    pv := Q!pv0; zv := Q!zv0; rv := Q!rv0;
    ok, fv, Rv, Qv, ellv, Fv := ChartData(pv, zv, rv);
    if not ok then continue; end if;
    okn, n1, n2 := NormsAt(pv, zv, rv);
    if not okn or n1 eq 0 or n2 eq 0 then continue; end if;
    // exact side
    L := 1;
    for i in [0..Degree(fv)] do L := LCM(L, Denominator(Coefficient(fv, i))); end for;
    fI := P2!(L^2*fv);
    C := HyperellipticCurve(fI);
    J := Jacobian(C);
    O := J!0;
    u4v := Qv/LeadingCoefficient(Qv);
    u6v := (Rv + xq - 1)/LeadingCoefficient(Rv + xq - 1);
    // Mumford v-polys scale by L
    try
        P4 := J![u4v, (L*Rv*ellv) mod u4v];
        P6 := J![u6v, (L*xq*Rv) mod u6v];
    catch e continue; end try;
    P12 := P4 + P6;
    if not (12*P12 eq O and &and[n*P12 ne O : n in [1..11]]) then continue; end if;
    nsamp +:= 1;
    div2, _ := IsDivisibleBy(P12, 2);
    if div2 then
        nhalvable +:= 1;
        if not (IsSquareQ(n1) and IsSquareQ(n2)) then
            nviol +:= 1;
            printf "VIOLATION p=%o z=%o r=%o: halvable but N1sq=%o N2sq=%o\n",
                pv, zv, rv, IsSquareQ(n1), IsSquareQ(n2);
        end if;
    end if;
    // track how discriminating the necessary conditions are
    if (IsSquareQ(n1) and IsSquareQ(n2)) eq div2 then agreeN +:= 1; end if;
end for; end for; end for;
printf "samples=%o exactly-halvable=%o necessary-condition violations=%o\n",
    nsamp, nhalvable, nviol;
printf "(exact <=> N1,N2 both square) agreement: %o/%o  (not expected to be 100%%:\n", agreeN, nsamp;
print "  N-conditions are necessary, not sufficient; second stage narrows the rest)";

// finite-field positive controls: over F_q, count points where P12 is
// divisible by 2 and confirm N1,N2 are squares there
for q0 in [11, 13] do
    Fq := GF(q0);
    PF<xf> := PolynomialRing(Fq);
    tested := 0; halv := 0; viol := 0;
    for pv in Fq do for zv in Fq do for rv in Fq do
        if pv eq 0 or zv eq 0 then continue; end if;
        sv := (zv^2 - 4*pv^2 + 1)/(2*zv);
        if sv^2 eq 1 then continue; end if;
        tv := (zv^2 + 4*pv^2 - 1)^2/(8*pv^2*zv);
        muv := ((sv^2 - 1)*(2*pv*rv + 1) - pv^2*(2*sv*tv - 4))/(4*pv^3);
        den := sv^2 - 1;
        lav := (4 - muv^2)*pv^2/den;
        if lav eq 0 then continue; end if;
        T1v := pv*xf + rv;
        Rv := (T1v^2 + xf - 1)/lav;
        ellv := sv*xf + tv;
        Qv := 2*T1v + muv*Rv;
        Fv := Rv*xf^2 + 4*(Rv + xf - 1)*(Rv - 1);
        fv := Rv*Fv;
        if Degree(fv) ne 6 or not IsSquarefree(fv) then continue; end if;
        u4v := Qv/LeadingCoefficient(Qv);
        u6v := (Rv + xf - 1)/LeadingCoefficient(Rv + xf - 1);
        Rmv := Rv/LeadingCoefficient(Rv);
        Fmv := Fv/LeadingCoefficient(Fv);
        if Resultant(Rmv, u4v) eq 0 or Resultant(Rmv, u6v) eq 0
           or Resultant(Fmv, u4v) eq 0 or Resultant(Fmv, u6v) eq 0 then continue; end if;
        try
            J := Jacobian(HyperellipticCurve(fv));
            P4 := J![u4v, (Rv*ellv) mod u4v];
            P6 := J![u6v, (xf*Rv) mod u6v];
        catch e continue; end try;
        P12 := P4 + P6;
        O := J!0;
        if not (12*P12 eq O and &and[n*P12 ne O : n in [1..11]]) then continue; end if;
        tested +:= 1;
        // divisibility on finite Jacobian
        G, phi := AbelianGroup(J);
        aG := P12 @@ phi;
        coords := Eltseq(aG); invs := Invariants(G);
        div2 := true;
        for i in [1..#coords] do
            g := GCD(2, invs[i]);
            if coords[i] mod g ne 0 then div2 := false; break; end if;
        end for;
        if not div2 then continue; end if;
        halv +:= 1;
        n1 := Resultant(Rmv, u4v)*Resultant(Rmv, u6v);
        n2 := Resultant(Fmv, u4v)*Resultant(Fmv, u6v);
        if not (IsSquare(n1) and IsSquare(n2)) then
            viol +:= 1;
            if viol le 3 then
                printf "F_%o VIOLATION p=%o z=%o r=%o n1sq=%o n2sq=%o\n",
                    q0, pv, zv, rv, IsSquare(n1), IsSquare(n2);
            end if;
        end if;
    end for; end for; end for;
    printf "F_%o: order-12 samples=%o halvable=%o necessary-condition violations=%o\n",
        q0, tested, halv, viol;
end for;
print "DONE";
quit;
