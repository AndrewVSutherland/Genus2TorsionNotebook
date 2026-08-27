// BLP [2,22] engine: rational points on the RESIDUAL locus at fixed (r,s).
//
// Background (results/blp22_locus_rm_test2.log + rmdebug2, 2026-08-01):
// the closure numerators N0, N1 in (q,d) share a d-free quartic gcd G whose
// points are DEGENERATE (disc = 0) -- a spurious component.  The residual
// system N0/G = N1/G = 0 is the true [1,1,2,2]+order-11 locus, and mod-p
// sampling shows its generic point is NOT RM-sqrt5 (34/36 incompatible).
// So rational residual points generically give NON-RM [2,22] curves.
//
// This script: for fixed rational (r,s), compute Res_q(N0/G, N1/G) as a
// polynomial in d, find its RATIONAL roots, then common q-roots, and for
// each resulting curve: disc/genus checks, CFOrder (independent order-11
// certificate), exact TorsionSubgroup, and print HIT22 lines.
//
// Run: magma -b Rnum:=1 Rden:=3 Snum:=2 Sden:=3 code/blp22_residual_solve.m

SetColumns(0);
SetSeed(1);
if not assigned Rnum then Rnum := 1; elif Type(Rnum) eq MonStgElt then Rnum := StringToInteger(Rnum); end if;
if not assigned Rden then Rden := 3; elif Type(Rden) eq MonStgElt then Rden := StringToInteger(Rden); end if;
if not assigned Snum then Snum := 2; elif Type(Snum) eq MonStgElt then Snum := StringToInteger(Snum); end if;
if not assigned Sden then Sden := 3; elif Type(Sden) eq MonStgElt then Sden := StringToInteger(Sden); end if;
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
r := Q!Rnum/Rden; sv := Q!Snum/Sden;
error if r eq sv, "need r != s";
printf "RS r=%o s=%o\n", r, sv;

load "code/blp22_locus_rm_test2.m0";   // DeriveT

K2<qv, dv> := RationalFunctionField(Q, 2);
P2 := PolynomialRing(Q, 2);
PQ<X> := PolynomialRing(Q);

rK := K2!r; sK := K2!sv;
den := 2*(sK^2 + dv) - sK*(sK - rK);
error if den eq 0, "degenerate (r,s)";
pf := (2*(rK-1)*(sK^2+dv) + (sK-rK)*(sK^2+qv)) / den;
cf := (rK - 1 - pf)/2;
t0, t1, ok := DeriveT(qv - rK*pf, -rK*qv + 2*cf*dv, cf, K2!dv);
error if not ok or t0 eq 0 or t1 eq 0, "derivation degenerate at this (r,s)";
N0 := P2!Numerator(t0); N1 := P2!Numerator(t1);
G := GCD(N0, N1);
printf "GCD degq=%o degd=%o\n", Degree(G, 1), Degree(G, 2);
N0r := N0 div G; N1r := N1 div G;
printf "RESIDUAL degq %o/%o degd %o/%o\n",
       Degree(N0r,1), Degree(N1r,1), Degree(N0r,2), Degree(N1r,2);
Rd := Resultant(N0r, N1r, P2.1);
error if Rd eq 0, "unexpected common factor in residual system";
Rdn := UnivariatePolynomial(Rd);
printf "RES_DEG_D %o\n", Degree(Rdn);
rts := Roots(Rdn);
printf "RATIONAL_D_ROOTS %o\n", [rt[1] : rt in rts];

// CFOrder for the independent D_inf check
function SqrtPolyPart6(f)
    P := Parent(f); xx := P.1;
    sp := xx^3;
    for k in [1..3] do
        dd := f - sp^2;
        if Degree(dd) le 2 then break; end if;
        sp := sp + (Coefficient(dd, 6-k)/(2*Coefficient(sp, 3)))*xx^(3-k);
    end for;
    return sp;
end function;
function CFOrderB(f, maxsteps, maxord)
    P := Parent(f);
    sp := SqrtPolyPart6(f);
    Pi := P!0; Qi := P!1; total := 0;
    for i in [0..maxsteps] do
        if Qi eq 0 then return 0; end if;
        ai := (Pi + sp) div Qi;
        total +:= Degree(ai);
        if total gt maxord then return 0; end if;
        Pn := ai*Qi - Pi;
        if (f - Pn^2) mod Qi ne 0 then return 0; end if;
        Qn := (f - Pn^2) div Qi;
        Pi := Pn; Qi := Qn;
        if i ge 1 and Degree(Qi) le 0 and Qi ne 0 then return total; end if;
    end for;
    return 0;
end function;
function ExactTorsionInvs(f)
    C := HyperellipticCurve(f);
    try
        Cm := ReducedMinimalWeierstrassModel(C);
        C := SimplifiedModel(Cm);
    catch e ; end try;
    return Invariants(TorsionSubgroup(Jacobian(C)));
end function;

for rt in rts do
    d0 := rt[1];
    if sv^2 + d0 eq 0 then continue; end if;
    A0 := UnivariatePolynomial(Evaluate(N0r, 2, d0));
    A1 := UnivariatePolynomial(Evaluate(N1r, 2, d0));
    if A0 eq 0 or A1 eq 0 then continue; end if;
    gg := GCD(PQ!A0, PQ!A1);
    printf "D_ROOT d=%o q_gcd_deg=%o\n", d0, Degree(gg);
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
        if Degree(f) ne 6 or Discriminant(f) eq 0 then continue; end if;
        cford := CFOrderB(f, 60, 45);
        dtype := {* Degree(fp[1])^^fp[2] : fp in Factorization(f) *};
        printf "CAND d=%o q=%o cford=%o type=%o\n", d0, q0, cford, dtype;
        if cford ne 0 and cford mod 11 ne 0 then continue; end if;
        invs := ExactTorsionInvs(f);
        printf "EXACT d=%o q=%o TORSION=%o\n", d0, q0, invs;
        if #invs ge 2 and invs[#invs] mod 11 eq 0 and invs[#invs-1] mod 2 eq 0 then
            printf "HIT22 r=%o s=%o d=%o q=%o TORSION=%o f=%o\n",
                   r, sv, d0, q0, invs, f;
        end if;
    end for;
end for;
printf "RESIDUAL_SOLVE_DONE r=%o s=%o\n", r, sv;
quit;
