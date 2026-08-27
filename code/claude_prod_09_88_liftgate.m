// (8,8) lane, steps 1.1(base)+1.2: per-base symbolic lift layer AND genus gate.
// At a fixed stage-1 base (m,n) (so s,t rational, v free), over K = Q(v):
//   1. solve the halving identity for the order-4 sections over T_j = [l_j]:
//      monic(lj)*(mu*x+nu)^2 - gamma_j*(other two l's) = kappa*u^2
//      (two closure equations in (mu,nu); rational solutions <=> sections),
//      u_D = x^2 + p(v) x + q(v);
//   2. VALIDATE the sections against exact J1 torsion at a rational v anchor;
//   3. lambda layer: val_j = u_D(theta_j) = Qt*theta' + Pt (theta'^2 = e_j),
//      norm nt_j = Pt^2 - e_j*Qt^2; if nt_j = w^2 in K the two lambda classes
//      are lambda_j^{+-} = 2(Pt +- w); D in 2*J1(Q) at v0 requires a COMMON
//      rational square class among {lambda_j^{+-}}_{j=1,2,3};
//   4. GENUS GATE: for each section and each sign combo, the lift locus is
//      the set of rational points on the curve
//        y12^2 = lambda_1 * lambda_2 ,   y13^2 = lambda_1 * lambda_3
//      over the v-line; compute its exact genus over Q(v) via function
//      fields.  Genus 0 => parametrize; genus 1 => rank; genus >= 2 => wall.
//
// Run: magma -b Mnum:=2 Mden:=1 Nnum:=1 Nden:=1 Lj:=3 Vanchor:=1 code/claude_prod_09_88_liftgate.m

SetColumns(0);
SetSeed(1);

if not assigned Mnum then Mnum := 2; elif Type(Mnum) eq MonStgElt then Mnum := StringToInteger(Mnum); end if;
if not assigned Mden then Mden := 1; elif Type(Mden) eq MonStgElt then Mden := StringToInteger(Mden); end if;
if not assigned Nnum then Nnum := 1; elif Type(Nnum) eq MonStgElt then Nnum := StringToInteger(Nnum); end if;
if not assigned Nden then Nden := 1; elif Type(Nden) eq MonStgElt then Nden := StringToInteger(Nden); end if;
if not assigned Lj then Lj := 3; elif Type(Lj) eq MonStgElt then Lj := StringToInteger(Lj); end if;
if not assigned Vanchor then Vanchor := 1; elif Type(Vanchor) eq MonStgElt then Vanchor := StringToInteger(Vanchor); end if;
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
// Scan:=1 -> skip anchor validation and genus; enumerate rational v of
// height <= Hv against the lambda-compatibility and print LIFTHIT lines.
if not assigned Scan then Scan := 0; elif Type(Scan) eq MonStgElt then Scan := StringToInteger(Scan); end if;
if not assigned Hv then Hv := 60; elif Type(Hv) eq MonStgElt then Hv := StringToInteger(Hv); end if;
SetMemoryLimit(MemGB*10^9);

// loaded at top level (load is a file-scope directive in Magma); provides
// SqfPart, IntSextic, Lambda334, StageOneST, ExactTorsion, P88
load "code/claude_prod_09_88_defs.m";

QQ := Rationals();
m0 := QQ!Mnum/Mden; n0 := QQ!Nnum/Nden;
t0 := (m0^2+1)/(2*m0); al0 := (m0^2-1)/(2*m0); s0 := (n0^2+al0^4)/(2*n0);
printf "BASE m=%o n=%o s=%o t=%o TARGETCLASS l%o\n", m0, n0, s0, t0, Lj;

K<vv> := RationalFunctionField(QQ);

// Lambda_334 data over K
A := s0^2 - t0^4 + t0^2;
den_u := -s0^2*t0*A*vv^2 + t0;
u334 := (-s0^2*A*vv^2 - 2*A*vv - 1) / den_u;
a := A/(1 - t0^2);
b := A/(u334^2*s0^2 + 1 - t0^2);
c := t0^2;
d2 := A * (s0^2*u334^2 + t0^4 - 2*t0^2 + 1)
        * (s0^4*u334^2 - s0^2*t0^2*u334^2 + s0^2*u334^2 - t0^6 + 3*t0^4 - 3*t0^2 + 1);

PR<mu,nu> := PolynomialRing(K, 2);
Px<x> := PolynomialRing(PR);
l1 := (-a+b+c-1)*x^2 + (2*a - 2*b*c)*x + (a*b*c - a*b - a*c + b*c);
l2 := -x^2 + b*c;
l3 := x^2 - a;
ljs := [l1, l2, l3];
lead123 := K!(LeadingCoefficient(l1)*LeadingCoefficient(l2)*LeadingCoefficient(l3));
// g1 = d2*f1 = (d2/lead123) * l1*l2*l3
gam := d2/lead123;

// halving of T_{Lj}: mlj*(mu*x+nu)^2 - gam*lc(lj)*(other two) = kappa*u^2
lj := ljs[Lj];
mlj := lj/LeadingCoefficient(lj);
others := &*[Px| ljs[k] : k in [1..3] | k ne Lj];
Q4 := mlj*(mu*x + nu)^2 - (Px!(gam*K!LeadingCoefficient(lj)))*others;
assert Degree(Q4) le 4;
e4 := Coefficient(Q4,4); e3 := Coefficient(Q4,3);
e2 := Coefficient(Q4,2); e1 := Coefficient(Q4,1); e0 := Coefficient(Q4,0);
F1 := 8*e4^2*e1 - e3*(4*e4*e2 - e3^2);
F2 := 64*e4^3*e0 - (4*e4*e2 - e3^2)^2;
R := Resultant(F1, F2, nu);
printf "RESULTANT degree_in_mu=%o\n", Degree(R, mu);
Pmu := PolynomialRing(K);
Rmu := &+[Pmu| K!Coefficient(R, mu, i)*Pmu.1^i : i in [0..Degree(R, mu)]];
Rmu := Rmu/LeadingCoefficient(Rmu);
fac := Factorization(Rmu);
printf "RES_FACTOR_DEGREES %o\n", [<Degree(fp[1]), fp[2]> : fp in fac];
mucands := [ -Coefficient(fp[1],0) : fp in fac | Degree(fp[1]) eq 1 ];
mucands := [ mv : mv in mucands | mv ne 0 ];

sections := [];
for muv in mucands do
    F1n := UnivariatePolynomial(Evaluate(F1, mu, muv));
    F2n := UnivariatePolynomial(Evaluate(F2, mu, muv));
    if F1n eq 0 then continue; end if;
    for r in Roots(F1n) do
        nuv := r[1];
        if Evaluate(F2n, nuv) ne 0 then continue; end if;
        ev := func< e | K!Evaluate(Evaluate(e, mu, muv), nu, nuv) >;
        k4 := ev(e4); k3 := ev(e3); k2 := ev(e2);
        if k4 eq 0 then continue; end if;
        pp := k3/(2*k4);
        qq := (4*k4*k2 - k3^2)/(8*k4^2);
        if forall{ sec : sec in sections | sec[1] ne pp or sec[2] ne qq } then
            Append(~sections, <pp, qq>);
        end if;
    end for;
end for;
printf "N_SECTIONS %o\n", #sections;
error if #sections eq 0, "no rational sections over this base";

// VALIDATION at v = Vanchor (skipped in scan mode)
if Scan eq 0 then
sv, tv := StageOneST(m0, n0);
h, g1m, am, bm, cm := Lambda334(sv, tv, QQ!Vanchor);
ljm := [P88| (-am+bm+cm-1)*P88.1^2 + (2*am-2*bm*cm)*P88.1 + (am*bm*cm-am*bm-am*cm+bm*cm),
             -P88.1^2 + bm*cm, P88.1^2 - am ];
mljm := ljm[Lj]/LeadingCoefficient(ljm[Lj]);
J1 := Jacobian(HyperellipticCurve(IntSextic(g1m)));
Tg, mp := TorsionSubgroup(J1);
assert Invariants(Tg) eq [4,4];
Tpt := J1![mljm, P88!0];
utrue := { mp(g)[1] : g in Tg | Order(g) eq 4 and 2*mp(g) eq Tpt };
uspec := { P88| P88.1^2 + Evaluate(sec[1], QQ!Vanchor)*P88.1 + Evaluate(sec[2], QQ!Vanchor)
              : sec in sections };
utrue2 := { P88![Coefficient(pol,i) : i in [0..2]] : pol in utrue };
printf "VALIDATION_TRUE_U %o\nVALIDATION_SPEC_U %o\n", utrue2, uspec;
assert uspec eq utrue2;
printf "SECTIONS_VALIDATED_AT v=%o\n", Vanchor;
end if;

// lambda layer + genus gate per section
for si in [1..#sections] do
    pp := sections[si][1]; qq := sections[si][2];
    lams := [* *];   // per j: list of lambda functions (1 or 2), in K or ext
    allrat := true;
    lamfns := [];    // when all nt_j are squares in K: sequence of [lam+, lam-]
    for j in [1..3] do
        ljK := ljs[j];  // in Px over PR; coefficients constant in mu,nu
        bet := K!(Coefficient(ljK,1)/Coefficient(ljK,2));
        gamj := K!(Coefficient(ljK,0)/Coefficient(ljK,2));
        ej := bet^2/4 - gamj;
        Qt := pp - bet;
        Pt := qq - gamj - bet*(pp - bet)/2;
        nt := Pt^2 - ej*Qt^2;
        issq, w := IsSquare(nt);
        printf "SECTION %o j=%o NORM_IDENTICALLY_SQUARE %o\n", si, j, issq;
        if not issq then allrat := false; break; end if;
        Append(~lamfns, [2*(Pt + w), 2*(Pt - w)]);
    end for;
    if not allrat then
        printf "SECTION %o: norm not identically square -- cover has an extra layer; SKIP (handle in follow-up)\n", si;
        continue;
    end if;
    if Scan eq 1 then
        // enumerate rational v of height <= Hv against lambda-compatibility
        nhits := 0;
        for hn in [1..Hv] do for hd in [1..Hv] do
            if GCD(hn, hd) ne 1 then continue; end if;
            for sg in [1, -1] do
                v0 := QQ!(sg*hn)/hd;
                lv := [];
                bad := false;
                for j in [1..3] do
                    lp := 0; lm := 0;
                    try
                        lp := Evaluate(lamfns[j][1], v0);
                        lm := Evaluate(lamfns[j][2], v0);
                    catch e bad := true; end try;
                    if bad or (lp eq 0 and lm eq 0) then bad := true; break; end if;
                    Append(~lv, [lp, lm]);
                end for;
                if bad then continue; end if;
                for s1 in [1,2] do for s2 in [1,2] do for s3 in [1,2] do
                    a1 := lv[1][s1]; a2 := lv[2][s2]; a3 := lv[3][s3];
                    if a1 eq 0 or a2 eq 0 or a3 eq 0 then continue; end if;
                    if IsSquare(a1*a2) and IsSquare(a1*a3) then
                        nhits +:= 1;
                        printf "LIFTHIT m=%o n=%o v=%o sec=%o signs=%o%o%o\n",
                               m0, n0, v0, si, s1, s2, s3;
                    end if;
                end for; end for; end for;
            end for;
        end for; end for;
        printf "SCAN_SECTION %o hits=%o Hv=%o\n", si, nhits, Hv;
        continue;
    end if;
    // genus of y12^2 = lam1*lam2, y13^2 = lam1*lam3 for all sign combos
    for s1 in [1,2] do for s2 in [1,2] do for s3 in [1,2] do
        A12 := lamfns[1][s1]*lamfns[2][s2];
        A13 := lamfns[1][s1]*lamfns[3][s3];
        if A12 eq 0 or A13 eq 0 then continue; end if;
        // squarefree polynomial representative of the square class of r in
        // K = Q(v).  NB Factorization drops the unit (lab convention #1) --
        // reconstruct it: it carries the square class of the constant part.
        sqA := function(r)
            nm := Numerator(r)*Denominator(r);   // in the poly ring under K
            fac := Factorization(nm);
            unit := nm div &*[Parent(nm)| fp[1]^fp[2] : fp in fac];
            assert Degree(unit) eq 0;
            uc := SqfPart(Rationals()!Coefficient(unit, 0));
            return uc * &*[Parent(nm)| fp[1]^(fp[2] mod 2) : fp in fac];
        end function;
        B12 := sqA(A12); B13 := sqA(A13);
        // a degree-0 condition is a CONSTANT square class: identically
        // satisfied iff that constant is a square (then no cover layer),
        // identically obstructed otherwise.
        con12 := Degree(B12) eq 0; con13 := Degree(B13) eq 0;
        if con12 and not IsSquare(QQ!Coefficient(B12,0)) then
            printf "GATE sec=%o signs=%o%o%o CONSTANT_OBSTRUCTION lam12=%o\n",
                   si, s1, s2, s3, B12; continue;
        end if;
        if con13 and not IsSquare(QQ!Coefficient(B13,0)) then
            printf "GATE sec=%o signs=%o%o%o CONSTANT_OBSTRUCTION lam13=%o\n",
                   si, s1, s2, s3, B13; continue;
        end if;
        if con12 and con13 then
            printf "GATE sec=%o signs=%o%o%o IDENTICALLY_SATISFIED (lift condition holds on the whole v-line!)\n",
                   si, s1, s2, s3;
            continue;
        end if;
        gtot := -1;
        try
            FKv<xv> := FunctionField(QQ);
            c12 := &+[FKv| (QQ!Coefficient(B12,i))*xv^i : i in [0..Degree(B12)]];
            c13 := &+[FKv| (QQ!Coefficient(B13,i))*xv^i : i in [0..Degree(B13)]];
            if Degree(B12) gt 0 then
                PF := PolynomialRing(FKv);
                L1 := FunctionField(PF.1^2 - c12);
            else
                L1 := FKv;
            end if;
            PL1 := PolynomialRing(L1);
            if Degree(B13) gt 0 then
                L2 := FunctionField(PL1.1^2 - L1!c13);
            else
                L2 := L1;
            end if;
            gtot := Genus(L2);
        catch e
            printf "GATE sec=%o signs=%o%o%o GENUS_ERROR %o\n", si, s1, s2, s3, e`Object;
            continue;
        end try;
        printf "GATE sec=%o signs=%o%o%o deg12=%o deg13=%o GENUS %o\n",
               si, s1, s2, s3, Degree(B12), Degree(B13), gtot;
    end for; end for; end for;
end for;
printf "LIFTGATE_DONE base m=%o n=%o class=l%o\n", m0, n0, Lj;
quit;
