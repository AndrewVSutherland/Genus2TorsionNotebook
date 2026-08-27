
//////////////////////////////////////////////////////////////////////
//  Hessian / weighted blowup of the degenerate stratum <1,0> for the
//  [4,16] second-halving cover at p=7.
//
//  Level 2 (code/agent_m18_416_p7_blowup_level2.m) could not resolve
//  (R,w) == (1,0) mod 7 because the aux-Jacobian J_aux vanishes
//  IDENTICALLY mod 7 there (nullity 5), so the linear Newton step
//  carries no information.  This script resolves the stratum using the
//  exact quadratic structure:
//
//  E416 is EXACTLY QUADRATIC in aux=(a,b,c,d,e), so over Z
//      E(aux0 + 7*sigma) = E(aux0) + 7*Grad(aux0).sigma + 49*Qf(sigma)
//  exactly, with Qf the pure aux-quadratic part (independent of aux0;
//  mod 7 it depends only on (R,w) mod 7 = (1,0), hence is CONSTANT on
//  the stratum).
//
//  Since Grad(aux0) == 0 mod 7 on the stratum:
//   - mod 49 solvability depends only on aux0:  E(R,w,aux0) == 0 mod 49,
//     a condition on (R,w) mod 49;
//   - mod 343 needs sigma in F_7^5 with
//         E(R,w,aux0)/49 + G'(aux0).sigma + Qf(sigma) == 0  (mod 7),
//     where G' = Grad(aux0)/7 mod 7.  The dependence of the constant
//     term on the third digits (rho2,omega2) of (R,w) is linear:
//         c(rho2,omega2) = c0 + dR*rho2 + dw*omega2,
//     so one sigma-scan per (rho1,omega1,aux0) survivor handles all 49
//     third digits at once.
//   - SMOOTHNESS: at a solution sigma, M = Grad(aux0+7*sigma)/7 mod 7.
//     If M is nonsingular, Newton converges at every further level, so
//     that 7-adic (R,w) branch carries a genuine Q_7 point of the
//     [4,16] cover ("depth-2 Hensel-smooth").
//
//  Output: mod-49 and mod-343 survivor counts for the stratum, and the
//  count of depth-2 smooth branches.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();
p := 7;
Fp := GF(p);

// ---- build cleared integer E416 ----
Rng<R,w,a,b,c,d,e> := PolynomialRing(Q, 7, "grevlex");
KF := FieldOfFractions(Rng); PX<x> := PolynomialRing(KF);
t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
A := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
B := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);
f := x*A*B; c4 := R + 2 + 4*t;
q416 := x^2 + a*x + b; ell416 := c*x^2 + d*x + e;
F416 := f - ell416^2 - c4*(x+R)*q416^2;
E416 := [];
for i in [0..4] do
    ci := Numerator(Coefficient(F416, i));
    den := LCM([Denominator(co) : co in Coefficients(ci)]); ci := den*ci;
    g := GCD([Z!co : co in Coefficients(ci)]);
    Append(~E416, Rng!(ci/g));
end for;
DerAux := [[Derivative(E416[i], 2+j) : j in [1..5]] : i in [1..5]];  // vars a..e = 3..7
DerR := [Derivative(E416[i], 1) : i in [1..5]];
DerW := [Derivative(E416[i], 2) : i in [1..5]];

// sanity: E416 is quadratic in aux
for i in [1..5] do
    degs := [Degree(E416[i], j) : j in [3..7]];
    assert Max(degs) le 2;
end for;
print "E416 aux-degree check passed (quadratic in a..e)";

// ---- pure aux-quadratic part Qf mod 7 at (R,w)=(1,0) ----
// Qf_i(sigma) = sum of aux-degree-2 monomials of E416_i, coefficients
// evaluated at (R,w)=(1,0), reduced mod 7.
S5<s1,s2,s3,s4b,s5b> := PolynomialRing(Fp, 5);
sig := [s1,s2,s3,s4b,s5b];
Qf := [S5!0 : i in [1..5]];
for i in [1..5] do
    for mon in Terms(E416[i]) do
        expv := Exponents(mon);
        auxdeg := &+[expv[j] : j in [3..7]];
        if auxdeg ne 2 then continue; end if;
        cf := Z!(MonomialCoefficient(E416[i], Monomial(Rng, expv)));
        // evaluate the (R,w)-part at (1,0): R^e1 * w^e2 -> 1 if e2=0 else 0
        if expv[2] ne 0 then continue; end if;   // w=0 kills it
        term := (Fp!cf);
        tm := S5!term;
        for j in [3..7] do
            tm *:= sig[j-2]^expv[j];
        end for;
        Qf[i] +:= tm;
    end for;
end for;

// lookup table Qvals over all 7^5 sigma
print "building Qf lookup table (7^5 sigma) ...";
sigmaList := [];
for cnt in [0..p^5-1] do
    x0 := cnt; sv := [];
    for _ in [1..5] do Append(~sv, x0 mod p); x0 := x0 div p; end for;
    Append(~sigmaList, sv);
end for;
Qvals := [];
for sv in sigmaList do
    fv := [Fp!s : s in sv];
    Append(~Qvals, [Evaluate(Qf[i], fv) : i in [1..5]]);
end for;
print "table done";

// ---- the 49 aux solutions mod 7 at (R,w)=(1,0) ----
aux7 := [];
for aa in [0..6] do for bb in [0..6] do for cv in [0..6] do
  for dv in [0..6] do for ev in [0..6] do
    v := [1,0,aa,bb,cv,dv,ev]; ok := true;
    for i in [1..5] do
        if (Z!Evaluate(E416[i], v)) mod p ne 0 then ok := false; break; end if;
    end for;
    if ok then Append(~aux7, [aa,bb,cv,dv,ev]); end if;
end for; end for; end for; end for; end for;
printf "aux solutions mod 7 at (1,0): %o\n", #aux7;

// verify J_aux == 0 mod 7 at each
for s in aux7 do
    base := [1,0] cat s;
    for i in [1..5] do for j in [1..5] do
        assert (Z!Evaluate(DerAux[i][j], base)) mod p eq 0;
    end for; end for;
end for;
print "J_aux == 0 mod 7 verified at all aux solutions";

// ---- stage 1: mod-49 gate ----
// (R,w) = (1+7*rho1, 7*omega1); survivor iff E(R,w,aux0) == 0 mod 49
// for some aux0 (independent of the aux second digit since J_aux==0 mod 7).
surv49 := [];   // entries <rho1,omega1,aux0index>
for rho1 in [0..6] do for om1 in [0..6] do
    Rv := 1 + 7*rho1; Wv := 7*om1;
    for ai in [1..#aux7] do
        base := [Rv, Wv] cat aux7[ai];
        ok := true;
        for i in [1..5] do
            if (Z!Evaluate(E416[i], base)) mod 49 ne 0 then ok := false; break; end if;
        end for;
        if ok then Append(~surv49, <rho1, om1, ai>); end if;
    end for;
end for; end for;
rw49 := Seqset([<s[1],s[2]> : s in surv49]);
printf "STAGE1 mod-49: survivor (rho1,omega1) pairs = %o / 49, (pair,aux0) branches = %o\n",
    #rw49, #surv49;

// ---- stage 2: mod-343 sigma-solve with linear (rho2,omega2) shift ----
live343 := {};        // <R mod 343, w mod 343>
smooth343 := {};      // subset with depth-2 Hensel-smooth branch
branchCount := 0;
for s49 in surv49 do
    branchCount +:= 1;
    rho1 := s49[1]; om1 := s49[2]; aux0 := aux7[s49[3]];
    Rv := 1 + 7*rho1; Wv := 7*om1;
    base := [Rv, Wv] cat aux0;
    c0 := [Fp!(((Z!Evaluate(E416[i], base)) div 49) mod p) : i in [1..5]];
    dR := [Fp!((Z!Evaluate(DerR[i], base)) mod p) : i in [1..5]];
    dw := [Fp!((Z!Evaluate(DerW[i], base)) mod p) : i in [1..5]];
    Gp := Matrix(Fp, 5, 5,
        [Fp!(((Z!Evaluate(DerAux[i][j], base)) div p) mod p) : i in [1..5], j in [1..5]]);
    // scan sigma; for each, solve dR*rho2 + dw*omega2 = -(Qvals + G'sigma + c0)
    for si in [1..#sigmaList] do
        sv := sigmaList[si];
        svf := Vector(Fp, [Fp!s : s in sv]);
        Gs := [ &+[Gp[i][j]*svf[j] : j in [1..5]] : i in [1..5] ];
        vvec := [ Qvals[si][i] + Gs[i] + c0[i] : i in [1..5] ];
        // 5 equations dR_i*r2 + dw_i*o2 = -vvec_i over F_7
        Mm := Matrix(Fp, 5, 2, &cat[[dR[i], dw[i]] : i in [1..5]]);
        rhs := Vector(Fp, [-vvec[i] : i in [1..5]]);
        cons, part := IsConsistent(Transpose(Mm), rhs);
        if not cons then continue; end if;
        kerM := KernelMatrix(Transpose(Mm));   // solutions (r2,o2) coset
        dk := Nrows(kerM);
        sols2 := [];
        for cnt2 in [0..p^dk-1] do
            x0 := cnt2; cf := [];
            for _ in [1..dk] do Append(~cf, x0 mod p); x0 := x0 div p; end for;
            pt := part;
            for r in [1..dk] do pt := pt + (Fp!cf[r])*kerM[r]; end for;
            Append(~sols2, pt);
        end for;
        // smoothness at this sigma: M = Grad(aux0+7*sigma)/7 mod 7
        aux1 := [aux0[j] + 7*sv[j] : j in [1..5]];
        base1 := [Rv, Wv] cat aux1;
        Msm := Matrix(Fp, 5, 5,
            [Fp!(((Z!Evaluate(DerAux[i][j], base1)) div p) mod p) : i in [1..5], j in [1..5]]);
        issmooth := Rank(Msm) eq 5;
        for pt in sols2 do
            r2 := Z!pt[1]; o2 := Z!pt[2];
            key := <(Rv + 49*r2) mod 343, (Wv + 49*o2) mod 343>;
            Include(~live343, key);
            if issmooth then Include(~smooth343, key); end if;
        end for;
    end for;
    if branchCount mod 25 eq 0 then
        printf "  ... branch %o/%o  live343=%o smooth=%o\n",
            branchCount, #surv49, #live343, #smooth343;
    end if;
end for;

printf "STAGE2 mod-343: live (R,w) residues = %o / 2401, depth-2 Hensel-smooth = %o\n",
    #live343, #smooth343;
if #live343 eq 0 then
    print "OBSTRUCTION: stratum <1,0> dies at mod 343.";
elif #smooth343 gt 0 then
    print "LIVE WITH Q_7 POINTS: depth-2 smooth branches exist; <1,0> carries genuine 7-adic [4,16] points.";
    // print a compact digest of smooth residues
    sm := Sort(Setseq(smooth343));
    printf "first smooth residues (R,w) mod 343: %o\n",
        [sm[i] : i in [1..Min(12,#sm)]];
else
    print "live but no depth-2 smooth branch found; deeper analysis needed.";
end if;
print "DONE";
quit;
