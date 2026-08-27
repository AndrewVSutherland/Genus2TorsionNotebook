
//////////////////////////////////////////////////////////////////////
//  [2,24] backward construction, step 2: p-adic Newton lift + rational
//  reconstruction.
//
//  Family f = (x^2-1)(x^2+ax+b)(x^2+cx+d).  On the generic CF cell
//  [3,1,1,...,1] (21 deg-1 quotients), the order-24 closure imposes
//  exactly 2 conditions: the x^2 and x^1 coefficients of Q_22 vanish.
//  Define G(a,b,c,d) = (coeff2(Q_22), coeff1(Q_22)) computed NUMERICALLY
//  over Z/p^k (no symbolic blowup).  For each F_p seed on the generic
//  cell: freeze (a,b) at integer lifts (sweep small shifts a0+p*s), and
//  Newton-solve (c,d) in Z_p via modified Newton (Jacobian mod p from
//  finite differences).  Rationally reconstruct (c,d); verify any
//  reconstruction exactly over Q; then torsion + simplicity.
//
//  Usage: magma -b p:=11 prec:=40 SweepS:=2 agent_a2_24_construct_lift.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB := 2; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
if not assigned p then p := 11; elif Type(p) eq MonStgElt then p := StringToInteger(p); end if;
if not assigned prec then prec := 40; elif Type(prec) eq MonStgElt then prec := StringToInteger(prec); end if;
if not assigned SweepS then SweepS := 2; elif Type(SweepS) eq MonStgElt then SweepS := StringToInteger(SweepS); end if;

Q := Rationals();
Z := Integers();
PQ<xq> := PolynomialRing(Q);

//------------------------------------------------------------------
// CF machinery over Z/p^k (unit-division only; generic cell)
//------------------------------------------------------------------
// returns ok, G1, G2: the x^2 and x^1 coefficients of Q_22 (in Z/p^k),
// following the generic pattern; ok=false if a division fails.
function ClosureG(av, bv, cv, dv, Rk)
    Pk := PolynomialRing(Rk); y := Pk.1;
    f := (y^2 - 1)*(y^2 + av*y + bv)*(y^2 + cv*y + dv);
    // s = sqrt poly part
    s := y^3;
    for k in [1..3] do
        dd := f - s^2;
        if Degree(dd) le 2 then break; end if;
        lc2 := Rk!2;
        co := Coefficient(dd, 6-k);
        s := s + (co * (lc2)^-1) * y^(3-k);
    end for;
    Pi := Pk!0; Qi := Pk!1;
    for i in [0..21] do
        lcq := Coefficient(Qi, Degree(Qi));
        if not IsUnit(lcq) then return false, Rk!0, Rk!0; end if;
        // ai = (Pi + s) div Qi  -- manual division by unit-lc poly
        num := Pi + s;
        ai := Pk!0;
        inv := lcq^-1;
        dq := Degree(Qi);
        while Degree(num) ge dq do
            t := Coefficient(num, Degree(num))*inv;
            ai := ai + t*y^(Degree(num)-dq);
            num := num - t*y^(Degree(num)-dq)*Qi;
        end while;
        Pn := ai*Qi - Pi;
        // Qn = (f - Pn^2) div Qi (exact); same manual division
        num2 := f - Pn^2;
        Qn := Pk!0;
        while Degree(num2) ge dq do
            t := Coefficient(num2, Degree(num2))*inv;
            Qn := Qn + t*y^(Degree(num2)-dq);
            num2 := num2 - t*y^(Degree(num2)-dq)*Qi;
        end while;
        if num2 ne 0 then return false, Rk!0, Rk!0; end if;   // not exact
        Pi := Pn; Qi := Qn;
    end for;
    // after 21 quotient steps past a_0, Qi is Q_22: closure wants deg <= 0
    return true, Coefficient(Qi, 2), Coefficient(Qi, 1);
end function;

// modified Newton on (c,d) with (a,b) fixed; returns ok, c, d in Z (mod p^prec)
function NewtonCD(av, bv, c0, d0)
    // Jacobian mod p via finite differences at precision p^2
    R2 := Integers(p^2);
    ok0, g1, g2 := ClosureG(R2!av, R2!bv, R2!c0, R2!d0, R2);
    if not ok0 then return false, 0, 0; end if;
    okc, g1c, g2c := ClosureG(R2!av, R2!bv, R2!(c0+p), R2!d0, R2);
    okd, g1d, g2d := ClosureG(R2!av, R2!bv, R2!c0, R2!(d0+p), R2);
    if not (okc and okd) then return false, 0, 0; end if;
    Fp := GF(p);
    J := Matrix(Fp, 2, 2, [
        Fp!(Z!(g1c - g1) div p), Fp!(Z!(g1d - g1) div p),
        Fp!(Z!(g2c - g2) div p), Fp!(Z!(g2d - g2) div p)]);
    if Rank(J) lt 2 then return false, 0, 0; end if;
    Jinv := J^-1;
    // iterate: gain one p-digit per step
    cc := Z!c0; dd := Z!d0;
    for m in [1..prec] do
        Rm := Integers(p^(m+1));
        okm, G1, G2 := ClosureG(Rm!av, Rm!bv, Rm!cc, Rm!dd, Rm);
        if not okm then return false, 0, 0; end if;
        // G should be 0 mod p^m; correction from digit m
        e1 := Z!G1; e2 := Z!G2;
        if e1 mod p^m ne 0 or e2 mod p^m ne 0 then return false, 0, 0; end if;
        v := Vector(Fp, [Fp!(e1 div p^m), Fp!(e2 div p^m)]);
        delta := v*Transpose(Jinv);
        cc := cc - (Z!delta[1])*p^m;
        dd := dd - (Z!delta[2])*p^m;
        cc := cc mod p^(m+1); dd := dd mod p^(m+1);
    end for;
    return true, cc mod p^prec, dd mod p^prec;
end function;

// exact CF order over Q (verification)
function CFOrderQ(f, ms)
    s := xq^3;
    for k in [1..3] do
        dd := f - s^2; if Degree(dd) le 2 then break; end if;
        s := s + (Coefficient(dd,6-k)/2)*xq^(3-k);
    end for;
    Pi := PQ!0; Qi := PQ!1; tot := 0;
    for i in [0..ms] do
        if Qi eq 0 then return 0; end if;
        ai := (Pi+s) div Qi; tot +:= Degree(ai);
        Pn := ai*Qi - Pi;
        if (f-Pn^2) mod Qi ne 0 then return 0; end if;
        Qi := (f-Pn^2) div Qi; Pi := Pn;
        if i ge 1 and Degree(Qi) le 0 and Qi ne 0 then return tot; end if;
    end for;
    return 0;
end function;

//------------------------------------------------------------------
// seeds (recomputed quickly over F_p, generic cell only)
//------------------------------------------------------------------
Fpf := GF(p);
Pf<xf> := PolynomialRing(Fpf);
function CFOrderPatternFp(f, ms)
    s := xf^3;
    for k in [1..3] do
        dd := f - s^2; if Degree(dd) le 2 then break; end if;
        s := s + (Coefficient(dd,6-k)/(Fpf!2))*xf^(3-k);
    end for;
    Pi := Pf!0; Qi := Pf!1; tot := 0; gen := true;
    for i in [0..ms] do
        if Qi eq 0 then return 0, false; end if;
        ai := (Pi+s) div Qi; tot +:= Degree(ai);
        if i ge 1 and Degree(ai) ne 1 then gen := false; end if;
        Pn := ai*Qi - Pi;
        if (f-Pn^2) mod Qi ne 0 then return 0, false; end if;
        Qi := (f-Pn^2) div Qi; Pi := Pn;
        if i ge 1 and Degree(Qi) le 0 and Qi ne 0 then return tot, gen; end if;
    end for;
    return 0, false;
end function;

printf "CONSTRUCT LIFT p=%o prec=%o sweep=+-%o\n", p, prec, SweepS;
seeds := [];
for a in Fpf do for b in Fpf do for c in Fpf do for d in Fpf do
    f := (xf^2-1)*(xf^2+a*xf+b)*(xf^2+c*xf+d);
    if not IsSquarefree(f) then continue; end if;
    ord, gen := CFOrderPatternFp(f, 40);
    if ord eq 24 and gen then Append(~seeds, [Z!a, Z!b, Z!c, Z!d]); end if;
end for; end for; end for; end for;
printf "generic order-24 seeds: %o\n", #seeds;

pk := p^prec;
tried := 0; lifted := 0; recon := 0; verified := 0;
for sd in seeds do
    a0 := sd[1]; b0 := sd[2]; c0 := sd[3]; d0 := sd[4];
    // centered lifts
    ac := a0 gt p div 2 select a0 - p else a0;
    bc := b0 gt p div 2 select b0 - p else b0;
    for sa in [-SweepS..SweepS] do for sb in [-SweepS..SweepS] do
        av := ac + p*sa; bv := bc + p*sb;
        tried +:= 1;
        ok, cc, dd := NewtonCD(av, bv, c0, d0);
        if not ok then continue; end if;
        lifted +:= 1;
        okc, crat := RationalReconstruction(Integers(pk)!cc);
        okd, drat := RationalReconstruction(Integers(pk)!dd);
        if not (okc and okd) then continue; end if;
        // small-height filter
        if Max(Abs(Numerator(crat)), Denominator(crat)) gt 10^12
           or Max(Abs(Numerator(drat)), Denominator(drat)) gt 10^12 then continue; end if;
        recon +:= 1;
        fq := (xq^2-1)*(xq^2+av*xq+bv)*(xq^2+crat*xq+drat);
        if Discriminant(fq) eq 0 then continue; end if;
        ordq := CFOrderQ(fq, 40);
        printf "RECON a=%o b=%o c=%o d=%o CForderQ=%o\n", av, bv, crat, drat, ordq;
        if ordq eq 24 then
            verified +:= 1;
            printf "!!!! ORDER24_CONSTRUCTED f=%o\n", fq;
            J := Jacobian(HyperellipticCurve(fq));
            inv := Invariants(TorsionSubgroup(J));
            printf "TORSION = %o\n", inv;
        end if;
    end for; end for;
end for;
printf "DONE tried=%o lifted=%o reconstructed=%o verified24=%o\n",
    tried, lifted, recon, verified;
quit;
