
//////////////////////////////////////////////////////////////////////
//  Exact 2-descent (x - T) formulation of the P_R-halving condition
//  in M_1(8,4): the main-sheet analysis done right.
//
//  For the odd-degree model y^2 = f = x*A*B (deg 5, lc = c4), monicize
//  F(T) = c4^4 f(T/c4) = T * At(T) * Bt(T),
//      At = c4^2*A(T/c4),  Bt = c4*B(T/c4)   (both monic quadratics),
//  X_R = -c4*R.  By Schaefer (odd degree, ker delta = 2J(Q) exactly):
//
//      P_R in 2J(Q)  <=>  u := X_R - T  is a square in A = Q[T]/F.
//
//  Componentwise (f(-R) = Y_R^2 != 0 so u is a unit):
//    C1: u0   = -c4*R          square in Q;
//    C2: N_A  = At(X_R) = c4^2*A(-R)  square in Q  (norm condition);
//    C3: N_B  = Bt(X_R) = c4*B(-R)   square in Q  (implied by C1*C2:
//        u0*N_A*N_B = F(X_R) = c4^4*Y_R^2 = square);
//    S_A, S_B: u must be an actual square in each quadratic component
//        (norm-square is necessary, not sufficient).
//
//  This script:
//   1. derives and FACTORS the conditions C1, C2, C3 over Q(R,w);
//   2. validates the full criterion against Magma's IsDivisibleBy on
//      random rational samples, and against finite-field Jacobian
//      divisibility over F_p (where halving points exist: positive
//      controls);
//   3. runs a fast necessary-condition sieve over (R,w) (heights far
//      beyond the old exact searches), exact-checking survivors.
//
//  Usage: magma -b height:=200 agent_m18_416_descent_conditions.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

if not assigned height then height := 100;
elif Type(height) eq MonStgElt then height := StringToInteger(height); end if;
if not assigned validate then validate := true;
elif Type(validate) eq MonStgElt then validate := validate in {"true","True","1","yes"}; end if;

// ---- symbolic conditions over Q(R,w) ----
K2<R,w> := RationalFunctionField(Q, 2);
PX<x> := PolynomialRing(K2);
t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
Apol := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
Bpol := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);
c4 := R + 2 + 4*t;
f := x*Apol*Bpol;

C1 := -c4*R;
C2 := Evaluate(Apol, -R);            // A(-R); c4^2 factor is a square
C3 := c4*Evaluate(Bpol, -R);         // c4*B(-R)

print "== squareclass conditions for P_R in 2J(Q) ==";
for pair in [<"C1 = -c4*R", C1>, <"C2 = A(-R)", C2>, <"C3 = c4*B(-R)", C3>] do
    nm := pair[1]; val := pair[2];
    num := Numerator(val); den := Denominator(val);
    printf "%o:\n  numerator factors: %o\n  denominator factors: %o\n",
        nm, Factorization(num), Factorization(den);
end for;

// consistency: C1*C2*C3 should be Y_R^2 * square
Qfac := R^2 - (K2!1/2)*R*w^2 + (K2!1/2)*R - w^2;
YR := -2*R*(R-1)^2*Qfac/(w^2-1);
prod := C1*C2*C3;
ratio := prod/YR^2;
rn := Factorization(Numerator(ratio));
rd := Factorization(Denominator(ratio));
printf "C1*C2*C3 / Y_R^2 factorization (should be all even exponents):\n num %o\n den %o\n", rn, rd;

// ---- numeric machinery ----
P<xq> := PolynomialRing(Q);

function FamilyData(Rv, wv)
    tv := (2*Rv^2 + (1-wv^2)*Rv - 2*wv^2)/(4*(wv^2-1));
    c4v := Rv + 2 + 4*tv;
    Av := xq^2 + (Rv^3 + 4*Rv^2*tv + Rv - 8*Rv*tv + 4*tv)*xq + Rv^4;
    Bv := c4v*xq^2 + (Rv^2 + 4*Rv + 1 + 8*tv)*xq + (2*Rv^2 + Rv + 4*tv);
    return xq*Av*Bv, Av, Bv, c4v;
end function;

function IsSquareQ(qv)
    qv := Q!qv;
    if qv le 0 then return qv eq 0; end if;
    return IsSquare(Numerator(qv)) and IsSquare(Denominator(qv));
end function;

// full exact criterion: u = X_R - T square in Q[T]/F (component-wise)
function DescentHalvable(Rv, wv)
    fq, Av, Bv, c4v := FamilyData(Rv, wv);
    if Degree(fq) ne 5 or Discriminant(fq) eq 0 then return false, "degenerate"; end if;
    XR := -c4v*Rv;
    // C1
    if not IsSquareQ(-c4v*Rv) then return false, "C1"; end if;
    // components: At, Bt monic transforms
    Atq := c4v^2*Evaluate(Av, xq/c4v);
    Btq := c4v*Evaluate(Bv, xq/c4v);
    Atq := P![Q!co : co in Coefficients(Atq)];
    Btq := P![Q!co : co in Coefficients(Btq)];
    for gpair in [Atq, Btq] do
        gq := gpair;
        dsc := Discriminant(gq);
        if IsSquareQ(dsc) then
            // split: two rational roots e1,e2; need XR - e_i square each
            rts := Roots(gq);
            if #rts ne 2 and not (#rts eq 1 and rts[1][2] eq 2) then return false, "rootfail"; end if;
            for rt in rts do
                if not IsSquareQ(XR - rt[1]) then return false, "split-comp"; end if;
            end for;
        else
            // quadratic field component: u_g = XR - theta in K = Q(theta)
            K<th> := NumberField(gq);
            ug := XR - th;
            if not IsSquare(ug) then return false, "quad-comp"; end if;
        end if;
    end for;
    return true, "all";
end function;

// ---- validation ----
if validate then
    print "\n== validation vs exact Jacobian divisibility ==";
    // (a) rational samples: criterion must match IsDivisibleBy(PR,2)
    function ExactHalvable(Rv, wv)
        fq, Av, Bv, c4v := FamilyData(Rv, wv);
        L := 1;
        for i in [0..Degree(fq)] do L := LCM(L, Denominator(Coefficient(fq, i))); end for;
        fI := P!(L^2*fq);
        J := Jacobian(HyperellipticCurve(fI));
        Qf := Rv^2 - (Q!1/2)*Rv*wv^2 + (Q!1/2)*Rv - wv^2;
        YRv := -2*Rv*(Rv-1)^2*Qf/(wv^2-1);
        if Evaluate(fI, -Rv) ne (L*YRv)^2 then return false, false; end if;
        PR := J![xq + Rv, P!(L*YRv)];
        ok, _ := IsDivisibleBy(PR, 2);
        return true, ok;
    end function;
    nchecked := 0; nagree := 0;
    for Rv in [Q!2, Q!3, Q!-2, Q!5/2, Q!-3/2, Q!7/3] do
        for wv in [Q!2, Q!3, Q!-2, Q!1/2, Q!5/3, Q!-4/3] do
            fq, Av, Bv, c4v := FamilyData(Rv, wv);
            if Degree(fq) ne 5 or Discriminant(fq) eq 0 then continue; end if;
            okmodel, exact := ExactHalvable(Rv, wv);
            if not okmodel then continue; end if;
            crit, why := DescentHalvable(Rv, wv);
            nchecked +:= 1;
            if crit eq exact then nagree +:= 1;
            else printf "MISMATCH R=%o w=%o exact=%o criterion=%o (%o)\n",
                Rv, wv, exact, crit, why;
            end if;
        end for;
    end for;
    printf "rational validation: %o/%o agree\n", nagree, nchecked;

    // (b) positive controls over F_p: compare with finite Jacobian
    for pp in [11,13] do
        Fp := GF(pp);
        PF<xf> := PolynomialRing(Fp);
        npos := 0; nneg := 0; agree := 0; tested := 0;
        for r0 in Fp do for w0 in Fp do
            if r0 in {Fp!0,Fp!1,Fp!-1} or w0 in {Fp!0,Fp!1,Fp!-1} then continue; end if;
            tv := (2*r0^2 + (1-w0^2)*r0 - 2*w0^2)/(4*(w0^2-1));
            c4v := r0 + 2 + 4*tv;
            if c4v eq 0 then continue; end if;
            Av := xf^2 + (r0^3 + 4*r0^2*tv + r0 - 8*r0*tv + 4*tv)*xf + r0^4;
            Bv := c4v*xf^2 + (r0^2 + 4*r0 + 1 + 8*tv)*xf + (2*r0^2 + r0 + 4*tv);
            fp := xf*Av*Bv;
            if Degree(fp) ne 5 or not IsSquarefree(fp) then continue; end if;
            Qf := r0^2 - (Fp!1/2)*r0*w0^2 + (Fp!1/2)*r0 - w0^2;
            YRp := -2*r0*(r0-1)^2*Qf/(w0^2-1);
            if Evaluate(fp, -r0) ne YRp^2 then continue; end if;
            // exact divisibility on finite Jacobian
            J := Jacobian(HyperellipticCurve(fp));
            PR := J![xf + r0, PF!YRp];
            G, phi := AbelianGroup(J);
            aG := PR @@ phi;
            coords := Eltseq(aG); invs := Invariants(G);
            exact := true;
            for i in [1..#coords] do
                g := GCD(2, invs[i]);
                if coords[i] mod g ne 0 then exact := false; break; end if;
            end for;
            // descent criterion over F_p
            XR := -c4v*r0;
            crit := IsSquare(-c4v*r0);
            if crit then
                for gq0 in [c4v^2*Evaluate(Av, xf/c4v), c4v*Evaluate(Bv, xf/c4v)] do
                    gq := PF![Fp!co : co in Coefficients(gq0)];
                    rr := Roots(gq);
                    if #rr ge 1 then
                        for rt in rr do
                            if not IsSquare(XR - rt[1]) then crit := false; break; end if;
                        end for;
                    else
                        Fq2 := GF(pp^2);
                        th := Roots(PolynomialRing(Fq2)!gq)[1][1];
                        if not IsSquare(Fq2!XR - th) then crit := false; end if;
                    end if;
                    if not crit then break; end if;
                end for;
            end if;
            tested +:= 1;
            if exact then npos +:= 1; else nneg +:= 1; end if;
            if crit eq exact then agree +:= 1;
            else printf "F_%o MISMATCH r=%o w=%o exact=%o crit=%o\n", pp, r0, w0, exact, crit;
            end if;
        end for; end for;
        printf "F_%o validation: tested=%o (halvable=%o) agree=%o\n",
            pp, tested, npos, agree;
    end for;
end if;

// ---- fast sieve over (R,w) ----
printf "\n== fast descent sieve, height=%o ==\n", height;
function RationalParametersOfHeight(Bd)
    vals := [];
    for den in [1..Bd] do for num in [-Bd..Bd] do
        if GCD(num, den) eq 1 then Append(~vals, Q!num/den); end if;
    end for; end for;
    return Sort(Setseq(Seqset(vals)));
end function;
params := RationalParametersOfHeight(height);
tested := 0; c1pass := 0; c2pass := 0; fullpass := 0;
survivors := [];
for Rv in params do
    if Rv in {Q!0, Q!1, Q!-1} then continue; end if;
    for wv in params do
        if wv in {Q!0, Q!1, Q!-1} then continue; end if;
        tested +:= 1;
        // C1: -c4*R square; c4 = (R^2-1)/(w^2-1)
        c4v := 2*(Rv^2-1)/(wv^2-1);   // NB factor 2 (Factorization drops units)
        if not IsSquareQ(-c4v*Rv) then continue; end if;
        c1pass +:= 1;
        // C2: A(-R) square
        tv := (2*Rv^2 + (1-wv^2)*Rv - 2*wv^2)/(4*(wv^2-1));
        A_mR := Rv^2 - (Rv^3 + 4*Rv^2*tv + Rv - 8*Rv*tv + 4*tv)*Rv + Rv^4;
        if not IsSquareQ(A_mR) then continue; end if;
        c2pass +:= 1;
        // full second-stage
        ok, why := DescentHalvable(Rv, wv);
        if not ok then continue; end if;
        fullpass +:= 1;
        printf "DESCENT_PASS R=%o w=%o  (P_R halvable!)\n", Rv, wv;
        Append(~survivors, <Rv, wv>);
    end for;
end for;
printf "sieve: tested=%o C1pass=%o C2pass=%o FULL=%o\n",
    tested, c1pass, c2pass, fullpass;

// exact verification + torsion for survivors
for sv in survivors do
    Rv := sv[1]; wv := sv[2];
    fq, Av, Bv, c4v := FamilyData(Rv, wv);
    L := 1;
    for i in [0..Degree(fq)] do L := LCM(L, Denominator(Coefficient(fq, i))); end for;
    fI := P!(L^2*fq);
    J := Jacobian(HyperellipticCurve(fI));
    Qf := Rv^2 - (Q!1/2)*Rv*wv^2 + (Q!1/2)*Rv - wv^2;
    YRv := -2*Rv*(Rv-1)^2*Qf/(wv^2-1);
    PR := J![xq + Rv, P!(L*YRv)];
    ok, half := IsDivisibleBy(PR, 2);
    printf "SURVIVOR R=%o w=%o exact_halvable=%o", Rv, wv, ok;
    if ok then
        invs := Invariants(TorsionSubgroup(J));
        printf " torsion=%o f=%o", invs, fI;
    end if;
    printf "\n";
end for;
print "DONE";
quit;
