// (8,8) lane, step 1.1 part A: NUMERIC lift-layer survey on Lambda_334.
// The lift condition for (8,4)/(8,8) is D in 2*J1(Q) for order-4 generators
// D of Sigma.  Via the even-degree x-T (Cassels) map on J1: for each
// quadratic factor l_j of f1, val_j := u_D(theta_j) in K_j = Q[x]/(l_j)
// must become a square after ONE COMMON rational twist lambda.
//
// Closed form used here (new to this pass): writing val_j = p + q*theta with
// theta^2 = e_j (after completing the square), a rational lambda with
// lambda*val_j in (K_j*)^2 exists iff w^2 := p^2 - e_j*q^2 = Norm(val_j) is a
// rational square (the norm presieve N_j), and then the two solutions are
//   lambda_j^{+-}  ==  2*(p +- w)   mod  Q*^2
// (their product is e_j mod squares, consistent with the <disc> ambiguity).
// So: D in 2*J1(Q)  ==>  the three 2-element sets {lambda_j^+, lambda_j^-}
// have a common square class.  This script computes, at the 13 recorded
// ground-truth members (J1 = [4,4] exactly at ALL of them), for each of the
// 12 order-4 torsion elements:
//   - the 2-torsion class it doubles to ([l1]/[l2]/[l3]),
//   - the norm presieve pattern (must reproduce liftdiag: N3 always square,
//     exactly 4 of 12 elements pass 111),
//   - for the 111-passers, the lambda sets and their intersection.
// Since J1 = [4,4] exactly, the intersection must be EMPTY for every element
// (a nonempty one at these members would be a LIFT_ANOMALY: either the
// even-degree kernel caveat of the x-T map or a bug -- printed loudly).
// The printed UD/LAM lines are the interpolation data for part B (symbolic).
//
// Run: nohup magma -b code/claude_prod_09_88_liftlocus.m > results/claude_prod_09_88_liftlocus_partA.log 2>&1 &

SetColumns(0);
SetSeed(1);
SetMemoryLimit(3*10^9);

load "code/claude_prod_09_88_defs.m";
x := P88.1;

// lambda square-class sets for val in K = Q[x]/(x^2 + beta*x + gamma).
// Returns: ok (norm presieve), and the set of squarefree lambda reps.
function LambdaSet(u_D, lj)
    lj := lj / LeadingCoefficient(lj);
    beta := Coefficient(lj, 1); gamma := Coefficient(lj, 0);
    e := beta^2/4 - gamma;                     // theta^2 = e after shift
    if IsSquare(e) then
        // split algebra: K = Q x Q, componentwise conditions
        rr := Roots(lj);
        if #rr lt 2 then return false, {}; end if;  // double root: degenerate
        v1 := Evaluate(u_D, rr[1][1]); v2 := Evaluate(u_D, rr[2][1]);
        if v1 eq 0 or v2 eq 0 then return false, {}; end if;
        if SqfPart(v1) ne SqfPart(v2) then return false, {}; end if;
        return true, { SqfPart(v1) };
    end if;
    K<th> := NumberField(P88!(x^2 - e));
    // val at the root theta0 = th - beta/2 of lj
    val := Evaluate(u_D, th - beta/2);
    if val eq 0 then return false, {}; end if;
    es := Eltseq(val);                          // val = p + q*th
    p := es[1]; q := es[2];
    nrm := p^2 - e*q^2;
    okn, w := IsSquare(nrm);
    if not okn then return false, {}; end if;   // norm presieve fails
    if q eq 0 then
        return true, { SqfPart(p), SqfPart(p*e) };
    end if;
    lams := {};
    for ww in [w, -w] do
        lam := 2*(p + ww);
        if lam ne 0 then Include(~lams, SqfPart(lam)); end if;
    end for;
    return true, lams;
end function;

members := [
    // [m, n, v] stage-1 rows (J1 = [4,4] exact, recorded)
    [2, 1, 1], [2, 1, -1], [2, 1, 2], [2, 3, 1], [3, 1, 1], [3, 2, -1],
    [1/2, 2, 1], [2/3, 1, 1], [2, 5, 1/2], [3, 1/3, 1],
    // double-stage-1 rows at base (3,1/3) (J1 = [4,4] exact, recorded)
    [3, 1/3, -729/17500], [3, 1/3, 26244/7975], [3, 1/3, 729/38425]
];

anomalies := 0;
for row in members do
    m := row[1]; n := row[2]; v := row[3];
    s, t := StageOneST(m, n);
    h, g1, a, b, c := Lambda334(s, t, v);
    // the three quadratic factors of f1 (same x-coordinate as g1)
    l1 := (-a+b+c-1)*x^2 + (2*a - 2*b*c)*x + (a*b*c - a*b - a*c + b*c);
    l2 := -x^2 + b*c;
    l3 := x^2 - a;
    ljs := [l1/LeadingCoefficient(l1), l2/LeadingCoefficient(l2), l3];
    J1 := Jacobian(HyperellipticCurve(IntSextic(g1)));
    Tg, mp := TorsionSubgroup(J1);
    inv := Invariants(Tg);
    printf "MEMBER m=%o n=%o v=%o J1=%o\n", m, n, v, inv;
    assert inv eq [4,4];
    Tcls := [J1![lj, P88!0] : lj in ljs];
    n111 := 0; nN3fail := 0;
    for g in Tg do
        if Order(g) ne 4 then continue; end if;
        D := mp(g);
        uD := D[1];
        DD := 2*D;
        cls := [j : j in [1..3] | DD eq Tcls[j]];
        assert #cls eq 1;
        Ns := [Resultant(lj, uD) : lj in ljs];
        pat := [IsSquare(Nj) select 1 else 0 : Nj in Ns];
        if pat[3] eq 0 then nN3fail +:= 1; end if;
        printf "UD m=%o n=%o v=%o class=l%o pat=%o%o%o u=%o\n",
               m, n, v, cls[1], pat[1], pat[2], pat[3], uD;
        if pat ne [1,1,1] then continue; end if;
        n111 +:= 1;
        lamsets := [];
        allok := true;
        for j in [1..3] do
            okj, lj_set := LambdaSet(uD, ljs[j]);
            if not okj then allok := false; break; end if;
            Append(~lamsets, lj_set);
        end for;
        if not allok then
            printf "LAM m=%o n=%o v=%o class=l%o DEGENERATE\n", m, n, v, cls[1];
            continue;
        end if;
        common := lamsets[1] meet lamsets[2] meet lamsets[3];
        printf "LAM m=%o n=%o v=%o class=l%o lam1=%o lam2=%o lam3=%o common=%o\n",
               m, n, v, cls[1], lamsets[1], lamsets[2], lamsets[3], common;
        if not IsEmpty(common) then
            anomalies +:= 1;
            printf "LIFT_ANOMALY m=%o n=%o v=%o class=l%o common=%o (J1 is [4,4] exactly -- x-T kernel caveat or bug)\n",
                   m, n, v, cls[1], common;
        end if;
    end for;
    printf "SUMMARY m=%o n=%o v=%o n111=%o nN3fail=%o\n", m, n, v, n111, nN3fail;
end for;
printf "LIFTLOCUS_PARTA_DONE anomalies=%o\n", anomalies;
quit;
