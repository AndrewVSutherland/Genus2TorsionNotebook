//////////////////////////////////////////////////////////////////////
//  [2,24] W-split contact cover: GENUS of curve slices.
//
//  Follow-up to notes/agent_a2_24_wsplit_contact_cover.md "Next exact move"
//  (#1): eliminate M on the strongest fixed slices and compute their actual
//  curve genera, to turn the point-poor evidence into geometric (Faltings)
//  evidence.  genus >= 2  =>  finitely many rational points on that slice.
//
//  Reuses TriangularSystem VERBATIM from
//  code/agent_a2_24_wsplit_contact_cover.m (Delta-fixed triangular reduction
//  removing Kappa,H2,H1,H0 from f - H^2 - Kappa*(x^2+U*x+V)^3, V=(U^2-D)/4).
//
//  Slices computed (from the note's beta-zero analysis, Delta=2, U=0, beta=0,
//  where W ~ (t-1)^2*(p-r)^2 splits into t=1 and p=r branches):
//     t=1  branch: note reports a degree-16 plane factor (squared)
//     p=r  branch: note reports a degree-18 plane factor (squared)
//
//  Usage: magma -b Delta:=2 code/agent_a2_24_wsplit_slice_genus.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals();
if not assigned Delta then Delta := 2; elif Type(Delta) eq MonStgElt then Delta := StringToInteger(Delta); end if;

// ---- VERBATIM from agent_a2_24_wsplit_contact_cover.m ----
function TriangularSystem(R, DeltaVal)
    r := R.1; p := R.2; t := R.3; beta := R.4; U := R.5; M := R.6;
    F := FieldOfFractions(R);
    P<X> := PolynomialRing(F);

    rF := F!r; pF := F!p; tF := F!t; bF := F!beta;
    UF := F!U; MF := F!M; DeltaF := F!DeltaVal;
    VF := (UF^2 - DeltaF)/4;

    e := tF^2 - 2*pF*tF/rF;
    d := e + 2*pF - rF^2;
    lambda := rF/tF;
    a := rF^2 - lambda;
    b := 2*rF*pF - 2*lambda*(pF + rF*tF) + 2*rF*lambda;
    c := pF^2 + 2*pF*rF^2 - rF^4 - rF^3*tF - rF*pF^2/tF
         - lambda*(rF^2 + e)
         + 2*lambda*(rF*pF + rF^2*tF - 3*pF*tF + rF*tF^2);
    Qpoly := X^2 + d;
    q := a*X^2 + b*X + c;
    f := q*(Qpoly^2 + q);
    fc := [Coefficient(f, i) : i in [0..6]];

    K := fc[7] - MF^2;
    H2 := (fc[6] - 3*K*UF)/(2*MF);
    H1 := (fc[5] - H2^2 - 3*K*(UF^2 + VF))/(2*MF);
    H0 := (fc[4] - 2*H2*H1 - K*(UF^3 + 6*UF*VF))/(2*MF);

    R2 := fc[3] - (H1^2 + 2*H2*H0) - K*(3*UF^2*VF + 3*VF^2);
    R1 := fc[2] - 2*H1*H0 - K*(3*UF*VF^2);
    R0 := fc[1] - H0^2 - K*VF^3;
    Wexpr := bF^2*(2*d + a + bF^2)^2 - 4*(d^2 + c)*bF^2 - b^2;

    return [R!Numerator(R2), R!Numerator(R1), R!Numerator(R0)],
           R!Numerator(Wexpr);
end function;
// ---- end verbatim ----

R<r,p,t,beta,U,M> := PolynomialRing(Q, 6);
resids, Wsplit := TriangularSystem(R, Delta);
printf "Delta=%o  resid degs: c2=%o c1=%o c0=%o ; wsplit deg=%o\n",
    Delta, TotalDegree(resids[1]), TotalDegree(resids[2]), TotalDegree(resids[3]),
    TotalDegree(Wsplit);

// confirm the beta=0 W-split collapse W ~ (t-1)^2*(p-r)^2 (up to r,t factors)
W0 := Evaluate(Wsplit, [r,p,t,0,U,M]);
printf "beta=0: wsplit factorization = %o\n", Factorization(W0);

// 2-var target ring for plane models
P2<a2,b2> := PolynomialRing(Q, 2);

// eliminate M from a pair of (r,p)-or-(r,t) branch residuals, return plane poly
procedure DoBranch(branchName, subvals, keep1slot, keep2slot)
    // subvals: sequence of length 6 giving the substitution (M kept as R.6)
    c2 := Evaluate(resids[1], subvals);
    c1 := Evaluate(resids[2], subvals);
    c0 := Evaluate(resids[3], subvals);
    printf "\n== branch %o ==\n", branchName;
    printf "  after substitution: deg c2=%o c1=%o c0=%o  (c1 zero? %o)\n",
        TotalDegree(c2), TotalDegree(c1), TotalDegree(c0), c1 eq 0;
    // strip the shared (M-dependent) boundary/spurious factor before eliminating
    g0 := GCD(c2, c0);
    printf "  GCD(c2,c0) total degree %o : %o\n", TotalDegree(g0),
        TotalDegree(g0) le 4 select g0 else "(suppressed)";
    if TotalDegree(g0) ge 1 then c2 := c2 div g0; c0 := c0 div g0; end if;
    printf "  after strip: deg c2=%o c0=%o\n", TotalDegree(c2), TotalDegree(c0);
    // eliminate M (=R.6)
    E := Resultant(c2, c0, M);
    if E eq 0 then printf "  resultant STILL vanished; deeper shared factor -- skip\n"; return; end if;
    // map to the plane ring using the two kept slots
    imgs := [P2 | 0,0,0,0,0,0];
    imgs[keep1slot] := a2; imgs[keep2slot] := b2;
    h := hom< R -> P2 | imgs >;
    plane := h(E);
    printf "  eliminant(M): total degree %o\n", TotalDegree(plane);
    fac := Factorization(plane);
    printf "  factors: %o\n", [<TotalDegree(g[1]), g[2]> : g in fac];
    for g in fac do
        G := g[1];
        if TotalDegree(G) lt 3 then continue; end if;   // skip boundary monomials/lines
        printf "  --- factor deg %o (mult %o): computing genus ---\n",
            TotalDegree(G), g[2];
        try
            C := Curve(AffineSpace(P2), G);
            gg := Genus(C);
            printf "  >>> GENUS = %o\n", gg;
        catch ee
            printf "  genus computation error: %o\n", ee`Object;
        end try;
    end for;
end procedure;

// t=1 branch: beta=0, U=0, t=1  -> curve in (r,p);   keep slots r=1, p=2
DoBranch("t=1 (Delta=" cat IntegerToString(Delta) cat ",U=0,beta=0)",
         [r, p, R!1, R!0, R!0, M], 1, 2);

// p=r branch: beta=0, U=0, p=r  -> curve in (r,t);   keep slots r=1, t=3
DoBranch("p=r (Delta=" cat IntegerToString(Delta) cat ",U=0,beta=0)",
         [r, r, t, R!0, R!0, M], 1, 3);

print "DONE";
quit;
