
//////////////////////////////////////////////////////////////////////
//  Depth-k smoothness scan for the [4,16] second-halving cover over
//  the p=7 boundary strata of M_1(8,4) with PARTIALLY degenerate
//  aux-Jacobian (rank 1..4), e.g. <0,0> and <5,0> (rank 3).
//
//  Method (implicit-function reduction + tower):
//  At a mod-7 aux solution aux0 with r = rank(J_aux mod 7) in 1..4,
//  choose r pivot equations / r pivot aux-variables with invertible
//  r x r minor mod 7.  By Hensel, the pivot variables are unique
//  Z_7-analytic functions of the remaining data.  Substituting leaves
//  a reduced system of (5-r) equations
//      G(R, w, u_1..u_m) = 0,   m = 5-r free aux variables,
//  evaluated numerically: Newton-lift the pivots mod 7^k, then evaluate
//  the residual equations.
//
//  Tower: nodes at level k are (R,w,u) mod 7^k with G == 0 mod 7^k.
//  At each node compute c = G/7^k mod 7 and the directional derivatives
//  D_j = (G(node + 7^k e_j) - G(node))/7^k mod 7 (exact for polynomials
//  composed with 7-adic implicit functions).  Then:
//    * if rank(D) = 5-r  ==>  SMOOTH: the node carries genuine Q_7
//      points of the [4,16] cover (implicit function theorem);
//    * else children are node + 7^k * delta for delta solving
//      D*delta = -c over F_7.
//
//  Consistency check: rank(J_aux) = r makes the free-column Schur
//  complement vanish mod 7, so at level 1 the u-columns of D are 0 and
//  the (R,w)-columns reproduce the level-1 cokernel conditions of
//  agent_m18_416_p7_blowup.m.
//
//  Usage:
//    magma -b strata:="0,0;5,0" maxlevel:=4 nodecap:=40000 \
//        agent_m18_416_p7_smooth_scan.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();
p := 7;

if not assigned strata then strata := "0,0;5,0"; end if;
strataList := [];
for s in Split(strata, ";") do
    xy := Split(s, ",");
    Append(~strataList, [StringToInteger(xy[1]), StringToInteger(xy[2])]);
end for;
if not assigned maxlevel then maxlevel := 4;
elif Type(maxlevel) eq MonStgElt then maxlevel := StringToInteger(maxlevel); end if;
if not assigned nodecap then nodecap := 40000;
elif Type(nodecap) eq MonStgElt then nodecap := StringToInteger(nodecap); end if;

// ---- cleared integer E416 ----
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
DerAux := [[Derivative(E416[i], 2+j) : j in [1..5]] : i in [1..5]];
Fp := GF(p);

function EvZ(pol, v) return Z!Evaluate(pol, v); end function;

function AuxMod7(R0, w0)
    sols := [];
    for aa in [0..6] do for bb in [0..6] do for cv in [0..6] do
      for dv in [0..6] do for ev in [0..6] do
        v := [R0,w0,aa,bb,cv,dv,ev]; ok := true;
        for i in [1..5] do
            if EvZ(E416[i], v) mod p ne 0 then ok := false; break; end if;
        end for;
        if ok then Append(~sols, [aa,bb,cv,dv,ev]); end if;
    end for; end for; end for; end for; end for;
    return sols;
end function;

// pivot selection: r-subsets of rows/cols with invertible minor mod 7
function PickPivots(Jm, r)
    rows := [1..5]; cols := [1..5];
    for rs in Subsets({1..5}, r) do
        for cs in Subsets({1..5}, r) do
            rl := Sort(Setseq(rs)); cl := Sort(Setseq(cs));
            M := Matrix(Fp, r, r, [Jm[i][j] : j in cl, i in rl]);
            if Rank(M) eq r then return true, rl, cl; end if;
        end for;
    end for;
    return false, [], [];
end function;

// Newton-lift pivot variables to precision 7^kk.
// vals layout: [Rv, Wv, a,b,c,d,e]; free/pivot positions given by index sets.
function LiftPivots(Rv, Wv, freeIdx, freeVals, pivIdx, pivSeed, pivRows, kk)
    pk := 7^kk;
    Rk := Integers(pk);
    piv := [pv mod pk : pv in pivSeed];
    r := #pivIdx;
    for iter in [1..kk+1] do
        vals := [Rv, Wv, 0,0,0,0,0];
        for j in [1..#freeIdx] do vals[2+freeIdx[j]] := freeVals[j]; end for;
        for j in [1..r] do vals[2+pivIdx[j]] := piv[j]; end for;
        Ep := [EvZ(E416[i], vals) mod pk : i in pivRows];
        if &and[ev eq 0 : ev in Ep] then return true, piv; end if;
        M := Matrix(Rk, r, r,
            [Rk!(EvZ(DerAux[i][pivIdx[j]], vals)) : j in [1..r], i in pivRows]);
        if not IsInvertible(M) then return false, piv; end if;
        corr := Vector(Rk, [Rk!ev : ev in Ep])*Transpose(M^-1);
        piv := [ (piv[j] - Z!corr[j]) mod pk : j in [1..r] ];
    end for;
    // final check
    vals := [Rv, Wv, 0,0,0,0,0];
    for j in [1..#freeIdx] do vals[2+freeIdx[j]] := freeVals[j]; end for;
    for j in [1..r] do vals[2+pivIdx[j]] := piv[j]; end for;
    ok := &and[ EvZ(E416[i], vals) mod pk eq 0 : i in pivRows ];
    return ok, piv;
end function;

// evaluate residual equations G at (Rv,Wv,freeVals) to precision 7^kk
function GEval(Rv, Wv, freeIdx, freeVals, pivIdx, pivSeed, pivRows, residRows, kk)
    ok, piv := LiftPivots(Rv, Wv, freeIdx, freeVals, pivIdx, pivSeed, pivRows, kk);
    if not ok then return false, [], piv; end if;
    pk := 7^kk;
    vals := [Rv, Wv, 0,0,0,0,0];
    for j in [1..#freeIdx] do vals[2+freeIdx[j]] := freeVals[j]; end for;
    for j in [1..#pivIdx] do vals[2+pivIdx[j]] := piv[j]; end for;
    G := [EvZ(E416[i], vals) mod pk : i in residRows];
    return true, G, piv;
end function;

printf "M18_416_P7_SMOOTH_SCAN strata=%o maxlevel=%o nodecap=%o\n",
    strataList, maxlevel, nodecap;

for st in strataList do
    R0 := st[1]; w0 := st[2];
    aux7 := AuxMod7(R0, w0);
    printf "\n=== stratum <%o,%o>: aux solutions mod 7 = %o ===\n", R0, w0, #aux7;
    branchNo := 0;
    for aux0 in aux7 do
        branchNo +:= 1;
        base := [R0, w0] cat aux0;
        Jm := [[Fp!(EvZ(DerAux[i][j], base)) : j in [1..5]] : i in [1..5]];
        JM := Matrix(Fp, 5, 5, &cat Jm);
        r := Rank(JM);
        printf "branch %o: aux0=%o rank(J_aux)=%o\n", branchNo, aux0, r;
        if r eq 0 then
            print "  rank 0: handled by agent_m18_416_p7_hessian_10.m; skipping";
            continue;
        end if;
        ok, pivRows, pivCols := PickPivots(Jm, r);
        if not ok then print "  ERROR: no invertible minor"; continue; end if;
        residRows := [i : i in [1..5] | i notin pivRows];
        freeIdx := [j : j in [1..5] | j notin pivCols];
        m := #freeIdx;
        nresid := #residRows;
        printf "  pivots: rows=%o cols(aux)=%o | resid rows=%o free aux=%o\n",
            pivRows, pivCols, residRows, freeIdx;

        // level-1 node
        nodes := [ <Z!R0, Z!w0,
                    [aux0[j] : j in freeIdx],
                    [aux0[j] : j in pivCols]> ];
        smoothFound := false;
        sampled := false;
        for k in [1..maxlevel-1] do
            pk := 7^k;
            children := [];
            rankHist := AssociativeArray();
            died := 0;
            for nd in nodes do
                Rv := nd[1]; Wv := nd[2]; fv := nd[3]; pseed := nd[4];
                okB, Gb, pivB := GEval(Rv, Wv, freeIdx, fv, pivCols, pseed,
                                       pivRows, residRows, k+1);
                if not okB then died +:= 1; continue; end if;
                // invariant: G == 0 mod 7^k
                if not &and[ g mod pk eq 0 : g in Gb ] then died +:= 1; continue; end if;
                cvec := [ Fp!((Gb[i] div pk) mod p) : i in [1..nresid] ];
                // directional derivatives: directions R, w, then free aux
                ndir := 2 + m;
                Dm := ZeroMatrix(Fp, nresid, ndir);
                okdir := true;
                for j in [1..ndir] do
                    Rj := Rv; Wj := Wv; fj := fv;
                    if j eq 1 then Rj := Rv + pk;
                    elif j eq 2 then Wj := Wv + pk;
                    else fj := [ fv[l] + (l eq j-2 select pk else 0) : l in [1..m] ];
                    end if;
                    okJ, Gj, _ := GEval(Rj, Wj, freeIdx, fj, pivCols, pivB,
                                        pivRows, residRows, k+1);
                    if not okJ then okdir := false; break; end if;
                    for i in [1..nresid] do
                        Dm[i][j] := Fp!(((Gj[i] - Gb[i]) div pk) mod p);
                    end for;
                end for;
                if not okdir then died +:= 1; continue; end if;
                rk := Rank(Dm);
                key := Sprint(rk);
                if IsDefined(rankHist,key) then rankHist[key] +:= 1;
                else rankHist[key] := 1; end if;
                if rk eq nresid then
                    smoothFound := true;
                    printf "  SMOOTH node at level %o: (R,w) == (%o,%o) mod %o, free=%o\n",
                        k, Rv, Wv, pk, fv;
                    printf "  ==> genuine Q_7 [4,16] points over this residue.\n";
                    break;
                end if;
                // children: solve Dm*delta = -c
                rhs := Vector(Fp, [-cvec[i] : i in [1..nresid]]);
                cons, part := IsConsistent(Transpose(Dm), rhs);
                if not cons then died +:= 1; continue; end if;
                ker := KernelMatrix(Transpose(Dm));
                dk := Nrows(ker);
                for cnt in [0..p^dk-1] do
                    x0 := cnt; cf := [];
                    for _ in [1..dk] do Append(~cf, x0 mod p); x0 := x0 div p; end for;
                    dl := part;
                    for rr in [1..dk] do dl := dl + (Fp!cf[rr])*ker[rr]; end for;
                    Rc := Rv + pk*(Z!dl[1]);
                    Wc := Wv + pk*(Z!dl[2]);
                    fc := [ fv[l] + pk*(Z!dl[2+l]) : l in [1..m] ];
                    Append(~children, <Rc, Wc, fc, pivB>);
                end for;
            end for;
            if smoothFound then break; end if;
            // report rank histogram and level stats
            rhstr := "";
            for key in Sort([kk2 : kk2 in Keys(rankHist)]) do
                rhstr cat:= Sprintf(" rank%o:%o", key, rankHist[key]);
            end for;
            rwproj := #Seqset([ <nd[1] mod 7^(k+1), nd[2] mod 7^(k+1)> : nd in children ]);
            printf "  level %o->%o: nodes %o -> children %o (died %o)%o  (R,w)-proj=%o\n",
                k, k+1, #nodes, #children, died, rhstr, rwproj;
            if #children eq 0 then
                printf "  BRANCH DEAD at level %o.\n", k+1;
                break;
            end if;
            if #children gt nodecap then
                sampled := true;
                // deterministic subsample
                step := (#children div nodecap) + 1;
                children := [children[i] : i in [1..#children by step]];
                printf "  (capped: subsampled to %o nodes; results beyond here are SAMPLED)\n",
                    #children;
            end if;
            nodes := children;
        end for;
        if not smoothFound then
            printf "  branch %o: NO smooth node through level %o%o\n",
                branchNo, maxlevel, sampled select " (sampled)" else " (exhaustive)";
        end if;
    end for;
end for;
print "DONE";
quit;
